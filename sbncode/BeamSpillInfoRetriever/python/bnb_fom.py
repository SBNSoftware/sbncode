#!/usr/bin/env python3
"""
bnb_fom.py -- standalone re-implementation of the SBN BNB beam-quality Figure of Merit.

Port of sbncode/BeamSpillInfoRetriever/getFOM.cpp (+ the FOM selection in
sbncode/CAFMaker/FillExposure.cxx), vectorized with numpy, with no ROOT/art
dependency.  It runs on the output of dump_bnb_quality.py (.root or .parquet),
which already carries every input the FOM needs (BPMs, BPM offsets, toroids,
multiwire database widths and the raw 96-channel multiwire profiles).

Two versions are computed side by side:

  legacy  -- reproduces the production C++ as written (including its bugs), so
             FOM_legacy should close against the FOM stored in the CAF.
  fixed   -- same physics with the bugs below corrected.

Bugs / issues found in the production code (getFOM.cpp unless noted)
---------------------------------------------------------------------
 B1  Missing-device handling is dead code.  Every device is pushed into a
     one-element vector, so the `.empty()` checks (TOR875 fallback, "return 2/3
     when BPM data missing") never fire.  Missing devices arrive as -999 and are
     used as real numbers: a -999 BPM or -999 BPM *offset* moves the beam ~1 m off
     target, giving FOM ~ 0 and the spill silently lands at FOM = -999 in the CAF.
 B2  TOR860 = -999 (x1e12) is never replaced by TOR875; a negative intensity
     makes the emittance negative -> sqrt(<0) -> NaN.
 B3  The "pre-fit" FOM (database M876/M875 widths) applies the chi2 cut from the
     *last multiwire profile fit* -- unrelated to the database widths.  If no
     profile was fitted (no multiwire data, e.g. Run 2), chi2x/chi2y are
     uninitialized (undefined behavior) and the outcome is random.
 B4  Multiwire profile <-> spill matching (ICARUSBNBRetriever matchMultiWireData)
     has no maximum time difference, and silently falls back to index 0 when no
     best match exists.  Profiles from spills tens of seconds away are used.
     (M875BB_spill_time_diff is stored but never checked.)
 B5  The fit result status is never checked (TFitResultPtr dereferenced
     blindly); histogram leaked when the profile is empty.
 B6  Units: angles are slopes in mm/m (= mrad) but go through atan() as if they
     were tangents.  Small effect for typical slopes (<1% at 0.2 mrad) but wrong.
 B7  No intensity requirement: beam-off spills (TOR860 <= 0) get FOM ~ 1.
 B8  FillExposure.cxx: validity is judged by value, FOM in (0,1].  A physical
     FOM = 0 (beam completely off target) is treated as "failed" and replaced by
     the next estimator; values from B1 garbage pass through the same way.
 B9  Cosmetic: pi = 3.14159 (normalization +8e-7, FOM saturates at exactly 1 for
     losses < 1e-6), MWRtoroidDelay commented as ms but is s, "24 cm" profile
     coverage is really 24 mm, VP873 used under the name "VTG1".

Fixed-mode choices (all switchable):
  * a device equal to -999 (or non-finite) is "missing"; TOR860 -> TOR875 fallback
  * missing BPM offsets can be filled from the nearest valid value in the same run
    (offsets are slow setpoints)                         [--fill-offsets]
  * multiwire profiles only used if |spill_time_diff| < --mw-max-dt (default 33 ms)
  * pre-fit FOM uses only the width window (there is no chi2 for database widths)
  * fit status / ndf checked; slopes not passed through atan
  * output FOM_best is the best available estimator, with fom_source and a
    status bitmask saying why (see STATUS_* below).  No +100 offset trick.

Usage
-----
  ./bnb_fom.py beamqual.root -o beamqual_fom.root -j 4     # both trees, 4 workers
  ./bnb_fom.py beamqual.root -o fom.root --tree bnbinfo --no-fill-offsets
  ./bnb_fom.py beamqual.root -o test.root --max-entries 50000   # quick look
  python3 -c "import bnb_fom; help(bnb_fom.compute_fom)"   # as a library

ROOT input is processed in streaming mode: multiwire branches are decompressed
block by block and everything else in chunks, so memory stays < ~1-2 GB even for
multi-GB inputs written as a single basket. Output is a TTree per input tree with
all scalar input columns plus the new FOM columns (multiwire arrays are not copied).
"""

import argparse
import os
import sys
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)

import numpy as np

# --------------------------------------------------------------------------- constants
MISSING = -999.0
Z = dict(vp873=191.153656, hp875=202.116104, vp875=202.3193205, hptg1=204.833267,
         vptg1=204.629608, hptg2=205.240662, vptg2=205.036835, tgt=206.870895)   # m
P875X = (0.431857, 0.158077, 0.00303551)
P875Y = (0.279128, 0.337048, 0.0)
P876X = (0.166172, 0.30999, -0.00630299)
P876Y = (0.13425, 0.580862, 0.0)
SMALL_SX, LARGE_SX, SMALL_SY, LARGE_SY, MAX_CHI2 = 0.5, 10.0, 0.3, 10.0, 20.0
NWIRE = 48

# beam optics (MiniBooNE DQ_BeamLine_*.F via MicroBooNE)
BX, AX, NX, NPX = 4.68, 0.0389, 0.0958, -0.0286
BY, AY, NY, NPY = 59.12, 2.4159, 0.4577, -0.0271
GX, GY = (1 + AX * AX) / BX, (1 + AY * AY) / BY
BEGTOCNT = np.array([
    [0.65954, 0.43311, 0.00321, 0.10786, 0.00000, 1.97230],
    [0.13047, 1.60192, 0.00034, 0.00512, 0.00000, 1.96723],
    [-0.00287, -0.03677, -0.35277, -4.68056, 0.00000, 0.68525],
    [-0.00089, -0.00430, -0.17722, -5.18616, 0.00000, 0.32300],
    [-0.00104, 0.00232, -0.00001, -0.00224, 1.00000, -0.00450],
    [0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 1.00000]])
CNTTOUPS, CNTTODNS, IDENT = np.eye(6), np.eye(6), np.eye(6)
CNTTOUPS[0, 1] = CNTTOUPS[2, 3] = -0.35710
CNTTODNS[0, 1] = CNTTODNS[2, 3] = +0.35710
BEGTOUPS, BEGTODNS = CNTTOUPS @ BEGTOCNT, CNTTODNS @ BEGTOCNT
W_UP, W_CNT, W_DN = 0.6347, 0.2812, 0.0841


def _target_grid():
    """Integration points of func_intbivar, built with the same float accumulation
    (x += 0.1) as the C++ so boundary points are included identically."""
    r, d = 4.75, 0.1
    n = int(round(2 * r / d))
    xs, x = [], -r
    for _ in range(n + 1):
        xs.append(x)
        x += d
    xs = np.array(xs)
    X, Y = np.meshgrid(xs, xs, indexing="ij")
    m = X * X + Y * Y < r * r
    return X[m], Y[m], d * d


GX_PTS, GY_PTS, DXDY = _target_grid()

# status bits for the fixed computation
STATUS = {
    "NO_TOR": 1,          # no valid, positive TOR860 or TOR875 (no beam)
    "NO_HBPM": 2,         # cannot form horizontal position/angle
    "NO_VBPM": 4,         # cannot form vertical position/angle
    "OFFSET_FILLED": 8,   # a BPM offset was filled from a neighboring spill
    "LOW_TOR": 16,        # intensity below --min-tor (beam off?)
    "MW_TIME_CUT": 32,    # >=1 multiwire profile dropped for |dt| > mw-max-dt
    "NO_MW_WIDTH": 64,    # no usable width -> fell back to nominal-width FOM
}
SOURCES = {-1: "none", 0: "fit_MMBTBB", 1: "fit_M876BB", 2: "fit_M875BB",
           3: "prefit_M876", 4: "prefit_M875", 5: "no_multiwire"}


# --------------------------------------------------------------------------- FOM core
def calc_fom(horpos, horang, verpos, verang, ppp, tgtsx=None, tgtsy=None, chunk=2000, pi=3.14159):
    """Vectorised calcFOM: returns log10 of weighted fractional beam loss (fom2).
    FOM = 1 - 10**fom2.  All inputs arrays of shape (N,); ppp in protons;
    positions in mm, angles in mrad, tgtsx/y in mm (NaN/-999 -> nominal width)."""
    horpos, horang, verpos, verang, ppp = (np.asarray(a, dtype=float) for a in
                                           (horpos, horang, verpos, verang, ppp))
    n = len(horpos)
    out = np.full(n, np.nan)
    if tgtsx is None:
        tgtsx = np.full(n, MISSING)
        tgtsy = np.full(n, MISSING)
    tgtsx, tgtsy = np.asarray(tgtsx, float), np.asarray(tgtsy, float)
    for s in range(0, n, chunk):
        sl = slice(s, s + chunk)
        out[sl] = _calc_fom_chunk(horpos[sl], horang[sl], verpos[sl], verang[sl],
                                  ppp[sl], tgtsx[sl], tgtsy[sl], pi)
    return out


def _calc_fom_chunk(hp, ha, vp, va, ppp, tsx, tsy, pi=3.14159):
    n = len(hp)
    p = ppp / 1e12
    ex = 0.1775e-06 + 0.1827e-07 * p
    ey = 0.1382e-06 + 0.2608e-08 * p
    dp = 0.4485e-03 + 0.6100e-04 * p
    dp2 = dp * dp
    S = np.zeros((n, 6, 6))
    S[:, 5, 5] = dp2
    S[:, 0, 0] = ex * BX + NX * NX * dp2
    S[:, 0, 1] = S[:, 1, 0] = -ex * AX + NX * NPX * dp2
    S[:, 1, 1] = ex * GX + NPX * NPX * dp2
    S[:, 0, 5] = S[:, 5, 0] = NX * dp2
    S[:, 1, 5] = S[:, 5, 1] = NPX * dp2
    S[:, 2, 2] = ey * BY + NY * NY * dp2
    S[:, 2, 3] = S[:, 3, 2] = -ey * AY + NY * NPY * dp2
    S[:, 3, 3] = ey * GY + NPY * NPY * dp2
    S[:, 2, 5] = S[:, 5, 2] = NY * dp2
    S[:, 3, 5] = S[:, 5, 3] = NPY * dp2
    C = np.zeros((n, 6))
    C[:, 0], C[:, 1], C[:, 2], C[:, 3] = hp, ha, vp, va

    def swim(xc, xs):
        c2 = C @ xc.T
        s2 = np.einsum("ij,njk,mk->nim", xs, S, xs)
        with np.errstate(invalid="ignore"):
            sx = np.sqrt(s2[:, 0, 0]) * 1000.0
            sy = np.sqrt(s2[:, 2, 2]) * 1000.0
            rho = s2[:, 0, 2] / np.sqrt(s2[:, 0, 0] * s2[:, 2, 2])
        return c2[:, 0], c2[:, 2], sx, sy, rho

    cx, cy, sx, sy, rho = swim(CNTTOUPS, BEGTOUPS)
    use = (tsx != MISSING) & (tsy != MISSING) & np.isfinite(tsx) & np.isfinite(tsy)
    with np.errstate(invalid="ignore", divide="ignore"):
        scx = np.where(use, tsx / sx, 1.0)
        scy = np.where(use, tsy / sy, 1.0)
    fa = intbivar(cx, cy, sx * scx, sy * scy, rho, pi)
    cx, cy, sx, sy, rho = swim(IDENT, BEGTOCNT)
    fb = intbivar(cx, cy, sx * scx, sy * scy, rho, pi)
    cx, cy, sx, sy, rho = swim(CNTTODNS, BEGTODNS)
    fc = intbivar(cx, cy, sx * scx, sy * scy, rho, pi)
    fom2 = fa * W_UP + fb * W_CNT + fc * W_DN
    guard = (fa <= -10000.) | (fb <= -10000.) | (fc <= -10000.)
    return np.where(guard, -10000., fom2)


def intbivar(cx, cy, sx, sy, rho, pi=3.14159):
    """Fraction of a bivariate Gaussian *outside* the r=4.75 mm target, as log10.
    Returns -10000 when the summed integral >= 1 (same guard as the C++)."""
    with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
        tx = (GX_PTS[None, :] - cx[:, None]) / sx[:, None]
        ty = (GY_PTS[None, :] - cy[:, None]) / sy[:, None]
        r = rho[:, None]
        z = tx * tx - 2.0 * r * tx * ty + ty * ty
        tot = np.exp(-z / (2.0 * (1.0 - r * r))).sum(axis=1)
        tot = tot * DXDY / (2.0 * pi * sx * sy * np.sqrt(1.0 - rho * rho))
        res = np.log10(1.0 - tot)
    return np.where(tot >= 1.0, -10000., res)


# --------------------------------------------------------------------------- multiwire fit
_XC = -12.0 + 0.5 * np.arange(NWIRE) + 0.25        # bin centers, mm


def process_profiles(mw, n_iter=200):
    """Vectorised processBNBprofile on profiles of shape (N, 48).
    Reproduces the C++: invert sign, subtract min, 20% threshold window, Gaussian
    chi2 fit with constant 2%-of-range errors (ROOT includes bins whose center is in
    [lo, hi]).  Returns mean, sigma, chi2/ndf, ok (fit converged & ndf>0).
    Profiles with zero sum get the C++ sentinel 99999."""
    mw = np.asarray(mw, dtype=float)
    n = mw.shape[0]
    inv = -mw
    minx = np.minimum(inv.min(axis=1), 0.0)
    maxx = np.maximum(inv.max(axis=1), 0.0)
    y = inv - minx[:, None]
    rng = maxx - minx
    thr, err = rng * 0.2, rng * 0.02
    above = y > thr[:, None]
    anyab = above.any(axis=1)
    first = np.where(anyab, above.argmax(axis=1), -1)
    last = np.where(anyab, NWIRE - above[:, ::-1].argmax(axis=1), -1)
    lo, hi = -12 + first * 0.5, -12 + last * 0.5
    inr = (_XC[None, :] >= lo[:, None]) & (_XC[None, :] <= hi[:, None])
    nb = inr.sum(axis=1)
    sumw = y.sum(axis=1) > 0

    # initial values: moments of the in-range histogram (as ROOT's gaus init does)
    w = np.where(inr, y, 0.0)
    ws = w.sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        mu = (w * _XC).sum(axis=1) / ws
        sg = np.sqrt(np.maximum((w * (_XC - mu[:, None]) ** 2).sum(axis=1) / ws, 0.0625))
    A = np.where(inr, y, 0).max(axis=1)
    p = np.stack([A, mu, sg], axis=1)
    lam = np.full(n, 1e-3)
    e2 = np.where(err > 0, err * err, np.nan)[:, None]

    def chi2_of(p):
        g = p[:, 0:1] * np.exp(-0.5 * ((_XC - p[:, 1:2]) / p[:, 2:3]) ** 2)
        return np.where(inr, (y - g) ** 2 / e2, 0.0).sum(axis=1)

    with np.errstate(all="ignore"):
        c2 = chi2_of(p)
        act = np.where(np.isfinite(c2) & (nb > 3))[0]     # rows still iterating
        for _ in range(n_iter):                           # Levenberg-Marquardt
            if len(act) == 0:
                break
            pa, ya, ma, sa = p[act], y[act], inr[act], np.sqrt(e2[act])
            u = (_XC - pa[:, 1:2]) / pa[:, 2:3]
            ex = np.exp(-0.5 * u * u)
            g = pa[:, 0:1] * ex
            J = np.stack([ex, g * u / pa[:, 2:3], g * u * u / pa[:, 2:3]], axis=2)
            J = np.where(ma[:, :, None], J, 0.0) / sa[:, :, None]
            rr = np.where(ma, ya - g, 0.0) / sa
            Jt = J.transpose(0, 2, 1)
            JTJ = Jt @ J
            JTr = (Jt @ rr[:, :, None])[:, :, 0]
            la = lam[act]
            H = JTJ + la[:, None, None] * JTJ * np.eye(3)[None]
            bad = ~np.isfinite(H).all(axis=(1, 2)) | ~np.isfinite(JTr).all(axis=1)
            H[bad] = np.eye(3)
            JTr[bad] = 0
            # row-wise pseudo-inverse: deterministic, independent of batch composition
            step = (np.linalg.pinv(H) @ JTr[:, :, None])[:, :, 0]
            pn = pa + step
            gn = pn[:, 0:1] * np.exp(-0.5 * ((_XC - pn[:, 1:2]) / pn[:, 2:3]) ** 2)
            c2n = np.where(ma, (ya - gn) ** 2 / e2[act], 0.0).sum(axis=1)
            c2a = c2[act]
            better = c2n < c2a
            p[act] = np.where(better[:, None], pn, pa)
            c2[act] = np.where(better, c2n, c2a)
            lam[act] = np.where(better, la * 0.3, la * 10.0)
            done = (better & (c2a - c2n < 1e-12 * c2a)) | (lam[act] > 1e10) | bad
            act = act[~done]
        ndf = nb - 3
        chi2ndf = c2 / ndf
    ok = sumw & (ndf > 0) & np.isfinite(chi2ndf) & np.isfinite(p).all(axis=1)
    x = np.where(sumw, p[:, 1], 99999.)
    sx = np.where(sumw, p[:, 2], 99999.)
    chi2ndf = np.where(sumw, chi2ndf, 99999.)
    return x, sx, chi2ndf, ok


def _poly(p, s):
    return p[0] + p[1] * s + p[2] * s * s


def _in_window(sx, sy):
    return (sx > SMALL_SX) & (sx < LARGE_SX) & (sy > SMALL_SY) & (sy < LARGE_SY)


# --------------------------------------------------------------------------- driver
def _profiles(arr, n):
    """per-spill multiwire arrays (awkward or list) -> (N,96) float with NaN rows when
    the profile is absent, plus presence mask."""
    import awkward as ak
    a = arr if isinstance(arr, ak.Array) else ak.Array([np.asarray(x).tolist() for x in arr])
    out = np.full((n, 2 * NWIRE), np.nan)
    if len(a) == 0:
        return out, np.zeros(n, bool)
    has = ak.to_numpy(ak.num(a)) == 2 * NWIRE
    if has.any():
        out[has] = ak.to_numpy(ak.flatten(a[has])).reshape(-1, 2 * NWIRE)
    return out, has


def _fit_profiles(P, has):
    """fit H (wires 0-47) and V (48-95) for spills with a profile."""
    n = len(has)
    res = {k: np.full(n, np.nan) for k in ("x", "sx", "chi2x", "y", "sy", "chi2y")}
    okx = np.zeros(n, bool)
    oky = np.zeros(n, bool)
    if has.any():
        idx = np.where(has)[0]
        xx, sx, cx, ox = process_profiles(P[idx, :NWIRE])
        yy, sy, cy, oy = process_profiles(P[idx, NWIRE:])
        for k, v in zip(("x", "sx", "chi2x", "y", "sy", "chi2y"), (xx, sx, cx, yy, sy, cy)):
            res[k][idx] = v
        okx[idx], oky[idx] = ox, oy
    return res, okx & oky


def _valid(a):
    a = np.asarray(a, float)
    return np.isfinite(a) & ~np.isclose(a, MISSING) & ~np.isclose(a, MISSING * 1e12)


def fill_offsets(df):
    """Fill -999 BPM offsets with the nearest-in-time valid value from the same run.
    Returns a copy and a boolean mask of spills where something was filled."""
    import pandas as pd
    df = df.copy()
    filled = np.zeros(len(df), bool)
    t = df["spill_time_sec"].to_numpy(float) + 1e-9 * df["spill_time_nsec"].to_numpy(float)
    for c in ["HP875Offset", "VP875Offset", "VP873Offset", "HPTG1Offset",
              "HPTG2Offset", "VPTG1Offset", "VPTG2Offset"]:
        v = df[c].to_numpy(float).copy()
        good = _valid(v)
        for run in np.unique(df["run"].to_numpy()):
            m = df["run"].to_numpy() == run
            gi = np.where(m & good)[0]
            bi = np.where(m & ~good)[0]
            if len(gi) == 0 or len(bi) == 0:
                continue
            order = np.argsort(t[gi])
            tg = t[gi][order]
            pos = np.clip(np.searchsorted(tg, t[bi]), 1, len(tg) - 1) if len(tg) > 1 else np.zeros(len(bi), int)
            if len(tg) > 1:
                left = np.abs(t[bi] - tg[pos - 1]) <= np.abs(t[bi] - tg[pos])
                pos = np.where(left, pos - 1, pos)
            v[bi] = v[gi][order][pos]
            filled[bi] = True
        df[c] = v
    return df, filled


def compute_fom(df, mw, mode="fixed", mw_max_dt=0.0333, min_tor=0.1e12, offsets_filled=None,
                fit_cache=None):
    """Compute the FOM for every row of `df` (one row per spill, SRBNBInfo column names).
    mw: dict {'M875BB','M876BB','MMBTBB'} -> sequence of per-spill arrays (len 0 or 96).
    fit_cache: optional dict; pass the same one to the legacy and fixed calls on the
    same spills so each multiwire profile is fitted only once.
    mode='legacy' reproduces getFOM.cpp + FillExposure.cxx; 'fixed' applies the fixes.
    Returns dict of numpy columns."""
    n = len(df)
    g = {c: df[c].to_numpy(float) for c in df.columns
         if c not in ("source_file",) and df[c].dtype.kind in "fiub"}
    legacy = mode == "legacy"
    status = np.zeros(n, np.int32)
    if offsets_filled is not None:
        status |= np.where(offsets_filled, STATUS["OFFSET_FILLED"], 0).astype(np.int32)

    # ---- intensity
    if legacy:
        tor = g["TOR860"]                                            # B2: no fallback
    else:
        t860, t875 = g["TOR860"], g["TOR875"]
        ok860, ok875 = _valid(t860) & (t860 > 0), _valid(t875) & (t875 > 0)   # <= 0: no beam
        tor = np.where(ok860, t860, np.where(ok875, t875, np.nan))
        status |= np.where(~np.isfinite(tor), STATUS["NO_TOR"], 0).astype(np.int32)
        status |= np.where(np.isfinite(tor) & (tor < min_tor), STATUS["LOW_TOR"], 0).astype(np.int32)

    # ---- position & angle at target center (straight line through two BPMs)
    def line(d0, z0, d1, z1):
        ang = (d1 - d0) / (z1 - z0)
        return ang, d0 + ang * (Z["tgt"] - z0)

    dh875 = g["HP875"] - g["HP875Offset"]
    dv875 = g["VP875"] - g["VP875Offset"]
    if legacy:                                                       # useHTG=useVTG=1
        hang, hpos = line(dh875, Z["hp875"], g["HPTG1"] - g["HPTG1Offset"], Z["hptg1"])
        vang, vpos = line(dv875, Z["vp875"], g["VP873"] - g["VP873Offset"], Z["vp873"])
        hang, vang = np.arctan(hang), np.arctan(vang)                # B6 kept
    else:
        okh875 = _valid(g["HP875"]) & _valid(g["HP875Offset"])
        okh1 = _valid(g["HPTG1"]) & _valid(g["HPTG1Offset"])
        okh2 = _valid(g["HPTG2"]) & _valid(g["HPTG2Offset"])
        a1, p1 = line(dh875, Z["hp875"], g["HPTG1"] - g["HPTG1Offset"], Z["hptg1"])
        a2, p2 = line(dh875, Z["hp875"], g["HPTG2"] - g["HPTG2Offset"], Z["hptg2"])
        hang, hpos = np.where(okh1, a1, a2), np.where(okh1, p1, p2)
        okh = okh875 & (okh1 | okh2)
        okv875 = _valid(g["VP875"]) & _valid(g["VP875Offset"])
        okv3 = _valid(g["VP873"]) & _valid(g["VP873Offset"])
        okv2 = _valid(g["VPTG2"]) & _valid(g["VPTG2Offset"])
        a1, p1 = line(dv875, Z["vp875"], g["VP873"] - g["VP873Offset"], Z["vp873"])
        a2, p2 = line(dv875, Z["vp875"], g["VPTG2"] - g["VPTG2Offset"], Z["vptg2"])
        vang, vpos = np.where(okv3, a1, a2), np.where(okv3, p1, p2)
        okv = okv875 & (okv3 | okv2)
        status |= np.where(~okh, STATUS["NO_HBPM"], 0).astype(np.int32)
        status |= np.where(~okv, STATUS["NO_VBPM"], 0).astype(np.int32)

    # ---- widths from multiwire profile fits
    tgtsx = np.full(n, np.nan)
    tgtsy = np.full(n, np.nan)
    src_fit = np.full(n, -1)
    last_chi2x = np.full(n, np.nan)          # legacy: chi2 of the *last* fit performed
    last_chi2y = np.full(n, np.nan)
    fitinfo = {}
    for code, dev, conv in ((0, "MMBTBB", None), (1, "M876BB", (P876X, P876Y)),
                            (2, "M875BB", (P875X, P875Y))):
        if fit_cache is not None and dev in fit_cache and len(fit_cache[dev]) == 3:
            P, has = None, fit_cache[dev][2]           # precomputed (streaming mode)
        else:
            P, has = _profiles(mw[dev], n)
        if not legacy:
            dt = g.get(dev + "_spill_time_diff", np.zeros(n))
            dt_bad = has & ~(np.abs(dt) < mw_max_dt)
            status |= np.where(dt_bad, STATUS["MW_TIME_CUT"], 0).astype(np.int32)
            has = has & ~dt_bad
        todo = has & (src_fit < 0)            # only fit if no better device succeeded
        if fit_cache is not None:
            if dev not in fit_cache:
                fit_cache[dev] = _fit_profiles(P, _profiles(mw[dev], n)[1])
            r_all, ok_all = fit_cache[dev][0], fit_cache[dev][1]
            r = {k: np.where(todo, v, np.nan) for k, v in r_all.items()}
            for k in ("chi2x", "chi2y", "sx", "sy"):
                r.setdefault(k, np.full(n, np.nan))
            ok = ok_all & todo
        else:
            r, ok = _fit_profiles(P, todo)
        fitinfo[dev] = r
        sx, sy = r["sx"], r["sy"]
        if conv is not None:
            sx, sy = _poly(conv[0], sx), _poly(conv[1], sy)
        good = todo & _in_window(sx, sy) & (r["chi2x"] < MAX_CHI2) & (r["chi2y"] < MAX_CHI2)
        if not legacy:
            good &= ok                                               # B5
        last_chi2x = np.where(todo, r["chi2x"], last_chi2x)
        last_chi2y = np.where(todo, r["chi2y"], last_chi2y)
        tgtsx, tgtsy = np.where(good, sx, tgtsx), np.where(good, sy, tgtsy)
        src_fit = np.where(good, code, src_fit)

    PI = 3.14159 if legacy else np.pi                                # B9
    fom_fit = np.full(n, -999.)
    m = src_fit >= 0
    if m.any():
        fom_fit[m] = 1 - 10 ** calc_fom(hpos[m], hang[m], vpos[m], vang[m], tor[m], tgtsx[m], tgtsy[m], pi=PI)

    # ---- "pre-fit" widths from the database (M876/M875 HS, VS)
    psx = np.full(n, np.nan)
    psy = np.full(n, np.nan)
    src_pre = np.full(n, -1)
    for code, dev, conv in ((3, "M876", (P876X, P876Y)), (4, "M875", (P875X, P875Y))):
        hs, vs = g[dev + "HS"], g[dev + "VS"]
        sx, sy = _poly(conv[0], hs), _poly(conv[1], vs)
        good = _in_window(sx, sy) & (src_pre < 0)
        if legacy:
            # B3: chi2 of the last profile fit; uninitialized if no fit ran -> assume pass
            c2ok = np.where(np.isfinite(last_chi2x), (last_chi2x < MAX_CHI2) & (last_chi2y < MAX_CHI2), True)
            good &= c2ok
        else:
            good &= _valid(hs) & _valid(vs)
        psx, psy = np.where(good, sx, psx), np.where(good, sy, psy)
        src_pre = np.where(good, code, src_pre)
    fom_pre = np.full(n, -999.)
    m = src_pre >= 0
    if m.any():
        fom_pre[m] = 1 - 10 ** calc_fom(hpos[m], hang[m], vpos[m], vang[m], tor[m], psx[m], psy[m], pi=PI)

    # ---- nominal-width FOM
    fom_nomw = 1 - 10 ** calc_fom(hpos, hang, vpos, vang, tor, pi=PI)

    # ---- combine
    if legacy:                                                       # FillExposure.cxx
        f32 = lambda a: a.astype(np.float32).astype(float)
        a, b, c = f32(fom_fit), f32(fom_pre), f32(fom_nomw)
        best = np.where((a > 0) & (a <= 1), a,
               np.where((b > 0) & (b <= 1), b,
               np.where((c > 0) & (c <= 1), 100. + c, -999.)))
        source = np.where((a > 0) & (a <= 1), src_fit,
                 np.where((b > 0) & (b <= 1), src_pre,
                 np.where((c > 0) & (c <= 1), 5, -1)))
        return {"FOM_legacy": best, "FOM_legacy_fit": fom_fit, "FOM_legacy_prefit": fom_pre,
                "FOM_legacy_nomw": fom_nomw, "fom_legacy_source": source.astype(np.int8)}

    fatal = (status & (STATUS["NO_TOR"] | STATUS["NO_HBPM"] | STATUS["NO_VBPM"])) != 0
    fom_fit[fatal] = np.nan
    fom_pre[fatal] = np.nan
    fom_nomw[fatal] = np.nan
    fom_fit[fom_fit == -999.] = np.nan
    fom_pre[fom_pre == -999.] = np.nan
    have_fit, have_pre = src_fit >= 0, src_pre >= 0
    best = np.where(have_fit, fom_fit, np.where(have_pre, fom_pre, fom_nomw))
    source = np.where(have_fit, src_fit, np.where(have_pre, src_pre, 5))
    source = np.where(fatal, -1, source)
    best = np.where(fatal, np.nan, best)
    status |= np.where(~fatal & ~have_fit & ~have_pre, STATUS["NO_MW_WIDTH"], 0).astype(np.int32)
    out = {"FOM_best": best, "FOM_fit": fom_fit, "FOM_prefit": fom_pre, "FOM_nomw": fom_nomw,
           "fom_source": source.astype(np.int8), "fom_status": status,
           "hpos_tgt": hpos, "hang_tgt": hang, "vpos_tgt": vpos, "vang_tgt": vang,
           "tgt_sigx": np.where(have_fit, tgtsx, psx), "tgt_sigy": np.where(have_fit, tgtsy, psy)}
    for dev in ("MMBTBB", "M876BB", "M875BB"):
        r = fitinfo[dev]
        for k in ("x", "sx", "chi2x", "y", "sy", "chi2y"):
            if k in r:
                out[f"{dev}_fit_{k}"] = r[k]
    return out


# --------------------------------------------------------------------------- I/O
def read_dump(path, tree):
    """Read a dump_bnb_quality.py output tree -> (DataFrame, {dev: awkward})."""
    import pandas as pd
    mwdevs = ("M875BB", "M876BB", "MMBTBB")
    if path.endswith((".parquet", ".pq")):
        base = path.rsplit("_spill", 1)[0].rsplit("_bnbinfo", 1)[0].rsplit(".", 1)[0]
        df = pd.read_parquet(f"{base}_{tree}.parquet")
        mw = {d: list(df.pop(d)) if d in df else [[]] * len(df) for d in mwdevs}
        return df, mw
    import uproot
    import awkward as ak
    t = uproot.open(path)[tree]
    keys = [k for k in t.keys() if not k.startswith("n")]
    arr = t.arrays(keys, library="ak")
    flat = [k for k in keys if k not in mwdevs]
    df = pd.DataFrame({k: ak.to_numpy(arr[k]) for k in flat})
    mw = {d: (arr[d] if d in keys else [[]] * len(df)) for d in mwdevs}
    return df, mw


def write_out(path, frames):
    if path.endswith((".parquet", ".pq")):
        base = path.rsplit(".", 1)[0]
        for name, df in frames.items():
            df.to_parquet(f"{base}_{name}.parquet", index=False)
        return
    import uproot
    with uproot.recreate(path) as f:
        for name, df in frames.items():
            cols = {c: df[c].to_numpy() for c in df.columns if df[c].dtype.kind in "fiub"}
            tree = f.mktree(name, {c: v.dtype for c, v in cols.items()})
            for s in range(0, len(df), 100000):
                tree.extend({c: v[s:s + 100000] for c, v in cols.items()})
        f["fom_status_bits"] = "\n".join(f"{v}:{k}" for k, v in STATUS.items())
        f["fom_sources"] = "\n".join(f"{k}:{v}" for k, v in SOURCES.items())


def summarize(name, df):
    p = lambda *a: print(*a, file=sys.stderr)
    p(f"\n=== {name}: {len(df)} spills ===")
    if "FOM" in df:
        caf = df["FOM"].to_numpy(float)
        leg = df["FOM_legacy"].to_numpy(float)
        same = np.isclose(caf, leg, atol=1e-4, rtol=0)
        p(f"closure  FOM_legacy vs CAF FOM : {same.sum()}/{len(df)} agree (|d|<1e-4)")
        bad = ~same
        if bad.any():
            cats = {"CAF -999": caf == -999, "CAF in (0,1]": (caf > 0) & (caf <= 1), "CAF >100": caf > 100}
            p("   disagreements by CAF category: " + ", ".join(f"{k}: {(bad & v).sum()}" for k, v in cats.items()))
    src = df["fom_source"].to_numpy()
    p("fixed FOM source: " + ", ".join(f"{SOURCES[k]}={np.sum(src == k)}" for k in SOURCES if np.sum(src == k)))
    st = df["fom_status"].to_numpy()
    p("status bits     : " + ", ".join(f"{k}={np.sum((st & v) != 0)}" for k, v in STATUS.items() if np.sum((st & v) != 0)))
    if "FOM" in df:
        caf = df["FOM"].to_numpy(float)
        rec = (caf == -999) & np.isfinite(df["FOM_best"].to_numpy())
        p(f"spills with CAF FOM=-999 that now have a valid FOM_best: {rec.sum()}")


# --------------------------------------------------------------------------- streaming (large files)
OFFSET_COLS = ["HP875Offset", "VP875Offset", "VP873Offset", "HPTG1Offset",
               "HPTG2Offset", "VPTG1Offset", "VPTG2Offset"]
FOM_INPUTS = ["TOR860", "TOR875", "HP875", "VP875", "VP873", "HPTG1", "HPTG2", "VPTG1", "VPTG2",
              "M875HS", "M875VS", "M876HS", "M876VS",
              "M875BB_spill_time_diff", "M876BB_spill_time_diff", "MMBTBB_spill_time_diff"] + OFFSET_COLS
ID_COLS = ["run", "subrun", "evt", "event", "spill_time_sec", "spill_time_nsec", "spill_time",
           "has_spill_info", "triggered", "subrun_pot", "file_idx", "THCURR", "FOM"]


def stream_jagged_int32(path, branch, chunk_values, counts):
    """Yield successive blocks (native int32 numpy arrays) of the flat content of a
    jagged, uproot-written TTree branch, decompressing each ROOT basket one
    compression block (<=16 MB) at a time, so a multi-GB basket never sits in memory.
    `counts` is the per-entry length (the counter branch). Each basket holds its
    entries' values followed by an entry-offset table, which is skipped."""
    import struct
    import zlib
    import lzma
    nb = branch.num_baskets
    seeks = [int(x) for x in branch.member("fBasketSeek")[:nb]]
    eoff = np.asarray(branch.entry_offsets)                 # entry boundaries of baskets
    csum = np.concatenate([[0], np.cumsum(counts)])
    with open(path, "rb") as f:
        for b, seek in enumerate(seeks):
            remaining = int(csum[eoff[b + 1]] - csum[eoff[b]])   # values in this basket
            if remaining == 0:
                continue
            f.seek(seek)
            nbytes, _v, objlen, _d, keylen = struct.unpack(">ihiIh", f.read(16))
            f.seek(seek + keylen)
            left = nbytes - keylen
            compressed = left != objlen
            buf = b""
            while left > 0 and remaining > 0:
                if not compressed:
                    raw = f.read(left)
                    left = 0
                else:
                    hdr = f.read(9)
                    algo = hdr[:2]
                    csz = int.from_bytes(hdr[3:6], "little")
                    usz = int.from_bytes(hdr[6:9], "little")
                    data = f.read(csz)
                    left -= 9 + csz
                    if algo == b"ZL":
                        raw = zlib.decompress(data)
                    elif algo == b"XZ":
                        raw = lzma.decompress(data)
                    elif algo == b"L4":
                        import lz4.block
                        raw = lz4.block.decompress(data[8:], uncompressed_size=usz)
                    elif algo == b"ZS":
                        import zstandard
                        raw = zstandard.ZstdDecompressor().decompress(data, max_output_size=usz)
                    else:
                        raise RuntimeError(f"unsupported ROOT compression {algo!r}")
                buf += raw
                nvals = min(len(buf) // 4, remaining)
                if nvals >= chunk_values or nvals == remaining or left <= 0:
                    if nvals:
                        yield np.frombuffer(buf[:nvals * 4], dtype=">i4").astype(np.int32)
                    buf = buf[nvals * 4:]
                    remaining -= nvals
            if remaining:
                raise RuntimeError(f"{branch.name}: basket {b} ended with {remaining} values missing")


def _fit_batch(prof):
    x, sx, cx, ox = process_profiles(prof[:, :NWIRE])
    y, sy, cy, oy = process_profiles(prof[:, NWIRE:])
    return sx, sy, cx, cy, ox & oy


def fit_device_streaming(path, tree, dev, jobs=1, batch=40000, max_entries=None):
    """Fit every 96-channel profile of `dev` without loading the branch.
    Returns ({sx, sy, chi2x, chi2y}, ok, has) arrays over all entries."""
    import concurrent.futures as cf
    import uproot
    t = uproot.open(path)[tree]
    counts = t["n" + dev].array(library="np").astype(np.int64)
    n = len(counts) if not max_entries else min(len(counts), max_entries)
    has = counts[:n] == 2 * NWIRE
    res = {k: np.full(n, np.nan) for k in ("sx", "sy", "chi2x", "chi2y")}
    ok = np.zeros(n, bool)
    if not has.any():
        return res, ok, has
    starts = np.concatenate([[0], np.cumsum(counts)[:-1]])     # value offset of each entry
    rows = np.where(has)[0]
    pending, pos, done = [], 0, 0                              # pos: values consumed so far
    carry = np.zeros(0, np.int32)
    ex = cf.ProcessPoolExecutor(jobs) if jobs > 1 else None

    def collect_one():
        nonlocal done
        idx, fut = pending.pop(0)
        sx, sy, cx, cy, o = fut.result() if ex else fut
        res["sx"][idx], res["sy"][idx], res["chi2x"][idx], res["chi2y"][idx] = sx, sy, cx, cy
        ok[idx] = o
        done += len(idx)
        print(f"\r[fit] {tree}/{dev}: {done}/{len(rows)} profiles", end="", file=sys.stderr, flush=True)

    r0 = 0                                                     # next row (index into rows) to emit
    for vals in stream_jagged_int32(path, t[dev], 96 * batch, counts):
        data = np.concatenate([carry, vals]) if len(carry) else vals
        base = pos - len(carry)                                # value offset of data[0]
        end = base + len(data)
        r1 = r0
        while r1 < len(rows) and starts[rows[r1]] + 2 * NWIRE <= end:
            r1 += 1
        if r1 > r0:
            idx = rows[r0:r1]
            take = (starts[idx] - base)[:, None] + np.arange(2 * NWIRE)[None, :]
            prof = data[take].astype(float)
            pending.append((idx, ex.submit(_fit_batch, prof) if ex else _fit_batch(prof)))
            keep_from = (starts[rows[r1]] - base) if r1 < len(rows) else len(data)
            r0 = r1
        else:
            keep_from = (starts[rows[r0]] - base) if r0 < len(rows) else len(data)
        carry = data[max(0, keep_from):].copy()
        pos += len(vals)
        while len(pending) > max(2, 2 * jobs):
            collect_one()
        if r0 >= len(rows):
            break                                              # all wanted rows queued
    while pending:
        collect_one()
    if ex:
        ex.shutdown()
    print(file=sys.stderr)
    if done != len(rows):
        raise RuntimeError(f"{dev}: fitted {done} of {len(rows)} profiles -- stream/offset mismatch")
    return res, ok, has


def fill_offsets_fast(run, t, offsets):
    """Nearest-in-time valid value within the same run for each -999 offset (vectorized).
    Returns (filled dict, mask of rows where anything was filled)."""
    import pandas as pd
    order = np.lexsort((t, run))
    rs, ts = run[order], t[order]
    filled_any = np.zeros(len(run), bool)
    out = {}
    for c, v in offsets.items():
        v = np.asarray(v, float)
        vs = v[order]
        good = _valid(vs)
        idx = np.arange(len(vs), dtype=float)
        gi = pd.Series(np.where(good, idx, np.nan))
        grp = pd.Series(rs)
        prev = gi.groupby(grp).ffill().to_numpy()
        nxt = gi.groupby(grp).bfill().to_numpy()
        tp = np.where(np.isfinite(prev), ts[np.nan_to_num(prev).astype(int)], np.inf)
        tn = np.where(np.isfinite(nxt), ts[np.nan_to_num(nxt).astype(int)], np.inf)
        use_prev = np.abs(ts - tp) <= np.abs(tn - ts)
        src = np.where(use_prev, prev, nxt)
        fill = ~good & np.isfinite(src)
        vs2 = vs.copy()
        vs2[fill] = vs[src[fill].astype(int)]
        res = np.empty_like(v)
        res[order] = vs2
        out[c] = res
        f = np.zeros(len(v), bool)
        f[order] = fill
        filled_any |= f
    return out, filled_any


def _stream_chunk(args):
    cols, fill_cols, filled, fits, mw_max_dt, min_tor = args
    import pandas as pd
    df = pd.DataFrame(cols)
    leg = compute_fom(df, None, mode="legacy", fit_cache={d: v for d, v in fits.items()})
    dff = df.copy()
    for c, v in fill_cols.items():
        dff[c] = v
    fix = compute_fom(dff, None, mode="fixed", mw_max_dt=mw_max_dt, min_tor=min_tor,
                      offsets_filled=filled, fit_cache={d: v for d, v in fits.items()})
    out = {**leg, **fix}
    for k in list(out):
        if out[k].dtype == np.float64 and k not in ("FOM_legacy",):
            out[k] = out[k].astype(np.float32) if not k.startswith(("FOM", "hpos", "vpos", "hang", "vang")) else out[k]
    return out


def process_tree_streaming(path, tree, writer, a):
    import concurrent.futures as cf
    import uproot
    t = uproot.open(path)[tree]
    n = t.num_entries if not a.max_entries else min(t.num_entries, a.max_entries)
    keys = set(t.keys())
    print(f"\n=== {tree}: {n} entries ===", file=sys.stderr)

    fits = {d: fit_device_streaming(path, tree, d, jobs=a.jobs, max_entries=n) for d in ("MMBTBB", "M876BB", "M875BB")
            if d in keys}
    for d in ("MMBTBB", "M876BB", "M875BB"):
        fits.setdefault(d, ({k: np.full(n, np.nan) for k in ("sx", "sy", "chi2x", "chi2y")},
                            np.zeros(n, bool), np.zeros(n, bool)))

    want = [c for c in FOM_INPUTS + ID_COLS if c in keys]
    cols = {c: t[c].array(library="np", entry_stop=n) for c in want}
    tt = cols["spill_time_sec"].astype(float) + 1e-9 * cols["spill_time_nsec"].astype(float)
    if a.fill_offsets:
        fill_cols, filled = fill_offsets_fast(cols["run"], tt, {c: cols[c] for c in OFFSET_COLS})
    else:
        fill_cols, filled = {c: cols[c] for c in OFFSET_COLS}, None

    edges = list(range(0, n, a.chunk)) + [n]

    def job(s, e):
        sl = slice(s, e)
        return ({c: v[sl] for c, v in cols.items()}, {c: v[sl] for c, v in fill_cols.items()},
                None if filled is None else filled[sl],
                {d: ({k: v[sl] for k, v in r.items()}, o[sl], h[sl]) for d, (r, o, h) in fits.items()},
                a.mw_max_dt, a.min_tor)

    ex = cf.ProcessPoolExecutor(a.jobs) if a.jobs > 1 else None
    futs, first, summary = [], True, []
    try:
        for i, (s, e) in enumerate(zip(edges[:-1], edges[1:])):
            futs.append((s, e, ex.submit(_stream_chunk, job(s, e)) if ex else _stream_chunk(job(s, e))))
            while len(futs) > (2 * a.jobs if ex else 0) or (i == len(edges) - 2 and futs):
                s0, e0, f0 = futs.pop(0)
                res = f0.result() if ex else f0
                out = {c: cols[c][s0:e0] for c in want}
                for c in ("has_spill_info", "triggered"):
                    if c in out:
                        out[c] = out[c].astype(np.int8)
                out.update(res)
                if first:      # explicit TTree (recent uproot would otherwise write an RNTuple)
                    wtree = writer.mktree(tree, {c: np.asarray(v).dtype for c, v in out.items()})
                    first = False
                wtree.extend(out)
                summary.append({k: out[k] for k in ("FOM", "FOM_legacy", "FOM_best", "fom_source", "fom_status")
                                if k in out})
                print(f"\r[fom] {tree}: {e0}/{n}", end="", file=sys.stderr, flush=True)
    finally:
        if ex:
            ex.shutdown()
    print(file=sys.stderr)
    import pandas as pd
    summarize(tree, pd.DataFrame({k: np.concatenate([x[k] for x in summary]) for k in summary[0]}))


def _one_chunk(args):
    df, dff, mw, filled, mw_max_dt, min_tor = args
    cache = {}
    leg = compute_fom(df, mw, mode="legacy", fit_cache=cache)
    fix = compute_fom(dff, mw, mode="fixed", mw_max_dt=mw_max_dt, min_tor=min_tor,
                      offsets_filled=filled, fit_cache=cache)
    return {**leg, **fix}


def _run_chunks(df, dff, mw, filled, a):
    """split rows into chunks (offset filling already done on the full table) and
    process them, optionally in parallel; returns concatenated columns."""
    import awkward as ak
    import concurrent.futures as cf
    n = len(df)
    mwa = {d: (v if isinstance(v, ak.Array) else ak.Array([np.asarray(x).tolist() for x in v]))
           for d, v in mw.items()}
    edges = list(range(0, n, a.chunk)) + [n]
    jobs = [(df.iloc[s:e], dff.iloc[s:e], {d: v[s:e] for d, v in mwa.items()},
             None if filled is None else filled[s:e], a.mw_max_dt, a.min_tor)
            for s, e in zip(edges[:-1], edges[1:])]
    if a.jobs > 1 and len(jobs) > 1:
        with cf.ProcessPoolExecutor(a.jobs) as ex:
            res = list(ex.map(_one_chunk, jobs))
    else:
        res = []
        for i, j in enumerate(jobs, 1):
            res.append(_one_chunk(j))
            print(f"\r[fom] chunk {i}/{len(jobs)}", end="", file=sys.stderr, flush=True)
        print(file=sys.stderr)
    return {k: np.concatenate([r[k] for r in res]) for k in res[0]}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help="dump_bnb_quality.py output (.root, or the *_spill.parquet)")
    ap.add_argument("-o", "--output", required=True, help="output (.root or .parquet)")
    ap.add_argument("--tree", choices=["spill", "bnbinfo", "both"], default="both")
    ap.add_argument("--mw-max-dt", type=float, default=0.0333,
                    help="max |multiwire - spill| time [s] for a profile to be used (fixed mode)")
    ap.add_argument("--min-tor", type=float, default=0.1e12, help="LOW_TOR flag threshold [protons]")
    ap.add_argument("--fill-offsets", action=argparse.BooleanOptionalAction, default=True,
                    help="fill missing BPM offsets from the nearest-in-time valid spill of the same run "
                         "(default: on; --no-fill-offsets to disable)")
    ap.add_argument("-j", "--jobs", type=int, default=1, help="parallel worker processes")
    ap.add_argument("--chunk", type=int, default=100000, help="spills per work unit")
    ap.add_argument("--max-entries", type=int, default=0, help="only the first N entries per tree (testing)")
    ap.add_argument("--in-memory", action="store_true",
                    help="load whole trees (small files); default for .root is the streaming path")
    a = ap.parse_args(argv)

    trees = ["spill", "bnbinfo"] if a.tree == "both" else [a.tree]
    if a.input.endswith(".root") and a.output.endswith(".root") and not a.in_memory:
        import uproot
        with uproot.recreate(a.output) as w:
            for tr in trees:
                process_tree_streaming(a.input, tr, w, a)
            w["fom_status_bits"] = "\n".join(f"{v}:{k}" for k, v in STATUS.items())
            w["fom_sources"] = "\n".join(f"{k}:{v}" for k, v in SOURCES.items())
        print(f"\nwrote {a.output}", file=sys.stderr)
        return
    frames = {}
    for tr in trees:
        df, mw = read_dump(a.input, tr)
        if len(df) == 0:
            continue
        dff, filled = (fill_offsets(df) if a.fill_offsets else (df, None))
        cols = _run_chunks(df, dff, mw, filled, a)
        for k, v in cols.items():
            df[k] = v
        frames[tr] = df
        summarize(tr, df)
    write_out(a.output, frames)
    print(f"\nwrote {a.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
