#!/usr/bin/env python3
"""
logger_fill.py -- fill the FOM for spills in the June 2023 ICARUS outage (runs 9941-9946)
using once-a-minute ACNet logger readings (node 96 "Backup") of the BPMs IFBeam lost.

In that outage HP875, VP875 and VP873 have no IFBeam reading, but the target BPMs HPTG1
and VPTG2 read on every spill.  ACNet kept one HP875 and one VP875 reading per minute;
each is the reading of one real spill (checked bit-for-bit against IFBeam).  For each
logger sample we take the target BPM of that same spill and form the beam angle:

    hang = (HPTG1 - HP875) / (z_HPTG1 - z_HP875)     vang = (VPTG2 - VP875) / (z_VPTG2 - z_VP875)

and every outage spill uses its own HPTG1 / VPTG2 plus the angle of the nearest sample
(<= 60 s away) to get the position at the target:

    hpos = HPTG1 + hang * (z_tgt - z_HPTG1)          vpos = VPTG2 + vang * (z_tgt - z_VPTG2)

Raw readings are used; the (constant) BPM offsets and the VPTG2 bias are one constant per
quantity, measured from good spills within 3 days of the outage (`calibrate`).  The same
code, run on good spills with one emulated sample per minute, is the closure test.

  ./logger_fill.py calibrate icarus_run2_full_fixed.root     # writes the calibration json
  (data and calibration: acnet/ next to this file; set BQ_ACNET_DIR to a writable directory
   for `calibrate` when the code is installed read-only -- files there take precedence)
  (the fill itself runs inside fom_fixups.py --logger-fill; bit 8192 LOGGER_FILL)

The March 2023 outage (runs 9648-9649) is NOT filled: its closure test accepted ~0.2% bad
spills (the angle has to be carried 4.75 m from HP875).
"""
import json
import os
import sys

import numpy as np
import uproot

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bnb_fom as ref                                          # noqa: E402

Z = ref.Z
OUTAGE = "A_2023-06-09to12"                     # the only outage filled
WINDOW = (1686344400.0, 1686592800.0)           # 2023-06-09 21:00 -> 2023-06-12 18:00 UTC
MAX_DT = 60.0                                   # spill -> logger sample, s
MATCH_DT = 1.5                                  # logger timestamp -> its spill (spill precedes it), s
EMPTY_MW = 4.0
BIT = 8192
HERE = os.path.dirname(os.path.abspath(__file__))
PKG_ACNET = os.path.join(HERE, "acnet")                  # data shipped with the code (read-only once installed)
ACNET = os.environ.get("BQ_ACNET_DIR", PKG_ACNET)        # where `calibrate` writes; read from first if present


def _acnet_file(name):
    p = os.path.join(ACNET, name)
    return p if os.path.exists(p) else os.path.join(PKG_ACNET, name)


CSV = _acnet_file("recovered_1min_backup_node96.csv")
CALIB = _acnet_file("logger_fill_calib.json")


def ok(v):
    v = np.asarray(v, float)
    return np.isfinite(v) & (v > -998)


def load_samples(csv=CSV):
    rows = [l.rstrip("\n").split(",") for l in open(csv)][1:]
    S = {}
    for dev, out, ms, _, val in rows:
        if out == OUTAGE:
            S.setdefault(dev, []).append((int(ms) / 1e3, float(val)))
    t1, h = np.array(S["E:HP875"]).T
    t2, v = np.array(S["E:VP875"]).T
    assert np.array_equal(t1, t2), "HP875 and VP875 samples are expected at the same times"
    return dict(t=t1, HP875=h, VP875=v)


def sample_angles(t, run, hptg1, vptg2, smp):
    """match each logger sample to its spill (last spill at or before it, within MATCH_DT,
    with HPTG1 and VPTG2 read) and return sample time, run and raw angles"""
    good = ok(hptg1) & ok(vptg2)
    gi = np.where(good)[0]
    o = gi[np.argsort(t[gi])]
    ts = t[o]
    j = np.searchsorted(ts, smp["t"], side="right") - 1
    hit = (j >= 0) & (smp["t"] - ts[np.clip(j, 0, None)] <= MATCH_DT)
    k = o[j[hit]]
    ha = (hptg1[k] - smp["HP875"][hit]) / (Z["hptg1"] - Z["hp875"])
    va = (vptg2[k] - smp["VP875"][hit]) / (Z["vptg2"] - Z["vp875"])
    return smp["t"][hit], run[k], ha, va, hptg1[k], vptg2[k]


def estimate(t, run, hptg1, vptg2, st, sr, ha, va, bias):
    """position/angle at the target for spills (t, run) from their own HPTG1/VPTG2 and the
    nearest sample of the same run within MAX_DT; NaN where there is none"""
    o = np.argsort(st)
    st, sr, ha, va = st[o], sr[o], ha[o], va[o]
    E = np.full((len(t), 4), np.nan)
    if len(st) == 0:
        return E
    i = np.searchsorted(st, t)
    best = np.full(len(t), -1)
    bdt = np.full(len(t), np.inf)
    for c in (np.clip(i - 1, 0, len(st) - 1), np.clip(i, 0, len(st) - 1)):
        dt = np.abs(st[c] - t)
        use = (sr[c] == run) & (dt <= MAX_DT) & (dt < bdt)
        best[use], bdt[use] = c[use], dt[use]
    m = best >= 0
    b = best[m]
    E[m, 1], E[m, 3] = ha[b], va[b]
    E[m, 0] = hptg1[m] + ha[b] * (Z["tgt"] - Z["hptg1"])
    E[m, 2] = vptg2[m] + va[b] * (Z["tgt"] - Z["vptg2"])
    return E + np.asarray(bias)


def widths(sx, sy, m876hs, m876vs):
    empty = (ok(m876hs) & (m876hs > EMPTY_MW)) | (ok(m876vs) & (m876vs > EMPTY_MW))
    have = np.isfinite(sx) & np.isfinite(sy) & ~empty
    return np.where(have, sx, np.nan), np.where(have, sy, np.nan), np.where(have, 3, 5)


def fom_of(E, tor, sx, sy):
    return 1 - 10 ** ref.calc_fom(E[:, 0], E[:, 1], E[:, 2], E[:, 3], tor, sx, sy, pi=np.pi)


# ------------------------------------------------------------------ calibration + closure
def calibrate(path, out=None, pad=3 * 86400, nmax=150_000):
    br = ["run", "spill_time_sec", "spill_time_nsec", "TOR860", "TOR875", "THCURR", "FOM_best", "fom_source",
          "HP875", "VP875", "VP873", "HPTG1", "VPTG2", "hpos_tgt", "hang_tgt", "vpos_tgt", "vang_tgt",
          "tgt_sigx", "tgt_sigy", "M876HS", "M876VS"]
    a, b = WINDOW
    if out is None:                          # BQ_ACNET_DIR, not the installed code directory
        out = os.path.join(ACNET, "logger_fill_calib.json")
        os.makedirs(ACNET, exist_ok=True)
    cache = os.path.join(os.path.dirname(out), "logger_fill_window_cache.npz")
    if os.path.exists(cache):
        d = dict(np.load(cache))
    else:                                      # one pass over the file, kept for re-runs
        parts = {k: [] for k in br}
        for d in uproot.open(path)["bnbinfo"].iterate(br, step_size=1_000_000, library="np"):
            m = (d["spill_time_sec"] >= a - pad) & (d["spill_time_sec"] <= b + pad) & (d["TOR860"] > 1e11)
            for k in br:
                parts[k].append(d[k][m])
        d = {k: np.concatenate(v) for k, v in parts.items()}
        np.savez(cache, **d)
        print(f"[logger_fill] cached the +-3 day window in {cache}; re-run to finish")
        return None
    t = d["spill_time_sec"] + 1e-9 * d["spill_time_nsec"]
    run = d["run"]
    allb = ok(d["HP875"]) & ok(d["VP875"]) & ok(d["VP873"]) & ok(d["HPTG1"]) & ok(d["VPTG2"])
    sel = (d["TOR875"] > 1e11) & (d["THCURR"] >= 173) & (d["THCURR"] <= 175)
    good = sel & allb & np.isfinite(d["FOM_best"]) & (d["fom_source"] >= 3) & ((t < a) | (t > b))
    # emulated logger: the first fully read spill of each minute (per run) is the sample
    key = run.astype(np.int64) * 10**8 + (t // 60).astype(np.int64)
    gi = np.where(allb)[0]
    o = gi[np.lexsort((t[gi], key[gi]))]
    first = o[np.r_[True, key[o][1:] != key[o][:-1]]]
    ha = (d["HPTG1"][first] - d["HP875"][first]) / (Z["hptg1"] - Z["hp875"])
    va = (d["VPTG2"][first] - d["VP875"][first]) / (Z["vptg2"] - Z["vp875"])
    E0 = estimate(t, run, d["HPTG1"], d["VPTG2"], t[first], run[first], ha, va, np.zeros(4))
    Y = np.column_stack([d["hpos_tgt"], d["hang_tgt"], d["vpos_tgt"], d["vang_tgt"]])
    good &= np.all(np.isfinite(E0), axis=1)
    sx, sy, _ = widths(d["tgt_sigx"], d["tgt_sigy"], d["M876HS"], d["M876VS"])
    tor = d["TOR860"]
    res = {}
    for name, cal, test in (("before->after", good & (t < a), good & (t > b)),
                            ("after->before", good & (t > b), good & (t < a))):
        bias = np.median(Y[cal] - E0[cal], axis=0)
        g = np.where(test)[0]
        g = np.random.default_rng(1).choice(g, min(nmax, len(g)), replace=False)
        fe = fom_of(E0[g] + bias, tor[g], sx[g], sy[g])
        ft = fom_of(Y[g], tor[g], sx[g], sy[g])
        acc, tru = (fe > 0.98) & (fe <= 1), (ft > 0.98) & (ft <= 1)
        res[name] = dict(test_spills=int(len(g)), false_accepts=int(np.sum(acc & ~tru)),
                         false_rejects=int(np.sum(~acc & tru)), true_pass=float(tru.mean()),
                         dfom_99=float(np.percentile(np.abs(fe - ft), 99)),
                         dfom_max=float(np.max(np.abs(fe - ft))),
                         pos_err_90_mm=[float(np.percentile(np.abs(E0[g, 0] + bias[0] - Y[g, 0]), 90)),
                                        float(np.percentile(np.abs(E0[g, 2] + bias[2] - Y[g, 2]), 90))])
        print(f"[logger_fill] closure {name}: {res[name]}")
    bias = np.median(Y[good] - E0[good], axis=0)
    margin = max(v["dfom_max"] for v in res.values())
    C = dict(outage=OUTAGE, window_utc=list(WINDOW), bias=[float(x) for x in bias], margin=margin,
             bias_order=["hpos_mm", "hang_mrad", "vpos_mm", "vang_mrad"], calibration_spills=int(good.sum()),
             closure=res, source_file=os.path.basename(path))
    json.dump(C, open(out, "w"), indent=1)
    print(f"[logger_fill] wrote {out}; bias {np.round(bias, 4)}")
    return C


# ------------------------------------------------------------------ the fill (called by fom_fixups)
def load(csv=CSV, calib=CALIB):
    if not os.path.exists(calib):
        sys.exit(f"[logger_fill] {calib} missing: run `logger_fill.py calibrate <X_fixed.root>` first")
    C = json.load(open(calib))
    tot = sum(v["false_accepts"] for v in C["closure"].values())
    if tot > 0:
        print(f"[logger_fill] WARNING: the closure test accepted {tot} bad spills", file=sys.stderr)
    return dict(samples=load_samples(csv), calib=C)


# stability guards and acceptance (validated across all of ICARUS Run 2 with burst_closure: no bad spill
# accepted in ~180k test spills with the angle taken 60 s away on both sides)
DANG, DTB, ACCEPT = 0.2, 0.5, 0.999


def _bracket(t, run, xs, xr):
    """index of the last x <= t and first x > t in the sorted sample times xs (same run, <= MAX_DT)"""
    i = np.searchsorted(xs, t, side="right")
    p, n = np.clip(i - 1, 0, len(xs) - 1), np.clip(i, 0, len(xs) - 1)
    okp = (i > 0) & (xr[p] == run) & (t - xs[p] <= MAX_DT)
    okn = (i < len(xs)) & (xr[n] == run) & (xs[n] - t <= MAX_DT) & (xs[n] > t)
    return p, n, okp & okn


def fill_group(d, fom, st, src, pos, sx, sy, tor, beam, S, L):
    """June 2023 outage: own target BPMs + the beam angle from the ACNet samples on both sides"""
    t, run = d["t"], d["run"]
    a, b = L["calib"]["window_utc"]
    cand = (t >= a) & (t <= b) & beam & ~np.isfinite(fom) & ((st & 1) == 0) & \
           ~ok(d["HP875"]) & ~ok(d["VP875"]) & ok(d["HPTG1"]) & ok(d["VPTG2"])
    S["logger_cand_pot"] = S.get("logger_cand_pot", 0.0) + float(tor[cand].sum())
    if not cand.any():
        return
    H1, V2 = d["HPTG1"].astype(float), d["VPTG2"].astype(float)
    xs, xr, ha, va, h1s, v2s = sample_angles(t, run, H1, V2, L["samples"])
    o = np.argsort(xs); xs, xr, ha, va, h1s, v2s = xs[o], xr[o], ha[o], va[o], h1s[o], v2s[o]
    r = np.where(cand)[0]
    if not len(xs):
        return
    p, n, both = _bracket(t[r], run[r], xs, xr)
    stable = both & (np.abs(ha[p] - ha[n]) < DANG) & (np.abs(va[p] - va[n]) < DANG) & \
        (np.abs(H1[r] - 0.5 * (h1s[p] + h1s[n])) < DTB) & (np.abs(V2[r] - 0.5 * (v2s[p] + v2s[n])) < DTB)
    r, p, n = r[stable], p[stable], n[stable]
    if not len(r):
        return
    h_a, v_a = 0.5 * (ha[p] + ha[n]), 0.5 * (va[p] + va[n])
    E = np.column_stack([H1[r] + h_a * (Z["tgt"] - Z["hptg1"]), h_a, V2[r] + v_a * (Z["tgt"] - Z["vptg2"]), v_a]) + \
        np.asarray(L["calib"]["bias"])
    wx, wy, wsrc = widths(sx[r], sy[r], d["M876HS"][r], d["M876VS"][r])
    f = fom_of(E, tor[r], wx, wy)
    acc = f > ACCEPT                       # anything less clear-cut stays without a FOM
    S["logger_rejected_pot"] = S.get("logger_rejected_pot", 0.0) + float(tor[r][~acc].sum())
    r, E, wx, wy, wsrc, f = r[acc], E[acc], wx[acc], wy[acc], wsrc[acc], f[acc]
    pos["hpos_tgt"][r], pos["hang_tgt"][r], pos["vpos_tgt"][r], pos["vang_tgt"][r] = E.T
    sx[r], sy[r] = wx, wy
    fom[r] = f
    src[r] = wsrc
    st[r] = (st[r] & ~6) | BIT | np.where(wsrc == 5, 64, 0)
    S["logger_pot"] = S.get("logger_pot", 0.0) + float(tor[r].sum())
    S["logger_spills"] = S.get("logger_spills", 0) + int(len(r))


# ------------------------------------------------------------------ 875-station burst fill (ICARUS)
BURST_BIT = 16384
B_FMIN, B_ACCEPT = 0.999, 0.995


def burst_fill_group(d, fom, st, src, pos, sx, sy, tor, beam, S, measured):
    """Short HP875+VP875 dropouts: own HPTG1/VPTG2 + position offset and angle of the nearest
    measured spill on each side (<= 60 s, same run).  Both neighbors must have FOM > 0.999 and
    agree in angle (0.2 mrad); the spill's target BPMs must be within 0.5 mm of theirs; the
    estimate must be > 0.995.  BPM offsets and the VPTG2 bias cancel in the differences."""
    t, run = d["t"], d["run"]
    H1, V2 = d["HPTG1"].astype(float), d["VPTG2"].astype(float)
    allb = ok(d["HP875"]) & ok(d["VP875"]) & ok(d["VP873"]) & ok(d["HPTG1"]) & ok(d["VPTG2"])
    good = measured & allb & np.isfinite(fom)
    cand = beam & ~np.isfinite(fom) & ((st & 1) == 0) & ~ok(d["HP875"]) & ~ok(d["VP875"]) & ok(d["HPTG1"]) & ok(d["VPTG2"])
    S["burst_cand_pot"] = S.get("burst_cand_pot", 0.0) + float(tor[cand].sum())
    gi, ci = np.where(good)[0], np.where(cand)[0]
    if not len(gi) or not len(ci):
        return
    o = np.lexsort((t[gi], run[gi])); gi = gi[o]
    key = lambda idx: run[idx].astype(np.float64) * 1e10 + t[idx]
    gk, ck = key(gi), key(ci)
    j = np.searchsorted(gk, ck)
    A, B = gi[np.clip(j - 1, 0, len(gi) - 1)], gi[np.clip(j, 0, len(gi) - 1)]
    c = (j > 0) & (j < len(gi)) & (run[A] == run[ci]) & (run[B] == run[ci]) & (t[ci] - t[A] <= MAX_DT) & (t[B] - t[ci] <= MAX_DT)
    P = pos
    c &= (fom[A] > B_FMIN) & (fom[B] > B_FMIN) & (np.abs(P["hang_tgt"][A] - P["hang_tgt"][B]) < DANG) & \
        (np.abs(P["vang_tgt"][A] - P["vang_tgt"][B]) < DANG) & \
        (np.abs(H1[ci] - 0.5 * (H1[A] + H1[B])) < DTB) & (np.abs(V2[ci] - 0.5 * (V2[A] + V2[B])) < DTB)
    r, A, B = ci[c], A[c], B[c]
    if not len(r):
        return
    E = np.column_stack([0.5 * (P["hpos_tgt"][A] - H1[A] + P["hpos_tgt"][B] - H1[B]) + H1[r], 0.5 * (P["hang_tgt"][A] + P["hang_tgt"][B]),
                         0.5 * (P["vpos_tgt"][A] - V2[A] + P["vpos_tgt"][B] - V2[B]) + V2[r], 0.5 * (P["vang_tgt"][A] + P["vang_tgt"][B])])
    wx, wy, wsrc = widths(sx[r], sy[r], d["M876HS"][r], d["M876VS"][r])
    f = fom_of(E, tor[r], wx, wy)
    acc = f > B_ACCEPT
    r, E, wx, wy, wsrc, f = r[acc], E[acc], wx[acc], wy[acc], wsrc[acc], f[acc]
    for k, col in zip(("hpos_tgt", "hang_tgt", "vpos_tgt", "vang_tgt"), E.T):
        P[k][r] = col
    sx[r], sy[r] = wx, wy
    fom[r] = f
    src[r] = wsrc
    st[r] = (st[r] & ~6) | BURST_BIT | np.where(wsrc == 5, 64, 0)
    S["burst_pot"] = S.get("burst_pot", 0.0) + float(tor[r].sum())
    S["burst_spills"] = S.get("burst_spills", 0) + int(len(r))


if __name__ == "__main__":
    if len(sys.argv) >= 3 and sys.argv[1] == "calibrate":
        calibrate(sys.argv[2])
    else:
        print(__doc__)
