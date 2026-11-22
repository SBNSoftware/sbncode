#!/usr/bin/env python3
"""
bnb_fom_fast.py -- drop-in, faster front end for bnb_fom.py.

bnb_fom.py is left untouched (it stays the reference).  This script imports it,
swaps its two hot spots for numba-compiled, multithreaded kernels that do the
same arithmetic element by element, and then runs bnb_fom.main() with the same
command-line options:

  intbivar          target integration (~7k grid points x 3 positions per spill)
  process_profiles  Levenberg-Marquardt Gaussian fit of each multiwire profile

Same algorithm, same initial values, same LM step/acceptance/stop rules, same
pseudo-inverse; only floating-point summation order differs (results agree with
bnb_fom.py to ~1e-12 or better; see --check).

  ./bnb_fom_fast.py beamqual.root -o beamqual_fom.root            # all cores via threads
  ./bnb_fom_fast.py beamqual.root -o beamqual_fom.root -j 4       # 4 processes, cores split
  ./bnb_fom_fast.py --check beamqual.root                           # compare kernels vs bnb_fom

Option --fixups: afterwards apply all post-fixes of fom_fixups.py (includes the neighbor fill).
Option --neighbor-fill: afterwards fill single-spill BPM gaps from the two adjacent spills
(neighbor_fill.py; filled spills get fom_status bit 128 NEIGHBOR_FILLED).

Needs numba (pip install --user numba).  Threads per process: BNB_FOM_THREADS;
work-sharing granularity: BNB_FOM_CHUNK (rows per hand-out, 0 = equal split)
(default: cores // jobs).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import numba as nb
import bnb_fom as ref

NWIRE = ref.NWIRE
_DXDY = float(ref.DXDY)
_XC = np.ascontiguousarray(ref._XC)
_ORIG = {"intbivar": ref.intbivar, "process_profiles": ref.process_profiles}

_JIT = dict(nopython=True, parallel=True, cache=True, error_model="numpy", fastmath=False)


# --------------------------------------------------------------------------- target integral
def _grid_runs():
    """GX_PTS/GY_PTS are a masked 96x96 meshgrid ('ij'): for each x column the kept
    points are one contiguous run of the y values.  Store (x, first y, last y+1)."""
    ys = np.unique(ref.GY_PTS)
    ux, st, en = [], [], []
    for v in np.unique(ref.GX_PTS):
        m = np.where(ref.GX_PTS == v)[0]
        j0 = np.searchsorted(ys, ref.GY_PTS[m[0]])
        j1 = j0 + len(m)
        assert np.all(np.diff(m) == 1) and np.array_equal(ys[j0:j1], ref.GY_PTS[m])
        ux.append(v); st.append(j0); en.append(j1)
    assert sum(e - s for s, e in zip(st, en)) == len(ref.GX_PTS)
    return np.array(ux), np.array(st, np.int64), np.array(en, np.int64), ys


_UX, _ST, _EN, _YS = _grid_runs()


@nb.jit(**_JIT)
def _intbivar_kernel(cx, cy, sx, sy, rho, pi, ux, st, en, ys, dxdy):
    # identical per-point arithmetic to bnb_fom.intbivar (tx, ty, z, exp(-z/den));
    # tx/ty are computed once per grid column/row instead of once per point.
    n = cx.shape[0]
    out = np.empty(n)
    for i in nb.prange(n):
        r = rho[i]
        den = 2.0 * (1.0 - r * r)
        ty = (ys - cy[i]) / sy[i]
        tot = 0.0
        for a in range(ux.shape[0]):
            tx = (ux[a] - cx[i]) / sx[i]
            for b in range(st[a], en[a]):
                z = tx * tx - 2.0 * r * tx * ty[b] + ty[b] * ty[b]
                tot += np.exp(-z / den)
        tot = tot * dxdy / (2.0 * pi * sx[i] * sy[i] * np.sqrt(1.0 - r * r))
        out[i] = -10000.0 if tot >= 1.0 else np.log10(1.0 - tot)
    return out


def intbivar(cx, cy, sx, sy, rho, pi=3.14159):
    f = lambda a: np.ascontiguousarray(a, dtype=np.float64)
    return _intbivar_kernel(f(cx), f(cy), f(sx), f(sy), f(rho), float(pi), _UX, _ST, _EN, _YS, _DXDY)


# --------------------------------------------------------------------------- profile fit
@nb.njit(cache=True, error_model="numpy")
def _inv3(H, out):
    """inverse of the symmetric 3x3 LM matrix by cofactors; falls back to the SVD
    pseudo-inverse (what bnb_fom uses) whenever the matrix is not comfortably
    non-singular, so the two agree to rounding."""
    a, b, c = H[0, 0], H[0, 1], H[0, 2]
    d, e, f = H[1, 1], H[1, 2], H[2, 2]
    A = d * f - e * e
    B = c * e - b * f
    C = b * e - c * d
    det = a * A + b * B + c * C
    nrm = abs(a) + abs(d) + abs(f) + 2.0 * (abs(b) + abs(c) + abs(e))
    if not (np.isfinite(det) and abs(det) > 1e-9 * nrm * nrm * nrm):
        return np.linalg.pinv(H)
    out[0, 0] = A / det
    out[0, 1] = out[1, 0] = B / det
    out[0, 2] = out[2, 0] = C / det
    out[1, 1] = (a * f - c * c) / det
    out[1, 2] = out[2, 1] = (b * c - a * e) / det
    out[2, 2] = (a * d - b * b) / det
    return out


@nb.jit(**_JIT)
def _profiles_kernel(mw, XC, n_iter):
    n, W = mw.shape
    ox = np.empty(n)
    osx = np.empty(n)
    oc = np.empty(n)
    ok = np.zeros(n, np.bool_)
    for i in nb.prange(n):
        y = np.empty(W)
        inr = np.zeros(W, np.bool_)
        mn = 0.0
        mx = 0.0
        for j in range(W):
            v = -mw[i, j]
            if v < mn:
                mn = v
            if v > mx:
                mx = v
        rng = mx - mn
        thr = rng * 0.2
        err = rng * 0.02
        sumy = 0.0
        first = -1
        last = -1
        for j in range(W):
            y[j] = -mw[i, j] - mn
            sumy += y[j]
            if y[j] > thr:
                if first < 0:
                    first = j
                last = j + 1
        lo = -12 + first * 0.5
        hi = -12 + last * 0.5
        nbin = 0
        ws = 0.0
        wx = 0.0
        A = 0.0
        for j in range(W):
            inr[j] = (XC[j] >= lo) and (XC[j] <= hi)
            if inr[j]:
                nbin += 1
                ws += y[j]
                wx += y[j] * XC[j]
                if y[j] > A:
                    A = y[j]
        mu = wx / ws
        wv = 0.0
        for j in range(W):
            if inr[j]:
                d = XC[j] - mu
                wv += y[j] * d * d
        v2 = wv / ws
        if v2 < 0.0625:          # np.maximum(v, 0.0625) keeps NaN
            v2 = 0.0625
        sg = np.sqrt(v2)
        e2 = err * err if err > 0 else np.nan
        sa = np.sqrt(e2)
        p0, p1, p2 = A, mu, sg
        c2 = 0.0
        for j in range(W):
            if inr[j]:
                u = (XC[j] - p1) / p2
                d = y[j] - p0 * np.exp(-0.5 * u * u)
                c2 += d * d / e2
        lam = 1e-3
        if np.isfinite(c2) and nbin > 3:
            H = np.empty((3, 3))
            Pi = np.empty((3, 3))
            for _ in range(n_iter):
                J00 = J01 = J02 = J11 = J12 = J22 = 0.0
                r0 = r1 = r2 = 0.0
                for j in range(W):
                    if inr[j]:
                        u = (XC[j] - p1) / p2
                        ex = np.exp(-0.5 * u * u)
                        g = p0 * ex
                        a0 = ex / sa
                        a1 = g * u / p2 / sa
                        a2 = g * u * u / p2 / sa
                        rr = (y[j] - g) / sa
                        J00 += a0 * a0; J01 += a0 * a1; J02 += a0 * a2
                        J11 += a1 * a1; J12 += a1 * a2; J22 += a2 * a2
                        r0 += a0 * rr; r1 += a1 * rr; r2 += a2 * rr
                H[0, 0] = J00 + lam * J00; H[1, 1] = J11 + lam * J11; H[2, 2] = J22 + lam * J22
                H[0, 1] = H[1, 0] = J01
                H[0, 2] = H[2, 0] = J02
                H[1, 2] = H[2, 1] = J12
                bad = not (np.isfinite(H).all() and np.isfinite(r0) and np.isfinite(r1) and np.isfinite(r2))
                if bad:
                    break        # vectorized version: step 0, not better, done
                Pi = _inv3(H, Pi)
                q0 = p0 + (Pi[0, 0] * r0 + Pi[0, 1] * r1 + Pi[0, 2] * r2)
                q1 = p1 + (Pi[1, 0] * r0 + Pi[1, 1] * r1 + Pi[1, 2] * r2)
                q2 = p2 + (Pi[2, 0] * r0 + Pi[2, 1] * r1 + Pi[2, 2] * r2)
                c2n = 0.0
                for j in range(W):
                    if inr[j]:
                        u = (XC[j] - q1) / q2
                        d = y[j] - q0 * np.exp(-0.5 * u * u)
                        c2n += d * d / e2
                better = c2n < c2
                c2a = c2
                if better:
                    p0, p1, p2 = q0, q1, q2
                    c2 = c2n
                    lam = lam * 0.3
                else:
                    lam = lam * 10.0
                if (better and (c2a - c2n < 1e-12 * c2a)) or lam > 1e10:
                    break
        ndf = nbin - 3
        chi2ndf = c2 / ndf
        sumw = sumy > 0
        ok[i] = sumw and ndf > 0 and np.isfinite(chi2ndf) and np.isfinite(p0) \
            and np.isfinite(p1) and np.isfinite(p2)
        ox[i] = p1 if sumw else 99999.
        osx[i] = p2 if sumw else 99999.
        oc[i] = chi2ndf if sumw else 99999.
    return ox, osx, oc, ok


def process_profiles(mw, n_iter=200):
    mw = np.ascontiguousarray(mw, dtype=np.float64)
    return _profiles_kernel(mw, _XC, int(n_iter))


# --------------------------------------------------------------------------- install
def _threads():
    t = os.environ.get("BNB_FOM_THREADS")
    return max(1, min(int(t), nb.config.NUMBA_NUM_THREADS)) if t else None


def install():
    t = _threads()
    if t:
        nb.set_num_threads(t)
    # BNB_FOM_CHUNK=k: threads take k rows at a time instead of one equal share each, so
    # fast cores are not left waiting for slow ones (Apple P/E cores, busy shared nodes).
    # Results are identical either way (rows are independent).
    c = int(os.environ.get("BNB_FOM_CHUNK", "0") or 0)
    if c > 0:
        nb.set_parallel_chunksize(c)
    ref.intbivar = intbivar
    ref.process_profiles = process_profiles


install()     # at import time, so spawned/forked worker processes are patched too


def _check(path, n=20000):
    """compare the kernels with the bnb_fom originals on real profiles."""
    import time
    import uproot
    t = uproot.open(path)["bnbinfo"]
    for dev in ("M875BB", "M876BB", "MMBTBB"):
        if dev not in t.keys():
            continue
        counts = t["n" + dev].array(library="np").astype(np.int64)
        vals = next(ref.stream_jagged_int32(path, t[dev], 96 * n, counts))
        c = counts[np.cumsum(counts) <= len(vals)]
        starts = np.concatenate([[0], np.cumsum(c)[:-1]])
        rows = starts[c == 2 * NWIRE]
        if len(rows) == 0:
            continue
        P = vals[rows[:, None] + np.arange(2 * NWIRE)[None, :]].astype(float)[:n]
        for half in (P[:, :NWIRE], P[:, NWIRE:]):
            t0 = time.time(); r0 = _ORIG["process_profiles"](half); t1 = time.time()
            r1 = process_profiles(half); t2 = time.time()
            same_ok = (r0[3] == r1[3]).mean()
            m = r0[3] & r1[3]
            dmax = max(np.nanmax(np.abs(r0[k][m] - r1[k][m]) / np.maximum(1, np.abs(r0[k][m]))) for k in range(3))
            print(f"{dev}: {len(half)} profiles  ok-agree {same_ok:.6f}  max rel diff {dmax:.2e}  "
                  f"numpy {t1-t0:.2f}s  numba {t2-t1:.2f}s", flush=True)
    rng = np.random.default_rng(1)
    k = 20000
    cx, cy = rng.normal(0, 1.5, k), rng.normal(0, 1.5, k)
    sx, sy, rho = rng.uniform(0.5, 3, k), rng.uniform(0.5, 3, k), rng.uniform(-0.3, 0.3, k)
    t0 = time.time()
    f0 = np.concatenate([_ORIG["intbivar"](cx[s:s + 2000], cy[s:s + 2000], sx[s:s + 2000], sy[s:s + 2000],
                                           rho[s:s + 2000]) for s in range(0, k, 2000)])
    t1 = time.time()
    f1 = intbivar(cx, cy, sx, sy, rho); t2 = time.time()
    print(f"intbivar: max |diff| {np.nanmax(np.abs(f0 - f1)):.2e}  numpy {t1-t0:.2f}s  numba {t2-t1:.2f}s")


if __name__ == "__main__":
    argv = sys.argv[1:]
    if argv and argv[0] == "--check":
        _check(argv[1])
        sys.exit(0)
    jobs = 1
    for i, x in enumerate(argv):
        if x in ("-j", "--jobs"):
            jobs = int(argv[i + 1])
        elif x.startswith("--jobs="):
            jobs = int(x.split("=")[1])
    if "BNB_FOM_THREADS" not in os.environ:
        os.environ["BNB_FOM_THREADS"] = str(max(1, (os.cpu_count() or 1) // max(1, jobs)))
    if jobs > 1:
        # numba's thread pool does not survive fork(); spawned workers re-import this
        # script, which re-applies the patch (see install() above)
        import multiprocessing
        multiprocessing.set_start_method("spawn", force=True)
    nf = "--neighbor-fill" in argv
    fx = "--fixups" in argv
    nf_groups = 1
    argv = [x for x in argv if x not in ("--neighbor-fill", "--fixups")]
    if "--neighbor-fill-groups" in argv:
        k = argv.index("--neighbor-fill-groups")
        nf_groups = int(argv[k + 1])
        del argv[k:k + 2]
    install()
    ref.main(argv)
    if fx:
        # all post-fixes (TOR875 back-up, fallback removal, empty-multiwire / MMBTBB widths,
        # neighbor fill); see fom_fixups.py
        import fom_fixups
        out = argv[argv.index("-o") + 1] if "-o" in argv else argv[argv.index("--output") + 1]
        fom_fixups.run(out, out + ".fx.root", nf_groups)
        os.replace(out + ".fx.root", out)
    elif nf:
        # fill single-spill BPM gaps from the two adjacent spills (see neighbor_fill.py);
        # the output is rewritten in place, filled spills carry fom_status bit 128
        import neighbor_fill
        out = argv[argv.index("-o") + 1] if "-o" in argv else argv[argv.index("--output") + 1]
        neighbor_fill.run(out, out + ".nf.root", nf_groups)
        os.replace(out + ".nf.root", out)
