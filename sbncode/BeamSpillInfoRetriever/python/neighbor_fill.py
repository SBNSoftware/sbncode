#!/usr/bin/env python3
"""
neighbor_fill.py -- fill FOM_best for beam spills that lost it to a single missing BPM
reading, using the two adjacent spills (one 15 Hz cycle away on each side).

A spill is filled only if ALL of these hold (defaults; see neighbor_fill_study.py for the
validation that motivates them):
  * it has beam (TOR > 1e11) and FOM_best is NaN because a BPM reading is missing
    (fom_status has NO_HBPM or NO_VBPM, not NO_TOR)
  * the spill immediately before and the spill immediately after (same run) both have a
    valid FOM_best and are each <= --max-gap (0.1 s: one spill time at 15 Hz) away
  * the spill's intensity is within 3 % of the neighbors' mean, and the neighbors agree:
    beam position at the target within 1 mm, angle within 0.4 mrad, width within 0.1 mm,
    FOM within 0.001
The filled spill gets the neighbors' mean position/angle at the target, its own intensity
and (if it has one) its own width; FOM_best is recomputed with the fixed algorithm.  Its
fom_status gains bit 128 NEIGHBOR_FILLED (drop those spills with (fom_status & 128) == 0) and
fom_source is set to the neighbors' width source.
Both trees are updated (the 'spill' tree entries are matched to bnbinfo by run, subrun and
spill time).  Expected effect: SBND Run 1 p1 +0.37 % of POT, ICARUS Run 2 +0.20 %.

  ./neighbor_fill.py sbnd_run1_part1_fixed.root -o sbnd_run1_part1_fixed_nf.root
  ./bnb_fom_fast.py dump.root -o fom.root --neighbor-fill        (same, done after the FOM)
"""
import argparse
import os
import sys

import numpy as np
import uproot

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
NEIGHBOR_FILLED = 128
COLS = ["run", "subrun", "spill_time_sec", "spill_time_nsec", "TOR860", "TOR875", "FOM_best", "fom_status", "fom_source",
        "hpos_tgt", "hang_tgt", "vpos_tgt", "vang_tgt", "tgt_sigx", "tgt_sigy"]
CUTS = dict(max_gap=0.1, dtor=0.03, dpos=1.0, dang=0.4, dsig=0.1, dfom=0.001)


def keyhash(run, sub, sec, nsec):
    t = sec.astype(np.int64) * 1_000_000_000 + nsec.astype(np.int64)
    rs = (run.astype(np.uint64) << np.uint64(20)) + sub.astype(np.uint64)
    with np.errstate(over="ignore"):
        return t.astype(np.uint64) ^ (rs * np.uint64(0x9E3779B97F4A7C15))


def find_fills(path, cuts=CUTS, groups=1, verbose=True):
    """returns dict: h (sorted keys), FOM_best, hpos_tgt, hang_tgt, vpos_tgt, vang_tgt, tgt_sigx, tgt_sigy"""
    import bnb_fom_fast as F
    t = uproot.open(path)["bnbinfo"]
    out = []
    stats = dict(n_cand=0, n_fill=0, pot_cand=0.0, pot_fill=0.0, pot_all=0.0)
    for g in range(groups):
        parts = []
        for d in t.iterate(COLS, step_size=4_000_000, library="np"):
            s = (d["run"] % groups) == g
            parts.append({k: v[s] for k, v in d.items()})
        d = {k: np.concatenate([p[k] for p in parts]) for k in COLS}
        T = d["spill_time_sec"].astype(np.int64) * 1_000_000_000 + d["spill_time_nsec"].astype(np.int64)
        o = np.lexsort((T, d["run"]))
        d = {k: v[o] for k, v in d.items()}
        T = T[o]
        keep = np.ones(len(T), bool)
        keep[1:] = ~((T[1:] == T[:-1]) & (d["run"][1:] == d["run"][:-1]))
        d = {k: v[keep] for k, v in d.items()}
        ts = T[keep] / 1e9
        n = len(ts)
        if n < 3:
            continue
        t860, t875 = d["TOR860"].astype(float), d["TOR875"].astype(float)
        tor = np.where(np.isfinite(t860) & (t860 > 0), t860, np.where(np.isfinite(t875) & (t875 > 0), t875, 0.0))
        fom = d["FOM_best"].astype(float)
        st = d["fom_status"].astype(np.int64)
        beam = tor > 1e11
        valid = np.isfinite(fom)
        cand = beam & ~valid & ((st & 6) != 0) & ((st & 1) == 0)
        stats["n_cand"] += int(cand.sum())
        stats["pot_cand"] += float(tor[cand].sum())
        stats["pot_all"] += float(tor[beam].sum())
        i = np.where(cand)[0]
        i = i[(i > 0) & (i < n - 1)]
        p, q = i - 1, i + 1
        run = d["run"]
        ok = (run[p] == run[i]) & (run[q] == run[i]) & valid[p] & valid[q]
        ok &= (ts[i] - ts[p] <= cuts["max_gap"]) & (ts[q] - ts[i] <= cuts["max_gap"])
        tm = 0.5 * (tor[p] + tor[q])
        with np.errstate(invalid="ignore", divide="ignore"):
            ok &= np.abs(tor[i] - tm) / tm <= cuts["dtor"]
        f = lambda k: d[k].astype(float)
        dd = lambda k: np.abs(f(k)[p] - f(k)[q])
        ok &= (np.maximum(dd("hpos_tgt"), dd("vpos_tgt")) <= cuts["dpos"])
        ok &= (np.maximum(dd("hang_tgt"), dd("vang_tgt")) <= cuts["dang"])
        ok &= np.nan_to_num(np.maximum(dd("tgt_sigx"), dd("tgt_sigy")), nan=0.0) <= cuts["dsig"]
        ok &= np.abs(fom[p] - fom[q]) <= cuts["dfom"]
        i, p, q = i[ok], p[ok], q[ok]
        if not len(i):
            continue
        av = {k: 0.5 * (f(k)[p] + f(k)[q]) for k in ("hpos_tgt", "hang_tgt", "vpos_tgt", "vang_tgt", "tgt_sigx", "tgt_sigy")}
        sx = np.where(np.isfinite(f("tgt_sigx")[i]), f("tgt_sigx")[i], av["tgt_sigx"])
        sy = np.where(np.isfinite(f("tgt_sigy")[i]), f("tgt_sigy")[i], av["tgt_sigy"])
        f2 = F.ref.calc_fom(av["hpos_tgt"], av["hang_tgt"], av["vpos_tgt"], av["vang_tgt"], tor[i], sx, sy, pi=np.pi)
        r = {"h": keyhash(run[i], d["subrun"][i], d["spill_time_sec"][i], d["spill_time_nsec"][i]),
             "FOM_best": 1 - 10 ** f2, "hpos_tgt": av["hpos_tgt"], "hang_tgt": av["hang_tgt"],
             "vpos_tgt": av["vpos_tgt"], "vang_tgt": av["vang_tgt"], "tgt_sigx": sx, "tgt_sigy": sy,
             # width source: the neighbors' (their widths agree within 0.1 mm), so measured /
             # nominal categories stay meaningful; bit 128 marks the spill as filled
             "fom_source": d["fom_source"][p]}
        stats["n_fill"] += len(i)
        stats["pot_fill"] += float(tor[i].sum())
        out.append(r)
    if not out:
        return None, stats
    R = {k: np.concatenate([r[k] for r in out]) for k in out[0]}
    o = np.argsort(R["h"])
    R = {k: v[o] for k, v in R.items()}
    if verbose:
        print(f"[neighbor-fill] {stats['n_fill']:,} of {stats['n_cand']:,} BPM-missing beam spills filled: "
              f"{stats['pot_fill']/max(stats['pot_all'],1):.3%} of POT "
              f"(of {stats['pot_cand']/max(stats['pot_all'],1):.3%} missing)", file=sys.stderr)
    return R, stats


def apply(path, out, R, step=500_000):
    """copy the FOM file, writing the filled values (flat trees only)"""
    tmp = out + ".partial.root"
    with uproot.open(path) as fin, uproot.recreate(tmp) as fo:
        for tr in ("spill", "bnbinfo"):
            if tr not in fin:
                continue
            t = fin[tr]
            types = {k: t[k].interpretation.numpy_dtype for k in t.keys()}
            wt = fo.mktree(tr, types)
            nfill = 0
            for d in t.iterate(list(types), step_size=step, library="np"):
                if R is not None and len(R["h"]):
                    h = keyhash(d["run"], d["subrun"], d["spill_time_sec"], d["spill_time_nsec"])
                    j = np.clip(np.searchsorted(R["h"], h), 0, len(R["h"]) - 1)
                    hit = (R["h"][j] == h) & np.isnan(d["FOM_best"].astype(float))
                    if hit.any():
                        jj = j[hit]
                        for k in ("FOM_best", "hpos_tgt", "hang_tgt", "vpos_tgt", "vang_tgt", "tgt_sigx", "tgt_sigy",
                                  "fom_source"):
                            if k in d:
                                d[k] = d[k].copy()
                                d[k][hit] = R[k][jj].astype(d[k].dtype)
                        d["fom_status"] = d["fom_status"].copy()
                        d["fom_status"][hit] |= NEIGHBOR_FILLED
                        nfill += int(hit.sum())
                wt.extend({k: d[k].astype(types[k]) for k in types})
            print(f"[neighbor-fill] {tr}: {nfill:,} entries filled", file=sys.stderr)
        for k in fin.keys(cycle=False):
            if k in ("spill", "bnbinfo"):
                continue
            s = str(fin[k])
            if k == "fom_status_bits" and "NEIGHBOR_FILLED" not in s:
                s += f"\n{NEIGHBOR_FILLED}:NEIGHBOR_FILLED"
            fo[k] = s
    os.replace(tmp, out)


def run(path, out, groups=1, cuts=CUTS):
    R, stats = find_fills(path, cuts, groups)
    apply(path, out, R)
    return stats


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help="bnb_fom(_fast).py output")
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--groups", type=int, default=1, help="process runs in this many groups (less memory)")
    ap.add_argument("--max-gap", type=float, default=CUTS["max_gap"], help="max time to each neighbour [s]")
    a = ap.parse_args(argv)
    if os.path.abspath(a.output) == os.path.abspath(a.input):
        ap.error("output must differ from input")
    run(a.input, a.output, a.groups, {**CUTS, "max_gap": a.max_gap})
    print(f"wrote {a.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
