#!/usr/bin/env python3
"""Spill-by-spill check of the C++ beam-quality fields written to a flat CAF.

Compares rec.hdr.bnbinfo of a CAF made with the new retriever + CAFMaker
(FOM_best, fom_source, fom_status, tor_used, beam at the target, fom_version)
with
  * the production CAF of the same input (legacy FOM, device readings),
  * the analysis-level reference: bnb_fom.py (Tier 1) and fom_fixups.py
    (Tier 2) run on the production CAF,
and checks that rec.hdr.spillbnbinfo agrees with bnbinfo. Optionally checks a
CAF made by the new CAFMaker from an old-format BNBSpillInfo (fom_version 0).

    compare_cpp_caf.py --new new.flat.caf.root --prod prod.flat.caf.root \\
        --py1 prod_fom.root --py2 prod_fom_fx.root [--oldfmt oldfmt.flat.caf.root] -o report.md
"""
import argparse
import numpy as np
import awkward as ak
import uproot

FIELDS = ["spill_time_sec", "spill_time_nsec", "TOR860", "TOR875", "FOM", "HP875", "VP875", "VP873",
          "HPTG1", "HPTG2", "VPTG2", "HP875Offset", "VP875Offset", "VP873Offset", "HPTG1Offset",
          "HPTG2Offset", "VPTG2Offset", "M876HS", "M876VS", "M875HS", "M875VS",
          "M875BB_spill_time_diff", "M876BB_spill_time_diff"]
NEW = ["FOM_best", "fom_source", "fom_status", "tor_used", "hpos_tgt", "hang_tgt", "vpos_tgt",
       "vang_tgt", "tgt_sigx", "tgt_sigy", "fom_version"]
STATE = ["hpos_tgt", "hang_tgt", "vpos_tgt", "vang_tgt", "tgt_sigx", "tgt_sigy", "tor_used"]
FILL_BITS = 128 | 1024 | 16384
CPP_UNSET_BITS = 32 | 4096 | 8192   # set only by the analysis-level code (technote tab. cpp_py)


def key(sec, nsec):
    return np.asarray(sec, np.int64) * 1_000_000_000 + np.asarray(nsec, np.int64)


def read_caf(path, prefix, fields):
    t = uproot.open(path)["recTree"]
    have = set(t.keys())
    names = [f"{prefix}.{f}" for f in fields if f"{prefix}.{f}" in have]
    a = t.arrays(names, library="ak")
    out = {n.split(".")[-1]: ak.to_numpy(ak.flatten(a[n], axis=None)) if a[n].ndim > 1
           else ak.to_numpy(a[n]) for n in names}
    out["key"] = key(out["spill_time_sec"], out["spill_time_nsec"])
    return out


def unique_by_key(d):
    _, idx = np.unique(d["key"], return_index=True)
    return {k: v[idx] for k, v in d.items()}


def read_py(path):
    t = uproot.open(path)["bnbinfo"]
    cols = ["spill_time_sec", "spill_time_nsec", "FOM_best", "fom_source", "fom_status"] + STATE
    d = t.arrays([c for c in cols if c in t.keys()], library="np")
    for c in cols:   # tor_used is written only by fom_fixups.py
        if c not in d:
            d[c] = np.full(t.num_entries, np.nan)
    d["key"] = key(d["spill_time_sec"], d["spill_time_nsec"])
    return unique_by_key(d)


def match(a, b):
    """indices into a and b of the common keys"""
    common, ia, ib = np.intersect1d(a["key"], b["key"], return_indices=True)
    return ia, ib


def same(x, y, tol=0.0):
    x = np.asarray(x, float); y = np.asarray(y, float)
    both_nan = np.isnan(x) & np.isnan(y)
    return both_nan | (np.abs(x - y) <= tol)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--new", required=True)
    ap.add_argument("--prod", required=True)
    ap.add_argument("--py1", required=True)
    ap.add_argument("--py2", required=True)
    ap.add_argument("--oldfmt")
    ap.add_argument("-o", "--output", default="compare_cpp_caf.md")
    a = ap.parse_args(argv)
    R = []
    p = R.append

    new_all = read_caf(a.new, "rec.hdr.bnbinfo", FIELDS + NEW)
    new = unique_by_key(new_all)
    prod = unique_by_key(read_caf(a.prod, "rec.hdr.bnbinfo", FIELDS))
    py1, py2 = read_py(a.py1), read_py(a.py2)
    n = len(new["key"])
    p("# C++ beam-quality fields in the CAF: spill-by-spill check\n")
    p(f"- new CAF: `{a.new}`\n- production CAF: `{a.prod}`\n- reference: `{a.py1}` (Tier 1), `{a.py2}` (Tier 2)\n")

    # ---- 1. content
    p("## 1. New fields present and filled\n")
    missing = [f for f in NEW if f not in new]
    p(f"- bnbinfo records: {len(new_all['key'])}, unique spills: {n} (production: {len(prod['key'])})")
    p(f"- missing fields: {missing or 'none'}")
    if not missing:
        v = new["fom_version"]
        p(f"- fom_version values: {dict(zip(*np.unique(v, return_counts=True)))}")
        p(f"- FOM_best finite: {np.isfinite(new['FOM_best']).sum()} / {n}; "
          f"in [0,1]: {((new['FOM_best'] >= 0) & (new['FOM_best'] <= 1)).sum()}; "
          f"any -999 left: {(new['FOM_best'] == -999).sum()}")
        src, cnt = np.unique(new["fom_source"], return_counts=True)
        p(f"- fom_source: {dict(zip(src.tolist(), cnt.tolist()))}")
        bits = {b: int(((new['fom_status'] & b) != 0).sum()) for b in [1 << i for i in range(17)]}
        p(f"- fom_status bits set (bit: spills): { {b: c for b, c in bits.items() if c} }")
        bad = np.isfinite(new["FOM_best"]) != (new["fom_source"] >= 0)
        p(f"- FOM_best finite <=> fom_source >= 0: {'OK' if not bad.any() else f'{bad.sum()} spills disagree'}")
    p("")

    # ---- 2. same spills and inputs as production
    ia, ib = match(new, prod)
    p("## 2. Inputs vs the production CAF (same input file)\n")
    p(f"- common spills: {len(ia)} (only in new: {n - len(ia)}, only in production: {len(prod['key']) - len(ib)})")
    for f in ["TOR860", "TOR875", "HP875", "VP875", "VP873", "HPTG1", "HPTG2", "VPTG2", "M876HS", "M876VS"]:
        ok = same(new[f][ia], prod[f][ib])
        p(f"- {f}: identical for {ok.sum()} / {len(ia)}")
    for f in ["HP875Offset", "VP875Offset", "VP873Offset", "HPTG1Offset", "HPTG2Offset", "VPTG2Offset"]:
        ok = same(new[f][ia], prod[f][ib])
        filled = (prod[f][ib] == -999) & (new[f][ia] != -999)
        p(f"- {f}: identical {ok.sum()}, filled from cache (prod -999) {filled.sum()}, other {(~ok & ~filled).sum()}")
    for f in ["M875BB_spill_time_diff", "M876BB_spill_time_diff"]:
        dropped = (np.abs(prod[f][ib]) > 0.0333) & (prod[f][ib] != -999) & (new[f][ia] == -999)
        ok = same(new[f][ia], prod[f][ib], 1e-6)
        p(f"- {f}: identical {ok.sum()}, profile dropped by the 33 ms cut {dropped.sum()}, other {(~ok & ~dropped).sum()}")
    lf_new, lf_old = new["FOM"][ia], prod["FOM"][ib]
    good = lambda x: (x > 0.98) & (x <= 1)
    p(f"- legacy `FOM` field: identical {same(lf_new, lf_old, 1e-6).sum()} / {len(ia)}; "
      f"production FOM in (0.98,1]: {good(lf_old).sum()}, new FOM_best in (0.98,1]: {good(new['FOM_best'][ia]).sum()}")
    tor = np.where(prod["TOR860"][ib] > 0, prod["TOR860"][ib], 0.)
    p(f"- POT-weighted pass fraction (0.98 < x <= 1): production FOM {tor[good(lf_old)].sum() / tor.sum():.4%}, "
      f"new FOM_best {tor[good(new['FOM_best'][ia])].sum() / tor.sum():.4%}")
    p("")

    # ---- 3. C++ vs Python reference
    p("## 3. C++ (new CAF) vs the analysis-level reference\n")
    for label, ref in [("Tier 1 (bnb_fom.py)", py1), ("Tier 2 (fom_fixups.py)", py2)]:
        i, j = match(new, ref)
        filled = ((new["fom_status"][i] & FILL_BITS) != 0) | ((ref["fom_status"][j] & FILL_BITS) != 0)
        sel = ~filled if label.startswith("Tier 1") else np.ones(len(i), bool)
        i, j = i[sel], j[sel]
        p(f"### {label}: {len(i)} spills" + (" (without fill bits in either)" if label.startswith("Tier 1") else ""))
        fb = same(new["FOM_best"][i], ref["FOM_best"][j], 1e-5)
        p(f"- FOM_best agree (|d| <= 1e-5 or both NaN): {fb.sum()} / {len(i)}; "
          f"max |d| {np.nanmax(np.abs(new['FOM_best'][i] - ref['FOM_best'][j])) if len(i) else 0:.2e}")
        p(f"- fom_source agree: {(new['fom_source'][i] == ref['fom_source'][j]).sum()} / {len(i)}")
        st_new = new["fom_status"][i] & ~CPP_UNSET_BITS
        st_ref = ref["fom_status"][j].astype(np.int64) & ~CPP_UNSET_BITS
        p(f"- fom_status agree (ignoring bits 32, 4096, 8192): {(st_new == st_ref).sum()} / {len(i)}")
        for f in STATE:
            x = np.where(new[f][i] == -999, np.nan, new[f][i])
            tol = 1e-3 * max(1.0, np.nanmax(np.abs(x))) if f == "tor_used" and len(i) else 1e-3
            p(f"- {f} agree (tol {tol:.0e}): {same(x, ref[f][j], tol).sum()} / {len(i)}")
        dis = np.where(~fb | (st_new != st_ref))[0][:15]
        if len(dis):
            p("\n| spill_time | C++ FOM_best | ref FOM_best | C++ src | ref src | C++ status | ref status |\n|---|---|---|---|---|---|---|")
            for k in dis:
                p(f"| {new['spill_time_sec'][i[k]]}.{new['spill_time_nsec'][i[k]]:09d} | {new['FOM_best'][i[k]]:.6f} | "
                  f"{ref['FOM_best'][j[k]]:.6f} | {new['fom_source'][i[k]]} | {ref['fom_source'][j[k]]} | "
                  f"{new['fom_status'][i[k]]} | {ref['fom_status'][j[k]]} |")
        p("")

    # ---- 4. triggered spill vs all spills
    sp = read_caf(a.new, "rec.hdr.spillbnbinfo", ["spill_time_sec", "spill_time_nsec"] + NEW)
    i, j = match(sp, new)
    agree = all(same(sp[f][i], new[f][j]).all() for f in NEW)
    p("## 4. rec.hdr.spillbnbinfo vs rec.hdr.bnbinfo\n")
    p(f"- triggered spills: {len(sp['key'])}, found in bnbinfo: {len(i)}; all new fields identical: {agree}\n")

    # ---- 5. backward compatibility
    if a.oldfmt:
        o = read_caf(a.oldfmt, "rec.hdr.bnbinfo", ["spill_time_sec", "spill_time_nsec", "FOM"] + NEW)
        p("## 5. New CAFMaker on an old-format BNBSpillInfo (ClassVersion 15)\n")
        p(f"- spills: {len(o['key'])}; fom_version: {dict(zip(*np.unique(o['fom_version'], return_counts=True)))}; "
          f"FOM_best all NaN: {np.isnan(o['FOM_best']).all()}; fom_source all -1: {(o['fom_source'] == -1).all()}; "
          f"tgt fields all NaN: {all(np.isnan(o[f]).all() for f in STATE)}")
        oo = unique_by_key(o)
        i, j = match(oo, prod)
        p(f"- legacy FOM identical to production: {same(oo['FOM'][i], prod['FOM'][j], 1e-6).sum()} / {len(i)}\n")

    open(a.output, "w").write("\n".join(R) + "\n")
    print("\n".join(R))


if __name__ == "__main__":
    main()
