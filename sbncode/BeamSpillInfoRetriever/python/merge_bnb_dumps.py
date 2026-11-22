#!/usr/bin/env python3
"""
merge_bnb_dumps.py -- merge dump_bnb_quality.py outputs (e.g. the --chunk-files pieces)
into one file, streaming, so memory stays small however big the inputs are.

  ./merge_bnb_dumps.py -o run2_all.root run2_chunk*.root
  ./merge_bnb_dumps.py -o run2_all.root -l run2_chunks.txt

Use this instead of hadd: hadd cannot fast-clone the uproot-written jagged multiwire
branches ("TObjArray::At index out of bounds"), and it would not renumber file_idx or
join the file_list strings.  Here file_idx is offset so it keeps pointing into the merged
file_list, and the trees are rewritten in the same layout dump_bnb_quality.py uses.
"""
import argparse
import os
import sys

import awkward as ak
import numpy as np
import uproot

TREES = ("spill", "bnbinfo")


def branch_types(t):
    """{branch: dtype or 'var * int32'}, without the counter branches uproot adds for
    jagged arrays (n<name>), which mktree recreates by itself."""
    jag = {k for k in t.keys() if "var" in str(t[k].typename) or "[]" in t[k].typename}
    out = {}
    for k in t.keys():
        if k.startswith("n") and k[1:] in jag:
            continue
        out[k] = "var * int32" if k in jag else t[k].interpretation.numpy_dtype
    return out, jag


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("inputs", nargs="*", help="dump files to merge, in order")
    p.add_argument("-l", "--list", help="text file listing the inputs (e.g. <out>_chunks.txt)")
    p.add_argument("-o", "--output", required=True)
    p.add_argument("--step", type=int, default=200000, help="entries per read/write batch")
    a = p.parse_args(argv)
    files = list(a.inputs)
    if a.list:
        files += [l.strip() for l in open(a.list) if l.strip()]
    if not files:
        p.error("no inputs")
    if os.path.abspath(a.output) in map(os.path.abspath, files):
        p.error("output is also an input")

    # file_list of each input -> offsets for file_idx
    lists = []
    for f in files:
        with uproot.open(f) as fh:
            s = fh["file_list"] if "file_list" in fh else ""
            lists.append([x for x in str(s).split("\n") if x])
    offsets = np.concatenate([[0], np.cumsum([len(x) for x in lists])[:-1]])

    with uproot.open(files[0]) as fh:
        types = {tr: branch_types(fh[tr]) for tr in TREES if tr in fh}

    tmp = a.output + ".partial.root"
    with uproot.recreate(tmp) as out:
        for tr, (ty, jag) in types.items():
            wt = out.mktree(tr, ty)
            n = 0
            for i, (f, off) in enumerate(zip(files, offsets)):
                with uproot.open(f) as fh:
                    t = fh[tr]
                    tyi, _ = branch_types(t)
                    if tyi.keys() != ty.keys():
                        sys.exit(f"ERROR: {f}:{tr} has different branches from {files[0]}")
                    for b in t.iterate(list(ty), step_size=a.step, library="ak"):
                        cols = {}
                        for k in ty:
                            v = b[k]
                            cols[k] = ak.values_astype(v, np.int32) if k in jag else ak.to_numpy(v).astype(ty[k])
                        if "file_idx" in cols:
                            cols["file_idx"] = (cols["file_idx"] + off).astype(np.int32)
                        wt.extend(cols)
                        n += len(b)
                print(f"\r[{tr}] {i + 1}/{len(files)} files, {n} entries", end="", file=sys.stderr, flush=True)
            print(file=sys.stderr)
        out["file_list"] = "\n".join(x for l in lists for x in l)
    os.replace(tmp, a.output)
    print(f"wrote {a.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
