#!/usr/bin/env python3
"""
update_flatcafs.py -- write the flat CAFs of a SAM definition (or a file list) with the
BNB beam-quality result of the full procedure (Tier 1 + Tier 2).

  1. dump the beam information of all the CAFs and merge it      [dump_bnb_quality.py]
  2. Tier 1 per spill, then Tier 2 over the whole set           [bnb_fom_fast.py, fom_fixups.py]
     (the Tier-2 fills need every spill of the period in time order, so all the files
     are processed before any is written)
  3. copy each flat CAF to the output directory with, for every spill of
     rec.hdr.bnbinfo and for rec.hdr.spillbnbinfo:
         FOM_best fom_source fom_status tor_used hpos_tgt hang_tgt vpos_tgt vang_tgt
         tgt_sigx tgt_sigy fom_version
     set from step 2 (added if the CAF predates these fields). Nothing else in the file
     changes; a TNamed "bnb_beam_quality" records how the file was updated.

Missing values are NaN (FOM_best, positions, angles, widths: nominal width), fom_source -1.

Examples
  update_flatcafs.py -d my_sbnd_flatcaf_def -e sbnd --detector sbnd -o updated/ -j 16
  update_flatcafs.py -d my_icarus_flatcaf_def -e icarus --detector icarus -o updated/ -j 16
  update_flatcafs.py -f "caf_dir/*.flat.caf.root" --detector sbnd -o updated/

Large definitions: --chunk-files N dumps N files at a time (resumable; see
dump_bnb_quality.py). Step 3 reads each input again and needs PyROOT; steps 1-2 need
numpy/pandas/awkward/uproot/numba (requirements.txt).
"""
import argparse
import concurrent.futures as cf
import datetime
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

# beam-quality fields of caf::SRBNBInfo: name, ROOT leaf type code, numpy type
FIELDS = [("FOM_best", "F", np.float32), ("fom_source", "I", np.int32), ("fom_status", "i", np.uint32),
          ("tor_used", "F", np.float32), ("hpos_tgt", "F", np.float32), ("hang_tgt", "F", np.float32),
          ("vpos_tgt", "F", np.float32), ("vang_tgt", "F", np.float32), ("tgt_sigx", "F", np.float32),
          ("tgt_sigy", "F", np.float32), ("fom_version", "i", np.uint32)]
TABLE_COLS = [f for f, _, _ in FIELDS if f != "fom_version"]
FOM_VERSION = 1          # same algorithm as sbn::FOMVersion (getFOM.h), with the Tier-2 recovery
FILLED = 128 | 8192 | 16384
DETECTOR_OPTIONS = {"sbnd": [], "icarus": ["--nb-min-fom", "0.998", "--burst-fill"]}
BNB, SPILL = "rec.hdr.bnbinfo.", "rec.hdr.spillbnbinfo."


# --------------------------------------------------------------------------- steps 1-2
def beam_quality_table(a):
    """Dump, Tier 1, Tier 2; returns (path of the Tier-2 file, list of input files)."""
    import dump_bnb_quality, merge_bnb_dumps, bnb_fom, bnb_fom_fast, fom_fixups
    os.makedirs(a.work, exist_ok=True)
    dump, merged = os.path.join(a.work, "dump.root"), os.path.join(a.work, "dump_all.root")
    fixed, fx = os.path.join(a.work, "tier1.root"), os.path.join(a.work, "tier2.root")
    filelist = os.path.join(a.work, "files.txt")
    src = (["-d", a.defname] if a.defname else ["-f", *a.files] if a.files else ["-l", a.filelist])
    d = src + ["-e", a.experiment, "-o", dump, "-j", str(a.jobs), "--save-filelist", filelist]
    if a.max_files:
        d += ["-n", str(a.max_files)]
    if a.chunk_files:
        d += ["--chunk-files", str(a.chunk_files)]
    print(f"\n##### [1/3] dump -> {dump}", file=sys.stderr)
    if not (os.path.exists(dump) or os.path.exists(merged)):
        dump_bnb_quality.main(d)
    if a.chunk_files:
        if not os.path.exists(merged):
            merge_bnb_dumps.main(["-l", os.path.splitext(dump)[0] + "_chunks.txt", "-o", merged])
        dump = merged
    print(f"\n##### [2/3] Tier 1 -> {fixed}, Tier 2 -> {fx}", file=sys.stderr)
    if not os.path.exists(fixed):
        bnb_fom_fast.install()               # numba kernels in place of bnb_fom's numpy ones
        bnb_fom.main([dump, "-o", fixed])
    if not os.path.exists(fx):
        opts = DETECTOR_OPTIONS[a.detector] + (["--logger-fill"] if a.logger_fill else [])
        fom_fixups.main([fixed, "-o", fx, "--groups", str(a.groups)] + opts)
    with open(filelist) as fh:
        files = [l.strip() for l in fh if l.strip() and not l.startswith("#")]
    return fx, files


def load_table(path):
    """Per-spill result, keyed by spill time [ns]; the bnbinfo tree first (all spills),
    then the triggered spills of the spill tree (blinded CAFs have no bnbinfo)."""
    import uproot
    parts = []
    with uproot.open(path) as f:
        for tree in ("bnbinfo", "spill"):
            if tree not in f or f[tree].num_entries == 0:
                continue
            d = f[tree].arrays(["run", "spill_time_sec", "spill_time_nsec"] + TABLE_COLS, library="np")
            parts.append(d)
    T = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
    key = T["spill_time_sec"].astype(np.int64) * 1_000_000_000 + T["spill_time_nsec"].astype(np.int64)
    key, first = np.unique(key, return_index=True)          # first occurrence: bnbinfo before spill
    T = {k: v[first] for k, v in T.items()}
    T["key"] = key
    # position/angle of a plane that could not be formed (and was not filled): NaN
    st = np.nan_to_num(T["fom_status"], nan=0).astype(np.int64)
    filled = (st & FILLED) != 0
    for bit, cols in ((2, ("hpos_tgt", "hang_tgt")), (4, ("vpos_tgt", "vang_tgt"))):
        bad = ((st & bit) != 0) & ~filled
        for c in cols:
            T[c] = np.where(bad, np.nan, T[c])
    T["fom_source"] = np.nan_to_num(T["fom_source"], nan=-1).astype(np.int32)
    T["fom_status"] = st.astype(np.uint32)
    return T


def lookup(T, run, sec, nsec):
    """Values for spills (arrays); fom_version 0 and empty values where not found."""
    key = sec.astype(np.int64) * 1_000_000_000 + nsec.astype(np.int64)
    i = np.clip(np.searchsorted(T["key"], key), 0, len(T["key"]) - 1)
    found = (T["key"][i] == key) & (T["run"][i] == run)
    out = {}
    for name, _, dt in FIELDS:
        if name == "fom_version":
            out[name] = np.where(found, FOM_VERSION, 0).astype(dt)
            continue
        missing = -1 if name == "fom_source" else 0 if name == "fom_status" else np.nan
        out[name] = np.where(found, T[name][i], missing).astype(dt)
    return out, int((~found).sum())


# --------------------------------------------------------------------------- step 3
def _copy_dir(src, dst):
    """Copy every object of a ROOT directory except recTree (highest cycle only)."""
    import ROOT
    done = set()
    for k in src.GetListOfKeys():
        name = k.GetName()
        if name in done or (src.GetName() == src.GetFile().GetName() and name == "recTree"):
            continue
        done.add(name)
        obj = k.ReadObj()
        if obj.InheritsFrom("TDirectory"):
            sub = dst.mkdir(name, obj.GetTitle())
            _copy_dir(obj, sub)
        elif obj.InheritsFrom("TTree"):
            dst.cd()
            obj.CloneTree(-1, "fast").Write(name)
        else:
            dst.cd()
            obj.Write(name)


_TABLE = None


def _init_worker(table_path):
    global _TABLE
    _TABLE = load_table(table_path)          # once per worker process


def update_file(src, dst, provenance):
    """Copy one flat CAF with the beam-quality fields set from the table."""
    import uproot
    T = _TABLE
    t = uproot.open(src)["recTree"]
    cols = ["rec.hdr.run", BNB + "spill_time_sec", BNB + "spill_time_nsec",
            SPILL + "spill_time_sec", SPILL + "spill_time_nsec"]
    a = t.arrays(cols, library="np")
    nb = np.array([len(x) for x in a[BNB + "spill_time_sec"]], dtype=np.int64)
    run_b = np.repeat(a["rec.hdr.run"], nb)
    sec_b = np.concatenate(a[BNB + "spill_time_sec"]) if nb.sum() else np.zeros(0, np.uint64)
    nsec_b = np.concatenate(a[BNB + "spill_time_nsec"]) if nb.sum() else np.zeros(0, np.uint64)
    vb, miss_b = lookup(T, run_b, sec_b, nsec_b)
    vs, miss_s = lookup(T, a["rec.hdr.run"], a[SPILL + "spill_time_sec"], a[SPILL + "spill_time_nsec"])
    # a record without a triggered spill (time 0) is not "missing"
    miss_s -= int(((a[SPILL + "spill_time_sec"] == 0) | (a[SPILL + "spill_time_sec"] > 4e9)).sum())
    off = np.concatenate([[0], np.cumsum(nb)])

    import ROOT
    fin = ROOT.TFile.Open(src)
    tin = fin.Get("recTree")
    tmp = dst + ".partial"
    fout = ROOT.TFile(tmp, "RECREATE", "", fin.GetCompressionSettings())
    _copy_dir(fin, fout)
    fout.cd()
    tout = tin.CloneTree(0)
    nmax = max(1, int(nb.max()) if len(nb) else 1)
    bufs, sbufs = {}, {}
    for name, code, dt in FIELDS:
        bufs[name] = np.zeros(nmax, dt)
        sbufs[name] = np.zeros(1, dt)
        for br, buf, leaf in ((BNB + name, bufs[name], f"{BNB}{name}[{BNB}.length]/{code}"),
                              (SPILL + name, sbufs[name], f"{SPILL}{name}/{code}")):
            if tout.GetBranch(br):
                tout.SetBranchAddress(br, buf)
            else:
                tout.Branch(br, buf, leaf)
    for e in range(tin.GetEntries()):
        tin.GetEntry(e)
        lo, hi = off[e], off[e + 1]
        for name, _, _ in FIELDS:
            bufs[name][:hi - lo] = vb[name][lo:hi]
            sbufs[name][0] = vs[name][e]
        tout.Fill()
    tout.Write("", ROOT.TObject.kOverwrite)
    prov = dict(provenance, input=src, spills=int(nb.sum()), spills_without_result=miss_b,
                triggered_spills_without_result=max(miss_s, 0))
    ROOT.TNamed("bnb_beam_quality", json.dumps(prov)).Write()
    fout.Close()
    fin.Close()
    os.replace(tmp, dst)
    return os.path.basename(dst), int(nb.sum()), miss_b, max(miss_s, 0)


def _local(path):
    """root://.../pnfs/... -> /pnfs/... when /pnfs is mounted (as the dump reads it)."""
    import dump_bnb_quality
    p = dump_bnb_quality.pnfs_local(path)
    return p if p and os.path.exists(p) else path


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("-d", "--defname", help="SAM definition of flat CAFs")
    src.add_argument("-f", "--files", nargs="+", help="flat CAF files/URLs (globs ok)")
    src.add_argument("-l", "--filelist", help="text file with one flat CAF path/URL per line")
    p.add_argument("-e", "--experiment", default=os.environ.get("SAM_EXPERIMENT", "sbnd"))
    p.add_argument("--detector", required=True, choices=sorted(DETECTOR_OPTIONS),
                   help="Tier-2 settings: icarus = neighbor fill only between spills with FOM > 0.998, "
                        "875-dropout fill on; sbnd = no neighbor-FOM requirement, no 875-dropout fill")
    p.add_argument("--logger-fill", action="store_true",
                   help="ICARUS Run 2 only: fill the June 2023 outage from the ACNet logger (acnet/)")
    p.add_argument("-o", "--outdir", required=True, help="directory for the updated flat CAFs")
    p.add_argument("--work", help="directory for the intermediate files (default <outdir>/bnb_work); "
                                  "rerunning reuses the finished steps")
    p.add_argument("-j", "--jobs", type=int, default=1, help="parallel readers / writers")
    p.add_argument("-n", "--max-files", type=int, default=0, help="only the first N files")
    p.add_argument("--chunk-files", type=int, default=0, metavar="N", help="dump N files at a time (large definitions)")
    p.add_argument("--groups", type=int, default=1, help="Tier 2 in this many groups of runs (memory)")
    p.add_argument("--table-version", default="", help="label recorded in each output file")
    a = p.parse_args(argv)
    a.work = a.work or os.path.join(a.outdir, "bnb_work")
    os.makedirs(a.outdir, exist_ok=True)

    table, files = beam_quality_table(a)
    prov = {"tool": "sbncode BeamSpillInfoRetriever/python/update_flatcafs.py",
            "sbncode": os.environ.get("SBNCODE_VERSION", "unknown"), "fom_version": FOM_VERSION,
            "detector": a.detector, "logger_fill": a.logger_fill, "table_version": a.table_version,
            "definition": a.defname or "", "date": datetime.datetime.utcnow().isoformat(timespec="seconds") + "Z"}
    print(f"\n##### [3/3] {len(files)} flat CAFs -> {a.outdir}/", file=sys.stderr)
    todo = [(_local(f), os.path.join(a.outdir, os.path.basename(f.split("?")[0]))) for f in files]
    todo = [(s, d) for s, d in todo if not os.path.exists(d)]
    nmiss = 0
    with cf.ProcessPoolExecutor(max_workers=max(1, min(a.jobs, len(todo) or 1)),
                                initializer=_init_worker, initargs=(table,)) as ex:
        for name, n, mb, ms in ex.map(update_file, [s for s, _ in todo], [d for _, d in todo],
                                      [prov] * len(todo)):
            nmiss += mb + ms
            print(f"[update] {name}: {n} spills" + (f", {mb} + {ms} triggered without a result" if mb or ms else ""),
                  file=sys.stderr)
    print(f"done: {len(todo)} files written to {a.outdir}" + (f" ({nmiss} spills without a result)" if nmiss else ""),
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
