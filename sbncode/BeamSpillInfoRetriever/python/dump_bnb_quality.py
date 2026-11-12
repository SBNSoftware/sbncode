#!/usr/bin/env python3
"""
dump_bnb_quality.py -- pull BNB beam-quality info out of SBN flat CAFs.

Given a SAM definition (or explicit files), this reads
    rec.hdr.spillbnbinfo  -> one row per CAF record (the spill matched to the trigger)
    rec.hdr.bnbinfo       -> one row per spill in the subrun (all spills, stored only
                             on the first record of each subrun)
and writes both as flat tables for quick study.

Output trees / tables
---------------------
  spill    : run, subrun, evt, spill_time, all SRBNBInfo scalars, M875BB/M876BB/MMBTBB
             (multiwire bunch-by-bunch arrays, may be empty), source_file, has_spill_info
  bnbinfo  : run, subrun, spill_time, all SRBNBInfo scalars + multiwire arrays,
             subrun_pot (rec.hdr.pot), triggered (spill matches a record in 'spill'),
             source_file

Examples
--------
  # SAM definition, read over xrootd (default)
  ./dump_bnb_quality.py -d my_bnb_caf_def -o beamqual.root

  # first 20 files, 4 parallel readers, parquet output
  ./dump_bnb_quality.py -d my_bnb_caf_def -n 20 -j 4 -o beamqual.parquet

  # files on disk (globs ok) or a text file list
  ./dump_bnb_quality.py -f Run2_Files/*.root Run4_Files/*.root -o test.root
  ./dump_bnb_quality.py -l filelist.txt -o test.csv

  # big definition: 5000-file chunks -> out_chunk000.root, out_chunk001.root, ...
  # (rerun the same command after a failure: finished chunks are skipped and a
  #  half-done chunk resumes from its completed worker parts)
  ./dump_bnb_quality.py -d big_def -e sbn -j 32 --chunk-files 5000 -o <output dir>/out.root

Requires: uproot, awkward, numpy, pandas (+ pyarrow for parquet).
For xrootd URLs uproot needs the XRootD python bindings or fsspec-xrootd
(both available in the usual SBN spack/ups environments). SAM lookup uses the
`samweb` CLI; set SAM_EXPERIMENT (icarus / sbnd) or pass -e.
"""

import argparse
import pickle
import concurrent.futures as cf
import glob
import os
import shutil
import subprocess
import sys

import awkward as ak
import numpy as np
import pandas as pd
import uproot

TREE = "recTree"
SPILL = "rec.hdr.spillbnbinfo."
SUBRUN = "rec.hdr.bnbinfo."
MULTIWIRE = ("M875BB", "M876BB", "MMBTBB")
HDR = ["rec.hdr.run", "rec.hdr.subrun", "rec.hdr.evt", "rec.hdr.first_in_subrun", "rec.hdr.pot"]


# --------------------------------------------------------------------------- SAM
def samweb(args, experiment):
    if shutil.which("samweb") is None:
        sys.exit("ERROR: samweb not found -- set up sam_web_client (or use -f/-l).")
    cmd = ["samweb"] + (["-e", experiment] if experiment else []) + args
    out = subprocess.run(cmd, capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError(f"{' '.join(cmd)} failed:\n{out.stderr.strip()}")
    return [l.strip() for l in out.stdout.splitlines() if l.strip()]


def _loc_to_input(name, loc, access):
    """SAM location ('dcache:/pnfs/...', 'enstore:/pnfs/...(tape label)') -> path/URL for uproot."""
    d = loc.split(":", 1)[-1].split("(")[0].rstrip("/")
    if "/pnfs/" not in d and not d.startswith("/pnfs"):
        return None
    path = f"{d}/{name}"
    if access == "local" or os.path.exists(path):
        return path                                   # /pnfs mounted (GPVM): read directly
    rel = path.split("/pnfs/", 1)[1]
    if not rel.startswith("fnal.gov/usr/"):
        rel = "fnal.gov/usr/" + rel
    return f"root://fndcadoor.fnal.gov:1094/pnfs/{rel}"


def _locate_bulk(names, experiment, access):
    """Bulk SAM location lookup through the samweb_client python API (a few hundred files per
    request instead of one CLI call per file). Returns {name: input} or None if unavailable."""
    try:
        import samweb_client
        client = samweb_client.SAMWebClient(experiment=experiment)
        if not hasattr(client, "locateFiles"):
            return None
    except Exception:
        return None
    out = {}
    for i in range(0, len(names), 500):
        batch = names[i:i + 500]
        try:
            res = client.locateFiles(batch)
        except Exception as e:
            print(f"\n[sam] bulk locate failed ({type(e).__name__}: {e}); falling back to per-file calls",
                  file=sys.stderr)
            return None
        for name in batch:
            for loc in (res.get(name) or []):
                full = loc.get("full_path") or loc.get("location") or "" if isinstance(loc, dict) else str(loc)
                inp = _loc_to_input(name, full, access)
                if inp:
                    out[name] = inp
                    break
        print(f"\r[sam] located {min(i + 500, len(names))}/{len(names)} files (bulk)", end="",
              file=sys.stderr, flush=True)
    print(file=sys.stderr)
    return out


def resolve_sam(defname, experiment, access, nmax, nproc):
    names = samweb(["list-files", f"defname:{defname}"], experiment)
    names.sort()
    print(f"[sam] {defname}: {len(names)} files", file=sys.stderr)
    if nmax:
        names = names[:nmax]

    if access in ("root", "local"):
        bulk = _locate_bulk(names, experiment, access)
        if bulk is not None:
            missing = [n for n in names if n not in bulk]
            if missing:
                print(f"[sam] {len(missing)} files without a /pnfs location; resolving those one by one",
                      file=sys.stderr)
            if not missing:
                return [bulk[n] for n in names]
        else:
            print("[sam] samweb_client python API not available: one samweb call per file "
                  "(slow for large definitions; use --save-filelist to do this only once)", file=sys.stderr)
            missing, bulk = names, {}
    else:
        missing, bulk = names, {}

    def one(name):
        if access == "local":
            for loc in samweb(["locate-file", name], experiment):
                inp = _loc_to_input(name, loc, "local")
                if inp and os.path.exists(inp):
                    return inp
            raise RuntimeError(f"no locally mounted location for {name}")
        urls = samweb(["get-file-access-url", "--schema=" + access, name], experiment)
        if not urls:
            raise RuntimeError(f"no {access} url for {name}")
        return urls[0]

    got = dict(bulk)
    with cf.ThreadPoolExecutor(max(1, min(32, nproc * 4))) as ex:
        futs = {ex.submit(one, n): n for n in missing}
        for i, fut in enumerate(cf.as_completed(futs), 1):
            name = futs[fut]
            try:
                got[name] = fut.result()
            except Exception as e:
                print(f"\n[sam] WARNING skipping {name}: {e}", file=sys.stderr)
            if i % 50 == 0 or i == len(missing):
                print(f"\r[sam] resolved {i}/{len(missing)} files", end="", file=sys.stderr, flush=True)
    print(file=sys.stderr)
    return [got[n] for n in names if n in got]


# --------------------------------------------------------------------------- reading
def _strip(name, prefix):
    return name[len(prefix):]


# Beam-quality result written by the production retriever (sbncode getFOM.h, CAFs
# with fom_version >= 1). Stored as cpp_<name>, so that it is kept next to the
# columns of the same name that bnb_fom.py / fom_fixups.py compute.
CPP_RESULT = ("FOM_best", "fom_source", "fom_status", "tor_used", "hpos_tgt", "hang_tgt",
              "vpos_tgt", "vang_tgt", "tgt_sigx", "tgt_sigy", "fom_version")


def _column(name, prefix):
    col = _strip(name, prefix)
    return "cpp_" + col if col in CPP_RESULT else col


def _scalar_branches(tree, prefix):
    """SRBNBInfo leaves, minus multiwire arrays and ROOT's '..' helper branches."""
    return [k for k in tree.keys() if k.startswith(prefix) and ".." not in k
            and _strip(k, prefix) not in MULTIWIRE]


def pnfs_local(url):
    """root://fndcadoor.fnal.gov:1094/pnfs/fnal.gov/usr/<exp>/... -> /pnfs/<exp>/... if mounted."""
    if "/pnfs/" not in url:
        return None
    path = "/pnfs/" + url.split("/pnfs/", 1)[1]
    path = path.replace("/pnfs/fnal.gov/usr/", "/pnfs/", 1)
    return path if os.path.exists(path) else None


def _have(mod):
    try:
        __import__(mod)
        return True
    except ImportError:
        return False


def open_tree(path):
    """Open recTree, picking whatever xrootd access actually works in this environment:
       1. local /pnfs mount (GPVMs)          2. fsspec-xrootd (uproot default)
       3. XRootD python bindings (pyxrootd, shipped with ROOT/xrootd in SBN envs)"""
    if not path.startswith(("root://", "xroot://")):
        return uproot.open(path, timeout=120)[TREE]
    if XROOT_MODE in ("auto", "pnfs"):
        local = pnfs_local(path)
        if local:
            return uproot.open(local)[TREE]
        if XROOT_MODE == "pnfs":
            raise RuntimeError(f"/pnfs path for {path} not mounted here")
    if _have("fsspec_xrootd"):
        return uproot.open(path, timeout=120)[TREE]
    if _have("XRootD"):
        return uproot.open(path, timeout=120, handler=uproot.source.xrootd.XRootDSource)[TREE]
    raise ImportError("no xrootd access: /pnfs not mounted, and neither fsspec-xrootd nor the "
                      "XRootD python bindings are importable (setup xrootd / ROOT, or "
                      "pip install --user fsspec-xrootd)")


XROOT_MODE = "auto"


def _set_mode(mode):
    global XROOT_MODE
    XROOT_MODE = mode


def _empty_lists(n):
    return ak.unflatten(np.zeros(0, np.int32), np.zeros(n, np.int64))


def read_file(path, read_mw=True):
    t = open_tree(path)        # ~90% of the time per file: uproot parses all ~2000 branch headers
    src = os.path.basename(path.split("?")[0])
    mw = MULTIWIRE if read_mw else ()

    # ---------------- per-record (triggered spill) info
    sbr = _scalar_branches(t, SPILL)
    a = t.arrays(HDR + sbr + [SPILL + m for m in mw], library="ak")
    spill = {"run": ak.to_numpy(a["rec.hdr.run"]),
             "subrun": ak.to_numpy(a["rec.hdr.subrun"]),
             "evt": ak.to_numpy(a["rec.hdr.evt"])}
    for b in sbr:
        spill[_column(b, SPILL)] = ak.to_numpy(a[b])
    mw_spill = {m: (a[SPILL + m] if read_mw else _empty_lists(len(a))) for m in MULTIWIRE}

    # ---------------- per-subrun spill list (only filled on first_in_subrun records)
    bbr = _scalar_branches(t, SUBRUN)
    b = t.arrays(bbr + [SUBRUN + m for m in mw] + [SUBRUN + m + "..length" for m in mw], library="ak")
    nsp = ak.num(b[bbr[0]])
    sub = {"run": np.repeat(spill["run"], nsp),
           "subrun": np.repeat(spill["subrun"], nsp),
           "subrun_pot": np.repeat(ak.to_numpy(a["rec.hdr.pot"]), nsp)}
    for br in bbr:
        sub[_column(br, SUBRUN)] = ak.to_numpy(ak.flatten(b[br]))
    # multiwire is stored flat per record; split it back into one array per spill
    mw_sub = {m: (ak.unflatten(ak.flatten(b[SUBRUN + m]), ak.flatten(b[SUBRUN + m + "..length"]))
                  if read_mw else _empty_lists(int(ak.sum(nsp))))
              for m in MULTIWIRE}

    return src, spill, mw_spill, sub, mw_sub


def to_frames(results):
    sp_df, sp_mw, sb_df, sb_mw = [], {m: [] for m in MULTIWIRE}, [], {m: [] for m in MULTIWIRE}
    for src, spill, mw_spill, sub, mw_sub in results:
        d = pd.DataFrame(spill); d["source_file"] = src; sp_df.append(d)
        d = pd.DataFrame(sub);   d["source_file"] = src; sb_df.append(d)
        for m in MULTIWIRE:
            sp_mw[m].append(mw_spill[m]); sb_mw[m].append(mw_sub[m])
    sp = pd.concat(sp_df, ignore_index=True)
    sb = pd.concat(sb_df, ignore_index=True)
    sp_mw = {m: ak.concatenate(v) for m, v in sp_mw.items()}
    sb_mw = {m: ak.concatenate(v) for m, v in sb_mw.items()}

    for df in (sp, sb):
        df["spill_time"] = df["spill_time_sec"].astype("f8") + 1e-9 * df["spill_time_nsec"].astype("f8")
    # CAFMaker leaves spillbnbinfo at its default if the event wasn't in the spill map
    sp["has_spill_info"] = ((sp["spill_time_sec"] > 0)
                            & (sp["spill_time_sec"] < np.iinfo(np.uint64).max)
                            & (sp["event"] != np.iinfo(np.uint32).max))

    # flag subrun spills that ended up as a triggered record
    key = ["run", "subrun", "spill_time_sec", "spill_time_nsec"]
    trig = sp.loc[sp.has_spill_info, key].drop_duplicates()
    trig["triggered"] = True
    sb = sb.merge(trig, on=key, how="left")      # left merge on unique keys keeps row order
    sb["triggered"] = sb["triggered"].eq(True)

    # tidy column order
    lead_sp = ["run", "subrun", "evt", "spill_time", "has_spill_info"]
    lead_sb = ["run", "subrun", "spill_time", "triggered", "subrun_pot"]
    sp = sp[lead_sp + [c for c in sp.columns if c not in lead_sp]]
    sb = sb[lead_sb + [c for c in sb.columns if c not in lead_sb]]
    return sp, sp_mw, sb, sb_mw


# --------------------------------------------------------------------------- writing
def write(out, sp, sp_mw, sb, sb_mw, keep_mw):
    ext = os.path.splitext(out)[1].lower()
    if ext == ".root":
        files = sorted(set(sp.source_file))
        fidx = {n: i for i, n in enumerate(files)}
        with uproot.recreate(out) as f:
            for name, df, mw in (("spill", sp, sp_mw), ("bnbinfo", sb, sb_mw)):
                cols = {c: df[c].to_numpy() for c in df.columns if c != "source_file"}
                cols["file_idx"] = df["source_file"].map(fidx).to_numpy(np.int32)
                for c in ("has_spill_info", "triggered"):
                    if c in cols:
                        cols[c] = cols[c].astype(np.int8)
                if keep_mw:
                    for m in MULTIWIRE:
                        cols[m] = ak.values_astype(mw[m], np.int32)
                write_ttree(f, name, cols)
            # string branches are awkward in ROOT; store file_idx -> name as a TObjString
            f["file_list"] = "\n".join(files)
    elif ext in (".parquet", ".pq"):
        base = os.path.splitext(out)[0]
        for name, df, mw in (("spill", sp, sp_mw), ("bnbinfo", sb, sb_mw)):
            df = df.copy()
            if keep_mw:
                for m in MULTIWIRE:
                    df[m] = ak.to_list(mw[m])
            df.to_parquet(f"{base}_{name}.parquet", index=False)
    elif ext == ".csv":
        base = os.path.splitext(out)[0]
        for name, df in (("spill", sp), ("bnbinfo", sb)):
            df.to_csv(f"{base}_{name}.csv", index=False)   # multiwire arrays not written to CSV
    elif ext in (".h5", ".hdf5"):
        sp.to_hdf(out, key="spill", mode="w")
        sb.to_hdf(out, key="bnbinfo", mode="a")
    else:
        sys.exit(f"ERROR: unknown output extension '{ext}' (use .root, .parquet, .csv, .h5)")


def write_ttree(f, name, cols, step=100000):
    """Write `cols` as a TTree (not RNTuple: recent uproot defaults to RNTuple for
    `f[name] = ...`) in chunks of `step` entries, so each branch gets many modest
    baskets instead of one huge one (which forces readers to decompress everything)."""
    types = {c: ("var * int32" if isinstance(v, ak.Array) else v.dtype) for c, v in cols.items()}
    tree = f.mktree(name, types)
    n = len(next(iter(cols.values())))
    for s in range(0, max(n, 1), step):
        tree.extend({c: v[s:s + step] for c, v in cols.items()})


def summarize(sp, sp_mw, sb):
    print("\n=== summary ===", file=sys.stderr)
    print(f"records (spill table)      : {len(sp)}  (with spill info: {sp.has_spill_info.sum()})", file=sys.stderr)
    print(f"subrun spills (bnbinfo)    : {len(sb)}  (triggered: {sb.triggered.sum()})", file=sys.stderr)
    print(f"runs                       : {sp.run.nunique()}  [{sp.run.min()}-{sp.run.max()}]", file=sys.stderr)
    pot = sb.drop_duplicates(["run", "subrun", "source_file"]).subrun_pot.sum()
    print(f"sum rec.hdr.pot            : {pot:.4e}", file=sys.stderr)
    s = sp[sp.has_spill_info]
    if len(s):
        print(f"TOR860 (triggered) median  : {np.median(s.TOR860):.3e}", file=sys.stderr)
        for lbl, m in (("FOM == -999 (not computed)", s.FOM == -999),
                       ("TOR860 <= 0", s.TOR860 <= 0),
                       ("multiwire BB empty", ak.to_numpy(ak.num(sp_mw["M875BB"]))[sp.has_spill_info.to_numpy()] == 0)):
            print(f"{lbl:<27}: {m.sum()} / {len(s)}", file=sys.stderr)


# --------------------------------------------------------------------------- main
def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = p.add_mutually_exclusive_group()
    src.add_argument("-d", "--defname", help="SAM definition name")
    src.add_argument("-f", "--files", nargs="+", help="CAF files/URLs (globs ok)")
    src.add_argument("-l", "--filelist", help="text file with one CAF path/URL per line")
    p.add_argument("-o", "--output", default="bnb_beam_quality.root",
                   help="output file; format from extension: .root (default), .parquet, .csv, .h5")
    p.add_argument("-e", "--experiment", default=os.environ.get("SAM_EXPERIMENT", "icarus"),
                   help="SAM experiment (default: $SAM_EXPERIMENT or icarus)")
    p.add_argument("-a", "--access", default="root", choices=["root", "https", "local"],
                   help="how to reach SAM files: xrootd (default), https, or local /pnfs path")
    p.add_argument("--xrootd", default="auto", choices=["auto", "pnfs", "stream"],
                   help="for root:// inputs: auto = use local /pnfs mount if present, else stream; "
                        "pnfs = require /pnfs; stream = always stream over xrootd")
    p.add_argument("-n", "--max-files", type=int, default=0, help="only process the first N files")
    p.add_argument("-j", "--jobs", type=int, default=1,
                   help="parallel worker processes (independent copies of this script, each reading "
                        "every N-th file; output identical to -j 1)")
    p.add_argument("--threads", action="store_true",
                   help="use threads in one process for -j instead (little gain: per-file cost is CPU-bound)")
    p.add_argument("--save-filelist", help="write the resolved file list/URLs here (reuse later with -l)")
    p.add_argument("--shard", help=argparse.SUPPRESS)          # internal: K/N, write raw results to -o
    p.add_argument("--workdir", help="where -j workers keep their part files "
                                     "(default: the output file's directory; /tmp is often too small)")
    src.add_argument("--resume", metavar="PARTDIR",
                     help="finish an interrupted -j run from its part-file directory: reuse the "
                          "completed parts, rerun only the failed/missing workers, then write -o")
    p.add_argument("--chunk-files", type=int, default=0, metavar="N",
                   help="process the inputs in consecutive blocks of N files, one output per block "
                        "(<out>_chunk000.root, ...; list in <out>_chunks.txt). Rerunning the same "
                        "command skips finished chunks and resumes a partly finished one")
    p.add_argument("--no-multiwire", action="store_true", help="drop M875BB/M876BB/MMBTBB arrays")
    p.add_argument("--dedup", action="store_true",
                   help="drop duplicate subrun spills (same run/subrun/spill time) from bnbinfo")
    args = p.parse_args(argv)

    global XROOT_MODE
    XROOT_MODE = args.xrootd
    if not (args.defname or args.files or args.filelist or args.resume):
        p.error("one of -d/--defname, -f/--files, -l/--filelist, --resume is required")
    if args.resume:
        results, bad = run_workers(None, args)
        if not results:
            sys.exit("ERROR: nothing read")
        finish(results, bad, args)
        return
    if args.defname:
        files = resolve_sam(args.defname, args.experiment, args.access, args.max_files, args.jobs)
    else:
        if args.files:
            files = [g for f in args.files for g in (sorted(glob.glob(f)) or [f])]
        else:
            with open(args.filelist) as fh:
                files = [l.strip() for l in fh if l.strip() and not l.startswith("#")]
        if args.max_files:
            files = files[:args.max_files]
    if not files:
        sys.exit("ERROR: no input files")
    if args.save_filelist and not args.shard:
        with open(args.save_filelist, "w") as fh:
            fh.write("\n".join(files) + "\n")
        print(f"[files] {len(files)} inputs saved to {args.save_filelist}", file=sys.stderr)
    if args.chunk_files and not args.shard:
        run_chunks(files, args)
        return

    read_mw = not args.no_multiwire
    if args.shard:                                           # ---- worker: read a slice, save raw results
        k, n = (int(x) for x in args.shard.split("/"))
        mine = files[k::n]
        results, bad = [], []
        prog = args.output + ".progress"
        for i, f in enumerate(mine, 1):
            try:
                results.append(read_file(f, read_mw))
            except Exception as e:
                bad.append(f)
                print(f"WARNING failed {f}: {type(e).__name__}: {e}", file=sys.stderr, flush=True)
            if i % 10 == 0 or i == len(mine):
                try:                      # dCache (/pnfs) refuses to overwrite: write anew each time
                    if os.path.exists(prog):
                        os.remove(prog)
                    with open(prog, "w") as fh:
                        fh.write(str(i))
                except OSError:
                    pass                  # progress display only
        with open(args.output, "wb") as fh:
            pickle.dump((results, bad), fh, protocol=pickle.HIGHEST_PROTOCOL)
        return
    if args.jobs > 1 and not args.threads:
        results, bad = run_workers(files, args)
    else:
        results, bad = read_threads(files, args, read_mw)
    if not results:
        sys.exit("ERROR: nothing read")
    finish(results, bad, args)


def run_workers(files, args, wd=None, hint=None):
    """Run N copies of this script, each on every N-th file, then collect their results.
    Separate processes (not multiprocessing): no GIL, and nothing for an ups/cvmfs PYTHONPATH
    to hang on.  Part files live in a work directory; if any worker fails the directory is
    kept and the run can be finished with --resume <dir> (only failed workers are rerun)."""
    import subprocess
    import tempfile
    import time
    if wd is None and args.resume:
        wd = args.resume
    resuming = wd is not None and os.path.exists(os.path.join(wd, "files.txt"))
    if resuming:                                              # reuse an existing part directory
        wd = os.path.abspath(wd)
        flist = os.path.join(wd, "files.txt")
        with open(flist) as fh:
            old = [l.strip() for l in fh if l.strip()]
        if files is not None and old != list(files):
            sys.exit(f"ERROR: {wd} holds parts for a different file list; remove it or restore "
                     "the original inputs/--chunk-files")
        files = old
        try:
            with open(os.path.join(wd, "njobs")) as fh:
                n = int(fh.read())
        except OSError:                                       # older runs: count the worker logs
            n = len(glob.glob(os.path.join(wd, "part*.log")))
        if n < 1:
            sys.exit(f"ERROR: {wd} does not look like a part directory")
        print(f"[resume] {wd}: {len(files)} files, {n} workers", file=sys.stderr)
    else:
        n = min(args.jobs, len(files))
        if wd is None:
            base = args.workdir or os.path.dirname(os.path.abspath(args.output))
            wd = tempfile.mkdtemp(prefix="bnbdump_", dir=base)
        else:
            os.makedirs(wd, exist_ok=True)
        flist = os.path.join(wd, "files.txt")
        with open(flist, "w") as fh:
            fh.write("\n".join(files) + "\n")
        with open(os.path.join(wd, "njobs"), "w") as fh:
            fh.write(str(n))
    part = lambda k: os.path.join(wd, f"part{k}.pkl")
    loaded = {}

    def try_load(k):
        try:
            with open(part(k), "rb") as fh:
                loaded[k] = pickle.load(fh)
            return True
        except Exception:
            return False

    todo = [k for k in range(n) if not (resuming and try_load(k))]
    if resuming:
        print(f"[resume] {n - len(todo)} parts complete, rerunning workers {todo}", file=sys.stderr)
    procs = {}
    for k in todo:
        for f in (part(k), part(k) + ".progress", os.path.join(wd, f"part{k}.log")):
            if os.path.exists(f):
                os.remove(f)
        cmd = [sys.executable, os.path.abspath(__file__), "-l", flist, "--shard", f"{k}/{n}", "-o", part(k),
               "--xrootd", args.xrootd] + (["--no-multiwire"] if args.no_multiwire else [])
        log = open(os.path.join(wd, f"part{k}.log"), "w")
        procs[k] = (subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT), log)
    print(f"[workers] {len(procs)} processes, part files in {wd}", file=sys.stderr)
    nfiles = sum(len(files[k::n]) for k in todo)
    t0 = time.time()
    while any(p.poll() is None for p, _ in procs.values()):
        done = 0
        for k in todo:
            try:
                with open(part(k) + ".progress") as fh:
                    done += int(fh.read() or 0)
            except (OSError, ValueError):
                pass
        el = time.time() - t0
        eta = (nfiles - done) * el / done if done else float("nan")
        print(f"\r[read] {done}/{nfiles} files  {el/60:.1f} min  ETA {eta/60:.0f} min ",
              end="", file=sys.stderr, flush=True)
        time.sleep(5)
    print(file=sys.stderr)
    failed = []
    for k, (p, log) in procs.items():
        log.close()
        if p.returncode != 0 or not try_load(k):
            failed.append(k)
            print(f"ERROR: worker {k} failed (exit {p.returncode}); see {os.path.join(wd, f'part{k}.log')}",
                  file=sys.stderr)
    if failed:
        import shutil
        free = shutil.disk_usage(wd).free / 1e9
        sys.exit(f"ERROR: {len(failed)} of {n} workers failed; the other parts are kept "
                 f"({free:.1f} GB free in {wd}).\n"
                 + (hint or f"Fix the cause, then finish with:\n  {sys.argv[0]} --resume {wd} -o {args.output}"
                    + (" --no-multiwire" if args.no_multiwire else "")))
    results, bad = [], []
    for k in range(n):
        r, b = loaded.pop(k)
        results += r
        bad += b
    for b in bad:
        print(f"WARNING failed {b}", file=sys.stderr)
    import shutil
    shutil.rmtree(wd, ignore_errors=True)
    return results, bad


def read_threads(files, args, read_mw):
    results, bad = [], []
    with cf.ThreadPoolExecutor(max(1, args.jobs)) as ex:
        futs = {ex.submit(read_file, f, read_mw): f for f in files}
        for i, fut in enumerate(cf.as_completed(futs), 1):
            f = futs[fut]
            try:
                results.append(fut.result())
            except Exception as e:
                bad.append(f)
                print(f"\nWARNING failed {f}: {type(e).__name__}: {e}", file=sys.stderr)
                if isinstance(e, ImportError) and not results and len(bad) >= min(5, len(files)):
                    for fu in futs:
                        fu.cancel()
                    sys.exit(f"ERROR: aborting -- {e}")
            print(f"\r[read] {i}/{len(files)} files", end="", file=sys.stderr, flush=True)
    print(file=sys.stderr)
    return results, bad


def run_chunks(files, args):
    """--chunk-files: consecutive blocks of files, one complete output each.  An output is
    written under a temporary name and renamed when done, so an existing chunk file is
    always complete and is skipped on a rerun."""
    import copy
    import math
    import time
    stem, ext = os.path.splitext(args.output)
    N = args.chunk_files
    nch = math.ceil(len(files) / N)
    outs = [f"{stem}_chunk{c:03d}{ext}" for c in range(nch)]
    with open(f"{stem}_chunks.txt", "w") as fh:
        fh.write("\n".join(outs) + "\n")
    print(f"[chunks] {len(files)} files -> {nch} chunks of <= {N}; outputs listed in {stem}_chunks.txt",
          file=sys.stderr)
    hint = "Fix the cause, then rerun the same command: finished chunks are skipped and this one resumes."
    t0 = time.time()
    for c, out in enumerate(outs):
        if os.path.exists(out):
            print(f"[chunk {c + 1}/{nch}] {out} exists, skipping", file=sys.stderr)
            continue
        part = files[c * N:(c + 1) * N]
        print(f"\n[chunk {c + 1}/{nch}] {len(part)} files -> {out}", file=sys.stderr)
        a = copy.copy(args)
        a.output = f"{stem}_chunk{c:03d}.partial{ext}"
        if args.jobs > 1 and not args.threads:
            results, bad = run_workers(part, a, wd=f"{stem}_chunk{c:03d}.parts", hint=hint)
        else:
            results, bad = read_threads(part, a, not args.no_multiwire)
        if not results:
            sys.exit(f"ERROR: chunk {c}: nothing read")
        finish(results, bad, a)
        os.replace(a.output, out)
        print(f"[chunk {c + 1}/{nch}] done -> {out}  ({(time.time() - t0) / 60:.1f} min total)",
              file=sys.stderr)
    print(f"\nall {nch} chunks done; outputs listed in {stem}_chunks.txt", file=sys.stderr)


def finish(results, bad, args):
    results.sort(key=lambda r: r[0])

    sp, sp_mw, sb, sb_mw = to_frames(results)
    if args.dedup:
        keep = ~sb.duplicated(["run", "subrun", "spill_time_sec", "spill_time_nsec"]).to_numpy()
        sb = sb[keep].reset_index(drop=True)
        sb_mw = {m: v[keep] for m, v in sb_mw.items()}

    write(args.output, sp, sp_mw, sb, sb_mw, keep_mw=not args.no_multiwire)
    summarize(sp, sp_mw, sb)
    print(f"\nwrote {args.output}" + (f"  ({len(bad)} files failed)" if bad else ""), file=sys.stderr)


if __name__ == "__main__":
    main()
