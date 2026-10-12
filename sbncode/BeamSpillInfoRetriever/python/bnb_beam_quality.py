#!/usr/bin/env python3
"""
bnb_beam_quality.py -- BNB beam-quality / FOM toolkit for SBN flat CAFs.

One entry point for the steps of the procedure (each is also usable on its own):

  dump      CAFs (SAM definition, files, file list) -> beamqual.root
            rec.hdr.spillbnbinfo (per record) and rec.hdr.bnbinfo (all spills
            in the subrun) as flat TTrees            [dump_bnb_quality.py]
  fom       beamqual.root -> beamqual_fom.root
            recomputes the FOM twice: "legacy" (exact copy of production,
            closes against the CAF) and "fixed" (bugs corrected)  [bnb_fom.py]
  plot      *_fom.root -> plots + summary.md, old vs new, with and without the
            +100 nominal-width convention                 [plot_fom_compare.py]
  all       dump -> fom -> plot in one go
  fixups    *_fom.root -> Tier-1 corrections and Tier-2 recovery steps
            (target BPM, TOR875 back-up, widths, neighbor / logger / 875-dropout
            fills)                                          [fom_fixups.py]
  validate  validation report of a *_fom.root file          [validate_fom.py]
  compare-cpp  spill-by-spill check of the production (C++) result in a CAF
            against this code                               [compare_cpp_caf.py]
  selftest  run dump + fom on a few CAF files (exit code 0 = pass): the Python
            Tier 1 must reproduce the C++ FOM_best of CAFs that carry it, or the
            legacy FOM of older CAFs

Examples
  ./bnb_beam_quality.py all -d my_icarus_bnb_caf_def -o run4 -j 8
  ./bnb_beam_quality.py all -f "Run4_Files/*.root" -o run4test --label "Run 4"
  ./bnb_beam_quality.py fom beamqual_run4.root -o beamqual_run4_fom.root -j 8
  ./bnb_beam_quality.py plot run2_fom.root run4_fom.root --labels "Run 2" "Run 4" -o fom_plots
  ./bnb_beam_quality.py selftest Run2_Files/*.root Run4_Files/*.root

Run `./bnb_beam_quality.py <step> -h` for the options of each step.
"""
import argparse
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def _all(argv):
    p = argparse.ArgumentParser(prog="bnb_beam_quality.py all",
                                description="dump -> fom -> plot. Writes <out>_beamqual.root, "
                                            "<out>_fom.root and <out>_plots/.")
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("-d", "--defname", help="SAM definition")
    src.add_argument("-f", "--files", nargs="+", help="CAF files/URLs (globs ok)")
    src.add_argument("-l", "--filelist", help="text file with one CAF per line")
    p.add_argument("-o", "--out", required=True, help="output prefix")
    p.add_argument("--label", help="label used in the plots (default: the prefix)")
    p.add_argument("-e", "--experiment", help="SAM experiment (icarus/sbnd)")
    p.add_argument("-j", "--jobs", type=int, default=1)
    p.add_argument("-n", "--max-files", type=int, default=0)
    p.add_argument("--cut", type=float, default=0.95)
    p.add_argument("--no-fill-offsets", action="store_true")
    a = p.parse_args(argv)

    import dump_bnb_quality, bnb_fom, plot_fom_compare
    beam, fom, plots = f"{a.out}_beamqual.root", f"{a.out}_fom.root", f"{a.out}_plots"
    d = (["-d", a.defname] if a.defname else ["-f", *a.files] if a.files else ["-l", a.filelist])
    d += ["-o", beam, "-j", str(a.jobs)]
    if a.experiment:
        d += ["-e", a.experiment]
    if a.max_files:
        d += ["-n", str(a.max_files)]
    print(f"\n##### [1/3] dump -> {beam}", file=sys.stderr)
    dump_bnb_quality.main(d)
    print(f"\n##### [2/3] fom -> {fom}", file=sys.stderr)
    bnb_fom.main([beam, "-o", fom, "-j", str(a.jobs)] + (["--no-fill-offsets"] if a.no_fill_offsets else []))
    print(f"\n##### [3/3] plot -> {plots}/", file=sys.stderr)
    plot_fom_compare.main([fom, "--labels", a.label or os.path.basename(a.out), "-o", plots, "--cut", str(a.cut)])


def _selftest(argv):
    p = argparse.ArgumentParser(prog="bnb_beam_quality.py selftest",
                                description="dump + fom on the given CAF files. CAFs with the production "
                                            "beam-quality result (fom_version >= 1): pass if the Python Tier-1 "
                                            "FOM_best and fom_source equal the C++ ones for >= 99.9%% of the spills "
                                            "the C++ did not fill from neighbors. Older CAFs: pass if the recomputed "
                                            "legacy FOM reproduces the CAF FOM for >= 99.9%% of spills.")
    p.add_argument("files", nargs="+", help="a few flat CAF files")
    p.add_argument("--min-closure", type=float, default=0.999)
    a = p.parse_args(argv)
    import numpy as np
    import uproot
    import dump_bnb_quality, bnb_fom
    FILLED = 128 | 1024 | 16384     # Tier-1W fills: depend on the window, not compared

    def key(d):                     # spill time [ns]: unique within the input
        return d["spill_time_sec"].astype(np.int64) * 1_000_000_000 + d["spill_time_nsec"].astype(np.int64)

    with tempfile.TemporaryDirectory() as tmp:
        beam, fom = os.path.join(tmp, "beam.root"), os.path.join(tmp, "fom.root")
        dump_bnb_quality.main(["-f", *a.files, "-o", beam])
        bnb_fom.main([beam, "-o", fom])
        ok = True
        print("\n=== selftest ===")
        for tree in ("spill", "bnbinfo"):
            t, td = uproot.open(fom)[tree], uproot.open(beam)[tree]
            ids = ["spill_time_sec", "spill_time_nsec"]
            new_caf = "cpp_fom_version" in td.keys() and td.num_entries > 0 \
                and td["cpp_fom_version"].array(library="np").max() >= 1
            if new_caf:
                # the C++ result is in the dump only (bnb_fom.py keeps its input and ID columns)
                d = t.arrays(ids + ["FOM_best", "fom_source"], library="np")
                c = td.arrays(ids + ["cpp_FOM_best", "cpp_fom_source", "cpp_fom_status"], library="np")
                _, i, j = np.intersect1d(key(d), key(c), return_indices=True)
                sel = (c["cpp_fom_status"][j].astype(np.int64) & FILLED) == 0
                i, j = i[sel], j[sel]
                x, y = d["FOM_best"][i], c["cpp_FOM_best"][j]
                agree = ((np.isnan(x) & np.isnan(y)) | (np.abs(x - y) <= 1e-5)) \
                    & (d["fom_source"][i] == c["cpp_fom_source"][j])
                n = len(i)
                frac = agree.mean() if n else 1.0
                good = frac >= a.min_closure and n > 0
                print(f"{tree:8s} {n:8,d} spills   Python Tier 1 == C++ (FOM_best, fom_source) {frac:.4%}   "
                      f"{'PASS' if good else 'FAIL'}")
            else:
                # older production wrote either the FillExposure choice of the three estimators
                # (FOM_legacy) or the fitted FOM alone, with 4 for "no multiwire width"
                d = t.arrays(["FOM", "FOM_legacy", "FOM_legacy_fit", "FOM_best"], library="np")
                n = len(d["FOM"])
                choice = np.isclose(d["FOM"], d["FOM_legacy"], atol=1e-4)
                fit_only = np.isclose(d["FOM"], d["FOM_legacy_fit"], atol=1e-4) \
                    | ((d["FOM"] == 4) & (d["FOM_legacy_fit"] == -999))
                (close, conv) = max((choice.mean(), "estimator choice"), (fit_only.mean(), "fitted FOM only")) \
                    if n else (1.0, "-")
                has_new = np.isfinite(d["FOM_best"]).mean() if n else 0.0
                good = close >= a.min_closure and n > 0
                print(f"{tree:8s} {n:8,d} spills   legacy closure {close:.4%} ({conv})   "
                      f"new FOM defined {has_new:.1%}   {'PASS' if good else 'FAIL'}")
            ok &= good
        print("SELFTEST", "PASSED" if ok else "FAILED")
        return 0 if ok else 1


def main():
    steps = {"dump": "dump_bnb_quality", "fom": "bnb_fom", "plot": "plot_fom_compare",
             "fixups": "fom_fixups", "validate": "validate_fom", "compare-cpp": "compare_cpp_caf"}
    if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
        print(__doc__)
        return 0
    step, argv = sys.argv[1], sys.argv[2:]
    if step in steps:
        mod = __import__(steps[step])
        sys.argv[0] = f"bnb_beam_quality.py {step}"
        return mod.main(argv)
    if step == "all":
        return _all(argv)
    if step == "selftest":
        return _selftest(argv)
    print(f"unknown step '{step}'\n{__doc__}", file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main() or 0)
