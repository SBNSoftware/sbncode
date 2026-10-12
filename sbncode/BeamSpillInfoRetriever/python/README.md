# BNB beam quality: analysis-level tools (Python)

The BNB beam-quality procedure has two parts.

* **Tier 1 / Tier 1W** run in production, in the BNB retriever modules of this
  directory (`getFOM.cpp`, `BNBFOMFill.cpp`). Each spill gets `FOM_best`,
  `fom_source`, `fom_status`, `tor_used`, the beam at the target and `fom_version`
  in `sbn::BNBSpillInfo` and in the CAF (`rec.hdr.bnbinfo`, `rec.hdr.spillbnbinfo`).
* **Tier 2** runs once per data period, on the merged `bnbinfo` record of all its
  CAFs. It recovers spills that Tier 1 could not give a FOM, using the spills on both
  sides and, for one ICARUS outage, ACNet logger readings. Every method has passed a
  closure test.

This directory holds the Python side:
* the **reference implementation** of Tier 1 (`bnb_fom.py`), against which the C++ is
  validated spill by spill;
* **Tier 2** (`fom_fixups.py`, `logger_fill.py`);
* the tools to dump CAFs, validate the result and compare it with the production (C++)
  result.

The method, its validation and the numbers are documented in the BNB beam-quality
technical note.

They are plain Python (numpy, pandas, awkward, uproot, numba), read flat CAFs with
uproot, and do not need art or ROOT. They are installed with the product under
`$SBNCODE_DIR/python/bnb_beam_quality/`.

## Environment

The art/UPS release environment does not provide uproot, awkward or numba; use a
virtual environment:

```bash
python3 -m venv bq-venv && source bq-venv/bin/activate
pip install -r $SBNCODE_DIR/python/bnb_beam_quality/requirements.txt
export PATH=$SBNCODE_DIR/python/bnb_beam_quality:$PATH
```

On the GPVMs, `root://` inputs on `/pnfs` are read through the local mount
automatically. Elsewhere install `fsspec-xrootd`.

## Processing a data period

```bash
export BQ_DATA=<directory for the per-spill files>  BQ_WORK=<scratch directory>
# 1. CAFs (SAM definition, file list or files) -> per-spill dump, in chunks; merge
dump_bnb_quality.py -d <definition> -e sbn -j 32 --chunk-files 5000 -o $BQ_DATA/<s>.root
merge_bnb_dumps.py -l $BQ_DATA/<s>_chunks.txt -o $BQ_DATA/<s>_all.root
# 2. Tier 1 per spill (fast version; --check compares it with the reference first)
bnb_fom_fast.py --check <small dump>.root
bnb_fom_fast.py $BQ_DATA/<s>_all.root -o $BQ_DATA/<s>_fixed.root
# 3. Tier-1 corrections and Tier-2 recovery
fom_fixups.py $BQ_DATA/<s>_fixed.root -o $BQ_DATA/<s>_fixed_fx.root                  # SBND
fom_fixups.py $BQ_DATA/<s>_fixed.root -o $BQ_DATA/<s>_fixed_fx.root --groups 12 \
    --nb-min-fom 0.998 --burst-fill [--logger-fill]                                   # ICARUS
# 4. validation report
validate_fom.py $BQ_DATA/<s>_fixed_fx.root --label "<label>" -o $BQ_WORK/fom_validation
```

* **Merge before step 2.** Offsets are filled from the nearest valid spill of the same
  run *within the input file*, so the chunks must be merged first.
* **Logger fill.** `--logger-fill` (ICARUS Run 2, June 2023 outage only) needs
  `logger_fill.py calibrate $BQ_DATA/<s>_fixed.root`, run twice: the first pass
  caches the outage spills. The calibration is written to `$BQ_ACNET_DIR` (default
  `acnet/` next to the code, which is read-only once installed).

`bnb_beam_quality.py` is a single entry point for all the steps:
`dump`, `fom`, `plot`, `all`, `fixups`, `validate`, `compare-cpp`, `selftest`.

## Checking production

```bash
# a few new CAFs: Python Tier 1 must equal the C++ FOM_best / fom_source (exit 0 = pass)
bnb_beam_quality.py selftest <caf files>
# spill-by-spill comparison of a CAF with the production result against a reference CAF
# of the same input and this code's Tier 1 / Tier 2 outputs
compare_cpp_caf.py --new <new.flat.caf.root> --prod <reference.flat.caf.root> \
    --py1 <..._fixed.root> --py2 <..._fixed_fx.root> [--oldfmt <caf from an older BNBSpillInfo>]
```

On CAFs made before the production code carried the result (`fom_version` 0 or
absent), `selftest` instead checks that the recomputed legacy FOM reproduces the CAF
`FOM`.

## Files

| file | role |
|---|---|
| `bnb_beam_quality.py` | entry point for all the steps, and `selftest` |
| `dump_bnb_quality.py`, `merge_bnb_dumps.py` | CAFs → per-spill dump (trees `spill`, `bnbinfo`, with multiwire profiles); parallel, chunked, resumable; merging (use it instead of `hadd`) |
| `bnb_fom.py` | **reference** implementation: the legacy algorithm (reproduces the old CAF `FOM`) and the Tier-1 algorithm. Change it only by a documented, versioned update, together with `getFOM.cpp` |
| `bnb_fom_fast.py` | same results with numba kernels, 5–8× faster; `--check` compares with the reference |
| `fom_fixups.py` | target-BPM choice, TOR875 back-up, backup-BPM removal, empty-chamber and MMBTBB widths, and the Tier-2 recovery steps (neighbor width and fill, logger fill, 875-dropout fill) |
| `neighbor_fill.py` | neighbor fill on its own (also used by `fom_fixups.py`) |
| `logger_fill.py`, `acnet/` | logger fill for the June 2023 ICARUS outage: recovered one-minute ACNet readings, calibration, and the D44 query script (`acnet/README.md`) |
| `validate_fom.py` | validation report and plots |
| `compare_cpp_caf.py` | production (C++) result in a CAF vs a reference CAF and this code, spill by spill |
| `spill_event_map.py` | spill → (run, subrun, event) map |
| `plot_fom_compare.py` | old vs new FOM plots and summary |

The studies behind the technical note are kept with the note's code archive, not here:
- the threshold and closure studies;
- the device timelines;
- the POT tables;
- the flux study;
- the figure scripts.

## Output columns

The trees of `bnb_fom.py` / `fom_fixups.py` keep the input and ID columns they use and add:

| column | meaning |
|---|---|
| `FOM_best` | FOM of the spill in [0, 1]; NaN if none can be formed. No +100 offset |
| `fom_source` | width used: 1/2 fitted M876/M875 profile, 3/4 M876/M875 database width, 5 nominal; 0 fitted MMBTBB (legacy only); −1 none |
| `fom_status` | bit mask (below) |
| `tor_used` | intensity used [protons]: TOR860, or TOR875 after a TOR860 drop-out (`fom_fixups.py`) |
| `hpos_tgt`, `hang_tgt`, `vpos_tgt`, `vang_tgt` | beam position [mm] and angle [mrad] at the target center |
| `tgt_sigx`, `tgt_sigy` | beam width at the target used for the FOM [mm]; NaN: nominal |
| `FOM_legacy` (+ `_fit`, `_prefit`, `_nomw`, `fom_legacy_source`) | the old production algorithm |
| `FOM_fit`, `FOM_prefit`, `FOM_nomw` | Tier-1 FOM with each width estimator |

The dump (`dump_bnb_quality.py`) also keeps the production result of CAFs that carry it, as
`cpp_<name>` (`cpp_FOM_best`, `cpp_fom_status`, …, `cpp_fom_version`), next to the device readings.

`fom_status` bits (`sbn::fomstatus` in `getFOM.h`):

| bit | meaning |
|---|---|
| 1 | no valid intensity |
| 2 / 4 | horizontal / vertical position or angle cannot be formed |
| 8 | BPM offset filled |
| 16 | intensity below 1e11 |
| 32 | multiwire profile rejected, more than 33 ms from the spill (analysis level only) |
| 64 | nominal width |
| 128 | neighbor fill |
| 256 | TOR875 used |
| 512 | empty chamber |
| 1024 | width from the neighbors |
| 2048 | backup target BPM not used |
| 4096 | MMBTBB width replaced (analysis level only) |
| 8192 | ACNet logger fill (analysis level only) |
| 16384 | 875-station drop-out fill |
| 32768 | horizontal projection through HPTG2 |

A spill has a FOM if `FOM_best` is finite. Filled spills keep bits 2/4, which record why
they were filled.
