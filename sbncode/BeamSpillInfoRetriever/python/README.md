# BNB beam quality: updating flat CAFs (Python)

`update_flatcafs.py` takes the flat CAFs of a SAM definition (or a file list) and writes a
copy of each, with the BNB beam-quality result of the full procedure in
`rec.hdr.bnbinfo` (every spill) and `rec.hdr.spillbnbinfo` (the triggered spill):

`FOM_best`, `fom_source`, `fom_status`, `tor_used`, `hpos_tgt`, `hang_tgt`, `vpos_tgt`,
`vang_tgt`, `tgt_sigx`, `tgt_sigy`, `fom_version`

These are the fields of `caf::SRBNBInfo` that the production retrievers fill (Tier 1 and the
end-of-subrun Tier 1W, `getFOM.h`, `BNBFOMFill.h`). This tool adds the Tier-2 recovery, which
needs every spill of a data period, and fills them in CAFs made before the retrievers had
them. Missing values are NaN and `fom_source` is −1. Nothing else in the file changes. A
`TNamed` called `bnb_beam_quality` records the tool, options, definition and date.

Unless `--no-bunch-rotation` is given, the Booster bunch-rotation settings `BRRMPS8` (scale factor, 0 = off), `BRRMPQ1` (frequency, kHz) and `BRRMPL9` (start time, µs) are also set for each spill (`bunch_rotation.py`). They come from IFBeam from 2025-10-09 on, and from the ACNet data logger before that (where `BRRMPL9` is not logged). This needs network access to `dbweb1.fnal.gov` and `www-bd.fnal.gov`. The production retrievers read them only when a bundle holding them is configured (`BunchRotationBundle`).

Steps (all files are processed before any is written, because Tier 2 needs the whole period):
1. dump the beam information of every CAF and merge it (`dump_bnb_quality.py`, `merge_bnb_dumps.py`);
2. compute Tier 1 per spill (`bnb_fom.py`, the reference implementation, sped up by `bnb_fom_fast.py`);
3. run Tier 2 over the whole set (`fom_fixups.py`, `neighbor_fill.py`, `logger_fill.py` + `acnet/`);
4. copy each CAF with the fields set from the resulting per-spill table.

## Environment

Steps 1–3 need numpy, pandas, awkward, uproot and numba (`requirements.txt`, in a virtual
environment). The copy needs ROOT with PyROOT. On the SBND GPVMs: the spack ROOT plus a venv,
with `samweb` on the `PATH` for `-d`.

## Usage

```bash
update_flatcafs.py -d <flat CAF definition> -e sbnd   --detector sbnd   -o <outdir> -j 16
update_flatcafs.py -d <flat CAF definition> -e icarus --detector icarus -o <outdir> -j 16 [--logger-fill]
update_flatcafs.py -f "dir/*.flat.caf.root" --detector sbnd -o <outdir>        # files instead of SAM
```

Detector settings:
* `--detector icarus` fills a spill from its neighbors only if both neighbors have FOM > 0.998, and turns on the 875-station dropout fill.
* `--detector sbnd` applies neither.
* `--logger-fill` is for ICARUS Run 2 only (the June 2023 outage); the calibration in `acnet/` is used.

Other options:
* `--chunk-files N` dumps large definitions N files at a time.
* `--groups N` runs Tier 2 in N groups of runs, to limit memory.
* `-n N` processes only the first N files.

Intermediate files go to `<outdir>/bnb_work/`, including the per-spill table `tier2.root` and its report `tier2_fixups.md`. A rerun reuses the finished steps and skips CAFs already written.

## Fields

| field | meaning |
|---|---|
| `FOM_best` | fraction of the beam on target, in [0, 1]; NaN if none can be formed. Cut on this |
| `fom_source` | width used: 1/2 fitted M876/M875 profile, 3/4 M876/M875 database width, 5 nominal; −1 none |
| `fom_status` | bit mask, `sbn::fomstatus` in `getFOM.h` |
| `tor_used` | intensity used for the FOM [protons]: TOR860, or TOR875 after a TOR860 drop-out |
| `hpos_tgt`, `hang_tgt`, `vpos_tgt`, `vang_tgt` | beam position [mm] and angle [mrad] at the target center (NaN for a plane that could not be formed) |
| `tgt_sigx`, `tgt_sigy` | beam width at the target [mm]; NaN: nominal width |
| `fom_version` | 1 (algorithm version, as `sbn::FOMVersion`); 0 for a spill not found in the table |
| `BRRMPS8`, `BRRMPQ1`, `BRRMPL9` | Booster bunch-rotation scale factor (0 = off), frequency [kHz] and start time [µs] in force at the spill; NaN if not available |

`fom_status` bits:

| bit | meaning |
|---|---|
| 1 | no valid intensity |
| 2 / 4 | horizontal / vertical position or angle missing |
| 8 | offset filled |
| 16 | intensity < 1e11 |
| 32 | multiwire profile rejected (> 33 ms from the spill) |
| 64 | nominal width |
| 128 | neighbor fill |
| 256 | intensity from TOR875 (TOR860 dropped out, ≤ 0 or missing) |
| 512 | empty chamber |
| 1024 | width from the neighbors |
| 2048 | backup target BPM not used |
| 4096 | MMBTBB width replaced |
| 8192 | ACNet logger fill |
| 16384 | 875-station dropout fill |
| 32768 | horizontal projection through HPTG2 |

Filled spills keep bits 2/4, which record why they were filled.

`bnb_fom.py` is the reference implementation of Tier 1 against which the C++ (`getFOM.cpp`) is
validated. Change the two together, by a documented, versioned update.
