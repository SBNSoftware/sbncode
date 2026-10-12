#!/usr/bin/env python3
"""
fom_fixups.py -- post-process a bnb_fom(_fast).py output with the fixes found in the
edge-case scan (bnb_fom.py itself is left untouched).  Applied in this order:

 0. target H-BPM choice  the horizontal target projection uses the target BPM that autotune was calibrated
                          on most recently (as sbnana getBNBFoM): HPTG2 before 2023-03-03 23:45 UTC, HPTG1
                          after.  bnb_fom.py always uses HPTG1, whose stale pre-March offset puts the
                          projected beam ~2.6 mm off-center.  Spills in an HPTG2 era get hpos/hang from
                          HP875 -> HPTG2 (missing offsets: nearest valid in the run, else the autotune
                          table) and their FOM recomputed; bit 32768 HPTG2_PRIMARY.  --hptg1-only disables.
 1. TOR875 back-up        TOR860 < 0.5 * TOR875 (TOR860 dropout, reads ~1e9-1e10) -> the
                          intensity is TOR875; the FOM is recomputed; bit 256 TOR875_USED.
                          New branch tor_used = intensity used for the FOM: sum it for POT.
 2. second-BPM fallback   spills whose FOM used HPTG2 / VPTG2 because HPTG1 / VP873 had no
                          reading (biased position) lose that FOM.  With step 0 the 'first' horizontal
                          BPM is the era's primary one (HPTG2 before March 2023, HPTG1 after) (bit 2048 FALLBACK_REMOVED,
                          NO_HBPM / NO_VBPM set) and become candidates for step 5.
 3. empty multiwire       database widths (M876/M875 HS, VS) > 4 mm mean the multiwire saw no
                          beam (beam: 1-2.5 mm, empty profile: ~7 mm).  Spills whose measured-
                          width FOM used such a width (sources 3, 4) get the neighbors' width
                          (step 4 rules) or else the nominal width (source 5, bit 64);
                          bit 512 MW_EMPTY.
 4. MMBTBB width          spills whose width came from an MMBTBB profile fit (source 0, not
                          reproducible; in SBND run 18576 it gives ~4.7 mm and FOM ~0.6) take
                          the neighbors' width when available, else the nominal width
                          (source 5, bits 64 + 4096 MMBTBB_DROPPED).
    Neighbor width (3, 4): both adjacent spills (<= 0.1 s, same run) have a valid FOM with a
    measured width from M876/M875 (sources 1-4, not empty), and their widths agree within
    0.1 mm in both planes; the spill gets their mean width and source; bit 1024 WIDTH_FROM_NB.
 5. neighbor fill        as neighbor_fill.py (bit 128): BPM-missing spills between two
                          adjacent (<= 0.1 s) stable spills get their mean position/angle.
 6. logger fill (opt.)    --logger-fill: spills in the June 2023 ICARUS outage (HP875/VP875/VP873
                          missing in IFBeam) get position/angle from their own target BPMs plus the
                          beam angle from once-a-minute ACNet readings (logger_fill.py; bit 8192).
                                                    Only spills between two ACNet samples (<= 60 s each side) whose angles agree
                          (0.2 mrad), with target BPMs within 0.5 mm of the samples', and estimated FOM
                          > 0.999 are filled.
 7. 875 dropouts (opt.)   --burst-fill: spills with HP875 and VP875 missing but HPTG1/VPTG2 read take the
                          position offset and angle of the measured spills on both sides (<= 60 s, both
                          FOM > 0.999, angles within 0.2 mrad, target BPMs within 0.5 mm); filled only if
                          the estimate is > 0.995 (bit 16384).
    --nb-min-fom X: step 5 only fills when both neighbors have FOM > X (ICARUS: 0.998).
FOMs are recomputed with the spill's own position/angle (steps 1, 3, 4) and the fixed
algorithm.  Both trees are written; the spill tree is matched by run/subrun/spill time.

  ./fom_fixups.py sbnd_run1_part1_fixed.root -o sbnd_run1_part1_fixed_fx.root
  ./fom_fixups.py icarus_run2_full_fixed.root -o icarus_run2_full_fixed_fx.root --groups 12 \
        --nb-min-fom 0.998 --logger-fill --burst-fill
  ./fom_fixups.py X_fixed.root --stats --groups 12 --group 0   (counts only, one run group)
"""
import argparse
import os
import sys

import numpy as np
import uproot

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from neighbor_fill import keyhash, CUTS as NB_CUTS          # noqa: E402

BITS = {128: "NEIGHBOR_FILLED", 256: "TOR875_USED", 512: "MW_EMPTY", 1024: "WIDTH_FROM_NB",
        2048: "FALLBACK_REMOVED", 4096: "MMBTBB_DROPPED", 8192: "LOGGER_FILL", 16384: "BURST_FILL",
        32768: "HPTG2_PRIMARY"}
NB_MIN_FOM = None      # --nb-min-fom: neighbor fill only if both neighbors have FOM above this
LOGGER = None          # --logger-fill: logger_fill.load() result
BURST = False          # --burst-fill: 875-station dropout fill (logger_fill.burst_fill_group)
HPTG_RULE = True       # step 0 (disable with --hptg1-only)
# autotune calibration change times (s, UTC) and values [mm] of the target-area BPMs, from sbnana
# SBNAna/Vars/getBNBFoM.cxx (Z. Pavlovic, IFBeam).  Used to pick the primary horizontal target BPM and as a
# last-resort offset when the CAF offset is -999 for a whole run.  The CAF offsets agree with these to
# < 0.005 mm in ICARUS Run 2 (fom_validation/sbnana_compare/offsets_vs_sbnana.md).
AUTOTUNE = {
    "HP875": ([1420092000, 1574190445, 1576014670, 1588603117, 1605634657, 1606510357, 1608593857, 1668802510,
               1711030691], [-3.40, -1.20, -1.20, -3.33, -1.20, 0.75, -2.40, -5.53, -3.42]),
    "HPTG1": ([1420092000, 1576014670, 1606510357, 1677887110, 1711030691], [0.46, 0.0, -2.15, -1.16, -0.03]),
    "HPTG2": ([1420092000, 1574190444, 1576014670, 1588603117, 1605706657, 1606510357, 1608593857, 1668802510],
              [1.07, -0.40, -0.40, 0.0, -0.40, -2.31, 0.0, -0.84]),
}
COLS = ["run", "subrun", "spill_time_sec", "spill_time_nsec", "TOR860", "TOR875", "FOM_best", "fom_status",
        "fom_source", "hpos_tgt", "hang_tgt", "vpos_tgt", "vang_tgt", "tgt_sigx", "tgt_sigy",
        "M876HS", "M876VS", "M875HS", "M875VS", "HP875", "HPTG1", "HPTG2", "VP875", "VP873", "VPTG2", "THCURR",
        "HP875Offset", "HPTG2Offset"]
UPD = ["FOM_best", "fom_status", "fom_source", "hpos_tgt", "hang_tgt", "vpos_tgt", "vang_tgt", "tgt_sigx", "tgt_sigy"]
MIN_TOR, BEAM_TOR, EMPTY_MW, DROP = 0.1e12, 1e11, 4.0, 0.5
# beam-quality selection used in the report: TOR860 > 1e11 and TOR875 > 1e11, 173 <= THCURR <= 175 kA,
# 0.98 < FOM <= 1 (LM875A/B/C > 1e-2 is not applied: the loss monitors are not in the dumps)
SEL_TOR, SEL_THCURR, SEL_FOM = 1e11, (173.0, 175.0), 0.98


def ok(v):
    v = np.asarray(v, float)
    return np.isfinite(v) & (v != -999)


def tor_used(t860, t875):
    t860, t875 = np.asarray(t860, float), np.asarray(t875, float)
    g860, g875 = ok(t860) & (t860 > 0), ok(t875) & (t875 > 0)
    drop = g860 & g875 & (t860 < DROP * t875)
    tor = np.where(g860 & ~drop, t860, np.where(g875, t875, 0.0))
    return tor, drop


def autotune(dev, tsec):
    """(value, calibration time) of the autotune offset in force at tsec."""
    T, V = AUTOTUNE[dev]
    i = np.clip(np.searchsorted(np.asarray(T), np.asarray(tsec, float), side="right") - 1, 0, len(T) - 1)
    return np.asarray(V)[i], np.asarray(T)[i]


def hptg2_primary(tsec):
    """sbnana rule: HPTG2 is the primary horizontal target BPM when its latest calibration is at least as
    recent as HPTG1's (true before 2023-03-03 23:45 UTC, false from then on)."""
    return autotune("HPTG2", tsec)[1] >= autotune("HPTG1", tsec)[1]


def fill_nearest(v, run, t, fallback):
    """-999/NaN entries -> nearest-in-time valid value in the same run (rows sorted by run, time);
    runs with no valid value take `fallback`."""
    v = np.asarray(v, float)
    n = len(v)
    good = ok(v)
    if good.all():
        return v
    idx = np.arange(n)
    last = np.maximum.accumulate(np.where(good, idx, -1))
    nxt = np.minimum.accumulate(np.where(good, idx, n)[::-1])[::-1]
    lc, nc = np.clip(last, 0, n - 1), np.clip(nxt, 0, n - 1)
    dl = np.where((last >= 0) & (run[lc] == run), t - t[lc], np.inf)
    dn = np.where((nxt < n) & (run[nc] == run), t[nc] - t, np.inf)
    val = np.where(dl <= dn, v[lc], v[nc])
    return np.where(good, v, np.where(np.isfinite(np.minimum(dl, dn)), val, fallback))


def load_group(path, groups, g):
    t = uproot.open(path)["bnbinfo"]
    keys = set(t.keys())
    cols = [c for c in COLS if c in keys]
    parts = []
    for d in t.iterate(cols, step_size=4_000_000, library="np"):
        s = (d["run"] % groups) == g
        parts.append({k: v[s] for k, v in d.items()})
        del d
    d = {k: np.concatenate([p[k] for p in parts]) for k in cols}
    T = d["spill_time_sec"].astype(np.int64) * 1_000_000_000 + d["spill_time_nsec"].astype(np.int64)
    o = np.lexsort((T, d["run"]))
    d = {k: v[o] for k, v in d.items()}
    T = T[o]
    keep = np.ones(len(T), bool)
    keep[1:] = ~((T[1:] == T[:-1]) & (d["run"][1:] == d["run"][:-1]))
    d = {k: v[keep] for k, v in d.items()}
    d["t"] = T[keep] / 1e9
    return d


def fix_group(d, S):
    import bnb_fom_fast as F
    n = len(d["t"])
    f = lambda k: d[k].astype(float).copy()
    fom, st, src = f("FOM_best"), d["fom_status"].astype(np.int32).copy(), d["fom_source"].astype(np.int8).copy()
    pos = {k: f(k) for k in ("hpos_tgt", "hang_tgt", "vpos_tgt", "vang_tgt")}
    sx, sy = f("tgt_sigx"), f("tgt_sigy")
    fom0, st0, src0 = fom.copy(), st.copy(), src.copy()
    tor, drop = tor_used(d["TOR860"], d["TOR875"])
    beam = tor > BEAM_TOR
    t860, t875 = d["TOR860"].astype(float), d["TOR875"].astype(float)   # intensity as bnb_fom.py (fixed mode) used it
    tor_old = np.where(ok(t860) & (t860 > 0), t860, np.where(ok(t875) & (t875 > 0), t875, 0.0))
    S["pot"] += tor[beam].sum()
    run = d["run"]
    i = np.arange(n)
    p, q = np.clip(i - 1, 0, n - 1), np.clip(i + 1, 0, n - 1)
    adj = (i > 0) & (i < n - 1) & (run[p] == run) & (run[q] == run) & \
          (d["t"] - d["t"][p] <= NB_CUTS["max_gap"]) & (d["t"][q] - d["t"] <= NB_CUTS["max_gap"])
    recompute = np.zeros(n, bool)
    prev = (st0 & 128) != 0                      # already neighbor-filled input (_nf file)

    # 0. primary horizontal target BPM (HPTG2 before 2023-03-03, HPTG1 after; see AUTOTUNE)
    use2 = hptg2_primary(d["spill_time_sec"]) if HPTG_RULE else np.zeros(n, bool)
    h875, h1, h2 = ok(d["HP875"]), ok(d["HPTG1"]), ok(d["HPTG2"])
    sw = use2 & np.isfinite(fom) & ~prev & h875 & h2       # re-project via HPTG2 (incl. bnb_fom's own fallback)
    if sw.any():
        Z = F.ref.Z
        tsec = d["spill_time_sec"].astype(float)
        o875 = fill_nearest(d.get("HP875Offset", np.full(n, -999.)), run, d["t"], autotune("HP875", tsec)[0])
        o2 = fill_nearest(d.get("HPTG2Offset", np.full(n, -999.)), run, d["t"], autotune("HPTG2", tsec)[0])
        x0 = d["HP875"].astype(float) - o875
        x2 = d["HPTG2"].astype(float) - o2
        ang = (x2 - x0) / (Z["hptg2"] - Z["hp875"])                       # mrad (mm / m)
        hp = x0 + ang * (Z["tgt"] - Z["hp875"])
        S["hptg2_dpos_sum"] += float(np.abs(hp - pos["hpos_tgt"])[sw & beam].sum())
        pos["hang_tgt"][sw], pos["hpos_tgt"][sw] = ang[sw], hp[sw]
        recompute |= sw
        st[sw] |= 32768
    S["hptg2_pot"] += float(tor[sw & beam].sum())
    S["hptg2_spills"] += int((sw & beam).sum())

    # 1. TOR875 back-up
    st[drop] |= 256
    st[drop] &= ~16
    st[drop & (tor < MIN_TOR)] |= 16
    recompute |= drop & np.isfinite(fom)
    S["tor875_spills"] += int((drop & beam).sum())
    S["tor875_pot_gain"] += float((tor - tor_old)[drop & beam].sum())

    # 2. second-BPM fallback -> no FOM (horizontal: the era's primary target BPM missing, the other read)
    fbh = np.isfinite(fom) & ~prev & h875 & np.where(use2, ~h2 & h1, ~h1 & h2)
    fbv = np.isfinite(fom) & ~prev & ok(d["VP875"]) & ~ok(d["VP873"]) & ok(d["VPTG2"])
    fb = fbh | fbv
    fom[fb] = np.nan
    src[fb] = -1
    st[fb] |= 2048
    st[fbh] |= 2
    st[fbv] |= 4
    S["fallback_pot"] += float(tor[fb & beam].sum())

    # 3./4. widths
    valid = np.isfinite(fom)
    db = lambda a, b: (ok(d[a]) & (d[a].astype(float) > EMPTY_MW)) | (ok(d[b]) & (d[b].astype(float) > EMPTY_MW))
    wide = valid & (((src == 3) & db("M876HS", "M876VS")) | ((src == 4) & db("M875HS", "M875VS")))
    mmb = valid & (src == 0)
    goodw = valid & (src >= 1) & (src <= 4) & ~wide
    nbw = adj & goodw[p] & goodw[q] & (np.abs(sx[p] - sx[q]) <= 0.1) & (np.abs(sy[p] - sy[q]) <= 0.1)
    for m, tag in ((wide, "wide"), (mmb, "mmb")):
        use = m & nbw
        sx[use], sy[use] = 0.5 * (sx[p] + sx[q])[use], 0.5 * (sy[p] + sy[q])[use]
        src[use] = src[p][use]
        st[use] |= 1024
        recompute |= use
        S[tag + "_pot"] += float(tor[m & beam].sum())
        S[tag + "_nb_pot"] += float(tor[use & beam].sum())
    mmb_nom = mmb & ~nbw          # whole-run MMBTBB use (e.g. SBND run 18576): no neighbor width
    st[mmb_nom] |= 4096
    nom = (wide | mmb_nom) & ~nbw
    sx[nom], sy[nom] = np.nan, np.nan
    src[nom] = 5
    st[nom] |= 64
    st[wide] |= 512
    recompute |= nom

    # recompute FOMs with the spill's own position
    r = np.where(recompute & np.isfinite(fom))[0]
    if len(r):
        f2 = F.ref.calc_fom(pos["hpos_tgt"][r], pos["hang_tgt"][r], pos["vpos_tgt"][r], pos["vang_tgt"][r],
                            tor[r], sx[r], sy[r], pi=np.pi)
        fom[r] = 1 - 10 ** f2

    # 5. neighbor fill
    valid = np.isfinite(fom)
    cand = beam & ~valid & ((st & 6) != 0) & ((st & 1) == 0)
    S["nan_bpm_pot"] += float(tor[cand].sum())
    c = adj & cand & valid[p] & valid[q]
    tm = 0.5 * (tor[p] + tor[q])
    with np.errstate(invalid="ignore", divide="ignore"):
        c &= np.abs(tor - tm) / tm <= NB_CUTS["dtor"]
    dd = lambda a: np.abs(a[p] - a[q])
    c &= np.maximum(dd(pos["hpos_tgt"]), dd(pos["vpos_tgt"])) <= NB_CUTS["dpos"]
    c &= np.maximum(dd(pos["hang_tgt"]), dd(pos["vang_tgt"])) <= NB_CUTS["dang"]
    c &= np.nan_to_num(np.maximum(dd(sx), dd(sy)), nan=0.0) <= NB_CUTS["dsig"]
    c &= dd(fom) <= NB_CUTS["dfom"]
    if NB_MIN_FOM is not None:
        c &= (np.nan_to_num(fom[p], nan=-1) > NB_MIN_FOM) & (np.nan_to_num(fom[q], nan=-1) > NB_MIN_FOM)
    r = np.where(c)[0]
    if len(r):
        for k in pos:
            pos[k][r] = 0.5 * (pos[k][p[r]] + pos[k][q[r]])
        own = np.isfinite(sx[r]) & np.isfinite(sy[r]) & ((st[r] & 2048) == 0)
        sx[r] = np.where(own, sx[r], 0.5 * (sx[p[r]] + sx[q[r]]))
        sy[r] = np.where(own, sy[r], 0.5 * (sy[p[r]] + sy[q[r]]))
        f2 = F.ref.calc_fom(pos["hpos_tgt"][r], pos["hang_tgt"][r], pos["vpos_tgt"][r], pos["vang_tgt"][r],
                            tor[r], sx[r], sy[r], pi=np.pi)
        fom[r] = 1 - 10 ** f2
        src[r] = src[p[r]]
        st[r] |= 128
    S["nb_fill_pot"] += float(tor[c].sum())
    S["nb_fill_fallback_pot"] += float(tor[c & ((st & 2048) != 0)].sum())

    # 6. logger fill (June 2023 ICARUS outage only; see logger_fill.py)
    if LOGGER is not None:
        import logger_fill
        logger_fill.fill_group(d, fom, st, src, pos, sx, sy, tor, beam, S, LOGGER)

    # 7. 875-station dropout fill (ICARUS; see logger_fill.burst_fill_group)
    if BURST:
        import logger_fill
        measured = (st & (128 | 8192)) == 0
        logger_fill.burst_fill_group(d, fom, st, src, pos, sx, sy, tor, beam, S, measured)

    # bookkeeping
    S["nan_pot_before"] += float(tor_old[(tor_old > BEAM_TOR) & ~np.isfinite(fom0)].sum())
    S["pot_before"] += float(tor_old[tor_old > BEAM_TOR].sum())
    S["nan_pot_after"] += float(tor[beam & ~np.isfinite(fom)].sum())
    S["p95_before"] += float(tor_old[(tor_old > BEAM_TOR) & (np.nan_to_num(fom0, nan=-1) > 0.95)].sum())
    S["p95_after"] += float(tor[beam & (np.nan_to_num(fom, nan=-1) > 0.95)].sum())
    # official beam-quality selection (TOR860 as stored, and with the TOR875 back-up)
    t860z, t875z = np.where(ok(t860), t860, 0.0), np.where(ok(t875), t875, 0.0)
    hc = d["THCURR"].astype(float) if "THCURR" in d else np.full(n, 174.0)
    hc_ok = (hc >= SEL_THCURR[0]) & (hc <= SEL_THCURR[1])
    sel = (t860z > SEL_TOR) & (t875z > SEL_TOR) & hc_ok
    sel_tu = (tor > SEL_TOR) & (t875z > SEL_TOR) & hc_ok
    win = lambda v: (np.nan_to_num(v, nan=-1) > SEL_FOM) & (np.nan_to_num(v, nan=-1) <= 1)
    S["sel_pot"] += float(t860z[sel].sum())
    S["sel_pass_before"] += float(t860z[sel & win(fom0)].sum())
    S["sel_pass_after"] += float(t860z[sel & win(fom)].sum())
    S["sel_nan_before"] += float(t860z[sel & ~np.isfinite(fom0)].sum())
    S["sel_nan_after"] += float(t860z[sel & ~np.isfinite(fom)].sum())
    S["sel_tu_pot"] += float(tor[sel_tu].sum())
    S["sel_tu_pass"] += float(tor[sel_tu & win(fom)].sum())
    ch = ~((fom == fom0) | (np.isnan(fom) & np.isnan(fom0))) | (st != st0) | (src != src0)
    k = np.where(ch)[0]
    out = {"h": keyhash(run[k], d["subrun"][k], d["spill_time_sec"][k], d["spill_time_nsec"][k]),
           "FOM_best": fom[k], "fom_status": st[k], "fom_source": src[k], "tgt_sigx": sx[k], "tgt_sigy": sy[k]}
    for kk in pos:
        out[kk] = pos[kk][k]
    return out


def find_updates(path, groups=1, only_group=None):
    S = {k: 0.0 for k in ("pot", "pot_before", "hptg2_pot", "hptg2_spills", "hptg2_dpos_sum", "tor875_spills", "tor875_pot_gain", "fallback_pot", "wide_pot", "wide_nb_pot",
                          "mmb_pot", "mmb_nb_pot", "nan_bpm_pot", "nb_fill_pot", "nb_fill_fallback_pot",
                          "nan_pot_before", "nan_pot_after", "p95_before", "p95_after",
                          "sel_pot", "sel_pass_before", "sel_pass_after", "sel_nan_before", "sel_nan_after",
                          "sel_tu_pot", "sel_tu_pass")}
    outs = []
    for g in range(groups):
        if only_group is not None and g != only_group:
            continue
        d = load_group(path, groups, g)
        if len(d["t"]) == 0:
            continue
        outs.append(fix_group(d, S))
        del d
        print(f"[fixups] group {g + 1}/{groups} done", file=sys.stderr, flush=True)
    U = {k: np.concatenate([o[k] for o in outs]) for k in outs[0]} if outs else None
    if U is not None:
        o = np.argsort(U["h"])
        U = {k: v[o] for k, v in U.items()}
    return U, S


def report(S, label=""):
    P, P0 = S["pot"], S["pot_before"]
    pc = lambda x: f"{100 * x / P:.3f}%"
    L = [f"## fixups {label}", "",
         f"primary target H-BPM (HPTG2 before 2023-03-03): {pc(S.get('hptg2_pot', 0))} of POT "
         f"({int(S.get('hptg2_spills', 0)):,} spills) re-projected via HPTG2, mean |change in h-pos| "
         f"{S.get('hptg2_dpos_sum', 0) / max(S.get('hptg2_spills', 0), 1):.2f} mm" + ("" if HPTG_RULE else " [disabled]"),
         f"beam POT {P0:.4e} -> {P:.4e} (+{(P - P0) / P0:.3%} from TOR875 back-up, {int(S['tor875_spills']):,} spills)",
         f"second-BPM fallback FOMs removed: {pc(S['fallback_pot'])} of POT; refilled from neighbours: {pc(S['nb_fill_fallback_pot'])}",
         f"empty-multiwire widths replaced: {pc(S['wide_pot'])} of POT (neighbour width {pc(S['wide_nb_pot'])}, rest nominal)",
         f"MMBTBB widths: {pc(S['mmb_pot'])} of POT, replaced by neighbour width {pc(S['mmb_nb_pot'])}, rest nominal",
         f"neighbour fill (BPM missing): {pc(S['nb_fill_pot'])} filled of {pc(S['nan_bpm_pot'])} missing"
         + (f" (neighbours required to have FOM > {NB_MIN_FOM})" if NB_MIN_FOM is not None else ""),
         f"logger fill (June 2023 outage): {pc(S.get('logger_pot', 0))} filled of {pc(S.get('logger_cand_pot', 0))} candidates, "
         f"{pc(S.get('logger_rejected_pot', 0))} left unfilled (estimate not > 0.999)",
         f"875-station dropout fill: {pc(S.get('burst_pot', 0))} filled of {pc(S.get('burst_cand_pot', 0))} candidates",
         f"POT without FOM: {S['nan_pot_before'] / P0:.3%} -> {S['nan_pot_after'] / P:.3%}",
         f"POT with FOM_best > 0.95: {S['p95_before'] / P0:.3%} -> {S['p95_after'] / P:.3%}"]
    if S.get("sel_pot", 0) > 0:
        Q = S["sel_pot"]
        L += ["", f"beam-quality selection (TOR860, TOR875 > 1e11; {SEL_THCURR[0]:g} <= THCURR <= {SEL_THCURR[1]:g} kA; "
              f"{SEL_FOM} < FOM <= 1; LM875 not applied): {Q:.4e} POT",
              f"  POT passing the FOM window: {S['sel_pass_before'] / Q:.3%} -> {S['sel_pass_after'] / Q:.3%}",
              f"  POT without FOM: {S['sel_nan_before'] / Q:.3%} -> {S['sel_nan_after'] / Q:.3%}",
              f"  intensity cut on tor_used instead of TOR860: {S['sel_tu_pot']:.4e} POT "
              f"(+{S['sel_tu_pot'] / Q - 1:.3%}), {S['sel_tu_pass'] / S['sel_tu_pot']:.3%} pass the FOM window"]
    print("\n".join(L))
    return "\n".join(L)


def apply(path, out, U, step=500_000):
    tmp = out + ".partial.root"
    with uproot.open(path) as fin, uproot.recreate(tmp) as fo:
        for tr in ("spill", "bnbinfo"):
            if tr not in fin:
                continue
            t = fin[tr]
            types = {k: t[k].interpretation.numpy_dtype for k in t.keys()}
            types["tor_used"] = np.dtype("f8")
            wt = fo.mktree(tr, types)
            nup = 0
            for d in t.iterate([k for k in types if k != "tor_used"], step_size=step, library="np"):
                d["tor_used"] = tor_used(d["TOR860"], d["TOR875"])[0]
                if U is not None and len(U["h"]):
                    h = keyhash(d["run"], d["subrun"], d["spill_time_sec"], d["spill_time_nsec"])
                    j = np.clip(np.searchsorted(U["h"], h), 0, len(U["h"]) - 1)
                    hit = U["h"][j] == h
                    if hit.any():
                        jj = j[hit]
                        for k in UPD:
                            if k in d:
                                d[k] = d[k].copy()
                                d[k][hit] = U[k][jj].astype(d[k].dtype)
                        nup += int(hit.sum())
                wt.extend({k: np.asarray(d[k]).astype(types[k]) for k in types})
            print(f"[fixups] {tr}: {nup:,} entries updated", file=sys.stderr)
        for k in fin.keys(cycle=False):
            if k in ("spill", "bnbinfo"):
                continue
            s = str(fin[k])
            if k == "fom_status_bits":
                for b, nm in BITS.items():
                    if nm not in s:
                        s += f"\n{b}:{nm}"
            fo[k] = s
    os.replace(tmp, out)


def run(path, out, groups=1, label=""):
    with uproot.open(path) as fh:
        if "fom_status_bits" in fh and "NEIGHBOR_FILLED" in str(fh["fom_status_bits"]):
            print("[fixups] WARNING: input was already post-processed (neighbour-filled); use the "
                  "plain bnb_fom output (*_fixed.root) as input for correct results", file=sys.stderr)
    U, S = find_updates(path, groups)
    txt = report(S, label)
    apply(path, out, U)
    with open(os.path.splitext(out)[0] + "_fixups.md", "w") as fh:
        fh.write(txt + "\n")
    return S


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input")
    ap.add_argument("-o", "--output")
    ap.add_argument("--groups", type=int, default=1, help="process runs in this many groups (memory)")
    ap.add_argument("--group", type=int, help="with --stats: only this group")
    ap.add_argument("--stats", action="store_true", help="only count what would change (no output file)")
    ap.add_argument("--label", default="")
    ap.add_argument("--nb-min-fom", type=float, help="neighbour fill only if both neighbours' FOM exceed this "
                    "(ICARUS: 0.998, closure false accepts 0.003%%)")
    ap.add_argument("--burst-fill", action="store_true", help="fill short HP875+VP875 dropouts from the "
                    "measured spills on both sides (ICARUS; validated with burst_closure)")
    ap.add_argument("--hptg1-only", action="store_true", help="old behaviour: horizontal projection always "
                    "HP875->HPTG1 (skip step 0)")
    ap.add_argument("--logger-fill", action="store_true", help="fill the June 2023 ICARUS outage from ACNet "
                    "logger readings (analysis/acnet/; run logger_fill.py calibrate first)")
    a = ap.parse_args(argv)
    global NB_MIN_FOM, LOGGER, BURST, HPTG_RULE
    NB_MIN_FOM = a.nb_min_fom
    HPTG_RULE = not a.hptg1_only
    BURST = a.burst_fill
    if a.logger_fill:
        import logger_fill
        LOGGER = logger_fill.load()
    if a.stats:
        import json
        U, S = find_updates(a.input, a.groups, a.group)
        print(json.dumps(S))
        return
    if not a.output or os.path.abspath(a.output) == os.path.abspath(a.input):
        ap.error("-o OUTPUT (different from the input) required")
    run(a.input, a.output, a.groups, a.label)
    print(f"wrote {a.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
