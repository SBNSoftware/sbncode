#!/usr/bin/env python3
"""
validate_fom.py -- streaming validation + old-vs-new comparison of the BNB FOM for
full-size bnb_fom(_fast).py outputs (tens of millions of spills, small memory).

  ./validate_fom.py icarus_run2_full_fixed.root --label "ICARUS Run 2" -o fom_validation
  ./validate_fom.py sbnd_run1_part1_fixed.root  --label "SBND Run 1 (part 1)" -o fom_validation

The bnbinfo tree is read in batches and only histograms / sums are kept.  Progress is
checkpointed after every batch (<outdir>/<tag>/state.pkl), so an interrupted run continues
where it stopped when the same command is repeated (--restart to start over).

Checks (report.md):
  closure      FOM_legacy (production algorithm re-run) == FOM stored in the CAF
  sanity       FOM_best in [0,1]; NaN exactly when fom_source == none; no FOM without beam
  duplicates   the same (run, subrun, spill time) stored more than once (inflates POT)
  trigger      FOM of each triggered spill in the 'spill' tree == its 'bnbinfo' entry
Comparisons (PNG + one PDF): category accounting and migration, FOM distributions,
per-spill change, POT passing FOM > cut per run and per day, why spills changed
(status bits), beam position/width of recovered spills vs spills that were always good,
nominal-width vs measured-width FOM.

Old = CAF FOM; categories as FillExposure: (0,1] measured width, >100 nominal width
(100 + FOM), anything else (-999; 0/4 in older SBND productions) no FOM.
New = FOM_best (no +100); fom_source 0-4 measured, 5 nominal, -1 none.
POT = TOR860 (TOR875 if TOR860 <= 0 or missing), spills with positive intensity only.
"""
import argparse
import glob
import os
import pickle
import sys
import time

import numpy as np
import uproot

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

COLS = ["run", "subrun", "spill_time_sec", "spill_time_nsec", "TOR860", "TOR875", "FOM",
        "FOM_legacy", "FOM_best", "FOM_fit", "FOM_prefit", "FOM_nomw", "fom_source",
        "fom_status", "hpos_tgt", "vpos_tgt", "hang_tgt", "vang_tgt", "tgt_sigx", "tgt_sigy",
        "triggered", "HPTG1", "VP873"]
CUTS = [0.90, 0.95, 0.98, 0.99]
VARIANTS = ["old", "new", "old_nom", "new_nom"]      # measured only / +100 removed
E_LIN = np.linspace(0, 1, 101)
E_ZOOM = np.linspace(0.9, 1, 101)
E_LOG = np.linspace(-8, 0, 161)                      # log10(1 - FOM)
E_D = np.linspace(-16, 0, 161)                       # log10 |new - old|
RMAX, DAY0, NDAY = 40000, 18500, 3000
BITS = {"NO_TOR": 1, "NO_HBPM": 2, "NO_VBPM": 4, "OFFSET_FILLED": 8, "LOW_TOR": 16,
        "MW_TIME_CUT": 32, "NO_MW_WIDTH": 64}
BEAM = {"hpos_tgt": (-3, 3, "horizontal position at target [mm]"),
        "vpos_tgt": (-3, 3, "vertical position at target [mm]"),
        "hang_tgt": (-1.5, 1.5, "horizontal angle [mrad]"),
        "vang_tgt": (-1.5, 1.5, "vertical angle [mrad]"),
        "tgt_sigx": (0, 5, "horizontal width at target [mm]"),
        "tgt_sigy": (0, 5, "vertical width at target [mm]"),
        "TOR860": (0, 6e12, "TOR860 [protons]"),
        "l1m": (-8, 0, "log10(1 - FOM_best)")}
GROUPS = ["always good", "recovered via filled BPM offset", "recovered, other reason", "was nominal → measured"]


def keyhash(run, sub, sec, nsec):
    t = sec.astype(np.int64) * 1_000_000_000 + nsec.astype(np.int64)
    rs = (run.astype(np.uint64) << np.uint64(20)) + sub.astype(np.uint64)
    with np.errstate(over="ignore"):
        return t.astype(np.uint64) ^ (rs * np.uint64(0x9E3779B97F4A7C15))


def l1m(v):
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.log10(np.clip(1.0 - v, 1e-8, 1.0))


def hist(v, edges, w=None):
    """histogram with under/overflow folded into the edge bins, NaN dropped"""
    m = np.isfinite(v)
    v = np.clip(v[m], edges[0], edges[-1])
    return np.histogram(v, edges, weights=None if w is None else w[m])[0].astype(float)


def new_acc():
    z = lambda *s: np.zeros(s)
    a = {"n": 0, "pot": 0.0, "batches": 0,
         "mig_n": z(3, 3), "mig_pot": z(3, 3), "src_n": z(7), "src_pot": z(7),
         "clos_real": 0, "clos_match": 0, "clos_999": 0, "clos_999_match": 0,
         "clos_other": {}, "clos_bad": [],
         "sanity": {"best_outside_01": 0, "valid_src_but_nan": 0, "none_src_but_value": 0,
                    "fom_without_beam": 0, "pot_le0": 0},
         "cut_pot": {v: z(len(CUTS)) for v in VARIANTS},
         "h_lin": {v: z(100) for v in VARIANTS}, "h_zoom": {v: z(100) for v in VARIANTS},
         "h_log": {v: z(160) for v in VARIANTS}, "h_logpot": {v: z(160) for v in VARIANTS},
         "h2": {k: z(80, 80) for k in ("meas", "nom")},
         "dpos": {k: z(160) for k in ("meas", "nom")}, "dneg": {k: z(160) for k in ("meas", "nom")},
         "dzero": {k: 0 for k in ("meas", "nom")},
         "run_pot": z(RMAX), "run_n": z(RMAX), "run_pass": {v: z(RMAX) for v in VARIANTS},
         "day_pot": z(NDAY), "day_pass": {v: z(NDAY) for v in VARIANTS},
         "bits": {g: z(len(BITS)) for g in ("recovered", "still", "all")},
         "bits_pot": {g: z(len(BITS)) for g in ("recovered", "still", "all")},
         "nset": {g: 0 for g in ("recovered", "still", "all")},
         "beam": {g: {k: z(100) for k in BEAM} for g in GROUPS}, "gsize": {g: 0 for g in GROUPS},
         "nomfit_pos": z(160), "nomfit_neg": z(160), "nomfit_zero": 0, "nomfit_n": 0,
         "trig_n": 0, "trig_found": 0, "trig_eq_new": 0, "trig_eq_old": 0, "trig_bad": []}
    return a


def process(d, a, cut, spill, hashdir, b):
    old = d["FOM"].astype(float)
    leg = d["FOM_legacy"].astype(float)
    new = d["FOM_best"].astype(float)
    src = d["fom_source"].astype(int)
    st = d["fom_status"].astype(np.int64)
    t860, t875 = d["TOR860"].astype(float), d["TOR875"].astype(float)
    tor = np.where(t860 > 0, t860, np.where(t875 > 0, t875, 0.0))
    pot = np.where(np.isfinite(tor), tor, 0.0)
    oc = np.where((old > 0) & (old <= 1), 0, np.where(old > 100, 1, 2))
    nc = np.where((src >= 0) & (src <= 4), 0, np.where(src == 5, 1, 2))
    val = {"old": np.where(oc == 0, old, np.nan), "new": np.where(nc == 0, new, np.nan),
           "old_nom": np.where(oc == 0, old, np.where(oc == 1, old - 100., np.nan)),
           "new_nom": np.where(nc <= 1, new, np.nan)}
    n = len(old)
    a["n"] += n
    a["pot"] += pot.sum()
    np.add.at(a["mig_n"], (oc, nc), 1)
    np.add.at(a["mig_pot"], (oc, nc), pot)
    a["src_n"] += np.bincount(src + 1, minlength=7)[:7]
    a["src_pot"] += np.bincount(src + 1, weights=pot, minlength=7)[:7]

    # closure
    real = oc <= 1
    m = real & np.isclose(old, leg, rtol=0, atol=1e-4)
    a["clos_real"] += int(real.sum())
    a["clos_match"] += int(m.sum())
    bad = np.where(real & ~m)[0]
    for i in bad[:max(0, 200 - len(a["clos_bad"]))]:
        a["clos_bad"].append((int(d["run"][i]), int(d["subrun"][i]), int(d["spill_time_sec"][i]),
                              int(d["spill_time_nsec"][i]), float(old[i]), float(leg[i])))
    c999 = old == -999
    a["clos_999"] += int(c999.sum())
    a["clos_999_match"] += int((c999 & (leg == -999)).sum())
    other = ~real & ~c999
    if other.any():
        lc = np.where((leg > 0) & (leg <= 1), "measured", np.where(leg > 100, "nominal", "none"))
        for code in np.unique(old[other]):
            s = other & (old == code)
            for k in ("measured", "nominal", "none"):
                key = f"{code:g} -> legacy {k}"
                a["clos_other"][key] = a["clos_other"].get(key, 0) + int((s & (lc == k)).sum())

    s = a["sanity"]
    fin = np.isfinite(new)
    s["best_outside_01"] += int((fin & ((new < 0) | (new > 1))).sum())
    s["valid_src_but_nan"] += int(((nc <= 1) & ~fin).sum())
    s["none_src_but_value"] += int(((nc == 2) & fin).sum())
    s["fom_without_beam"] += int((fin & (pot <= 0)).sum())
    s["pot_le0"] += int((pot <= 0).sum())

    for v in VARIANTS:
        x = val[v]
        xx = np.nan_to_num(x, nan=-1.0)
        a["cut_pot"][v] += np.array([pot[xx > c].sum() for c in CUTS])
        a["h_lin"][v] += hist(x, E_LIN)
        a["h_zoom"][v] += hist(x, E_ZOOM)
        a["h_log"][v] += hist(l1m(x), E_LOG)
        a["h_logpot"][v] += hist(l1m(x), E_LOG, pot)
    for k, (ov, nv) in {"meas": ("old", "new"), "nom": ("old_nom", "new_nom")}.items():
        both = np.isfinite(val[ov]) & np.isfinite(val[nv])
        o, w = val[ov][both], val[nv][both]
        a["h2"][k] += np.histogram2d(np.clip(l1m(o), -8, 0), np.clip(l1m(w), -8, 0), [E_LOG[::2], E_LOG[::2]])[0]
        dl = w - o
        a["dzero"][k] += int((dl == 0).sum())
        with np.errstate(divide="ignore"):
            a["dpos"][k] += hist(np.log10(dl[dl > 0]), E_D)
            a["dneg"][k] += hist(np.log10(-dl[dl < 0]), E_D)

    run = d["run"].astype(np.int64)
    if run.max() >= RMAX:
        sys.exit(f"run {run.max()} >= RMAX")
    a["run_pot"] += np.bincount(run, weights=pot, minlength=RMAX)
    a["run_n"] += np.bincount(run, minlength=RMAX)
    day = np.clip(d["spill_time_sec"].astype(np.int64) // 86400 - DAY0, 0, NDAY - 1)
    a["day_pot"] += np.bincount(day, weights=pot, minlength=NDAY)
    for v in VARIANTS:
        ok = np.nan_to_num(val[v], nan=-1.0) > cut
        a["run_pass"][v] += np.bincount(run, weights=pot * ok, minlength=RMAX)
        a["day_pass"][v] += np.bincount(day, weights=pot * ok, minlength=NDAY)

    sets = {"recovered": (oc == 2) & (nc != 2), "still": nc == 2, "all": np.ones(n, bool)}
    for g, sel in sets.items():
        a["nset"][g] += int(sel.sum())
        for j, bit in enumerate(BITS.values()):
            hb = sel & ((st & bit) != 0)
            a["bits"][g][j] += hb.sum()
            a["bits_pot"][g][j] += pot[hb].sum()

    gsel = {GROUPS[0]: (oc == 0) & (nc == 0) & ((st & 8) == 0),
            GROUPS[1]: (oc == 2) & (nc == 0) & ((st & 8) != 0),
            GROUPS[2]: (oc == 2) & (nc == 0) & ((st & 8) == 0),
            GROUPS[3]: (oc == 1) & (nc == 0)}
    # second-BPM fallback (HPTG2 / VPTG2 used because HPTG1 / VP873 has no reading)
    miss = lambda k: ~(np.isfinite(d[k].astype(float)) & (d[k] != -999))
    for k, bpm in (("fb_h", "HPTG1"), ("fb_v", "VP873")):
        if bpm in d:
            m = fin & miss(bpm) & ((st & 128) == 0)      # neighbor-filled spills do not use the fallback
            a.setdefault(k, 0)
            a.setdefault(k + "_pot", 0.0)
            a[k] += int(m.sum())
            a[k + "_pot"] += float(pot[m].sum())
    vals = {k: d[k].astype(float) for k in BEAM if k in d}
    vals["l1m"] = l1m(new)
    for g, sel in gsel.items():
        a["gsize"][g] += int(sel.sum())
        for k, (lo, hi, _) in BEAM.items():
            a["beam"][g][k] += hist(vals[k][sel], np.linspace(lo, hi, 101))

    fit = np.where(nc == 0, new, np.nan)
    nom = d["FOM_nomw"].astype(float)
    both = np.isfinite(fit) & np.isfinite(nom)
    dl = nom[both] - fit[both]
    a["nomfit_n"] += int(both.sum())
    a["nomfit_zero"] += int((dl == 0).sum())
    with np.errstate(divide="ignore"):
        a["nomfit_pos"] += hist(np.log10(dl[dl > 0]), E_D)
        a["nomfit_neg"] += hist(np.log10(-dl[dl < 0]), E_D)

    h = keyhash(d["run"], d["subrun"], d["spill_time_sec"], d["spill_time_nsec"])
    np.save(os.path.join(hashdir, f"h{b:05d}.npy"), h)
    np.save(os.path.join(hashdir, f"p{b:05d}.npy"), pot.astype(np.float32))

    if spill is not None:
        trg = d["triggered"].astype(bool)
        ht = h[trg]
        i = np.clip(np.searchsorted(spill["h"], ht), 0, len(spill["h"]) - 1)
        f = spill["h"][i] == ht
        a["trig_n"] += int(trg.sum())
        a["trig_found"] += int(f.sum())
        nb, ns = new[trg][f], spill["best"][i[f]]
        eqn = (np.isnan(nb) & np.isnan(ns)) | (np.abs(nb - ns) < 1e-9)
        eqo = np.abs(old[trg][f] - spill["fom"][i[f]]) < 1e-6
        a["trig_eq_new"] += int(eqn.sum())
        a["trig_eq_old"] += int(eqo.sum())
        for j in np.where(~eqn)[0][:max(0, 50 - len(a["trig_bad"]))]:
            a["trig_bad"].append((float(nb[j]), float(ns[j])))


def load_spill(path):
    t = uproot.open(path)["spill"]
    want = ["run", "subrun", "spill_time_sec", "spill_time_nsec", "FOM", "FOM_legacy", "FOM_best", "has_spill_info"]
    d = t.arrays([w for w in want if w in t.keys()], library="np")
    s = d["has_spill_info"].astype(bool) if "has_spill_info" in d else np.ones(len(d["FOM"]), bool)
    h = keyhash(d["run"][s], d["subrun"][s], d["spill_time_sec"][s], d["spill_time_nsec"][s])
    o = np.argsort(h)
    with np.errstate(invalid="ignore"):
        old, leg = d["FOM"][s].astype(float), d["FOM_legacy"][s].astype(float)
    real = ((old > 0) & (old <= 1)) | (old > 100)
    info = {"n_records": len(s), "n_with_spill": int(s.sum()), "closure_real": int(real.sum()),
            "closure_match": int((real & np.isclose(old, leg, rtol=0, atol=1e-4)).sum())}
    return {"h": h[o], "best": d["FOM_best"][s].astype(float)[o], "fom": old[o]}, info


# --------------------------------------------------------------------------- plots
def plots(a, label, cut, outdir, tag):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.colors import LogNorm
    import plot_fom_compare as P                       # same palette / style / names
    OLD, NEW, INK, INK2 = P.OLD, P.NEW, P.INK, P.INK2
    figs = []
    ctr = lambda e: 0.5 * (e[1:] + e[:-1])

    # 1 accounting
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8))
    for ax, M, unit in ((axes[0], a["mig_n"], "spills"), (axes[1], a["mig_pot"], "POT")):
        oc, nc = M.sum(1), M.sum(0)
        y = np.arange(3)
        ax.barh(y - 0.19, oc, height=0.36, color=OLD, label="old (CAF)")
        ax.barh(y + 0.19, nc, height=0.36, color=NEW, label="new (fixed)")
        tot = oc.sum()
        for yy, o, nn in zip(y, oc, nc):
            ax.text(max(o, nn) * 1.01, yy, f"{o/tot:.1%} → {nn/tot:.1%}", va="center", color=INK2, fontsize=9)
        ax.set_yticks(y, P.CATS)
        ax.invert_yaxis()
        ax.set_xlabel(unit)
        ax.set_xlim(0, max(oc.max(), nc.max()) * 1.32)
        ax.grid(axis="y", visible=False)
    axes[0].set_title(f"{label}: where spills land")
    axes[1].set_title("same, weighted by POT")
    axes[0].legend(loc="lower right")
    fig.tight_layout()
    figs.append(("common", "accounting", fig))

    # 2 migration
    m, mp = a["mig_n"], a["mig_pot"]
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    norm = LogNorm(vmin=max(1, m[m > 0].min()), vmax=max(m.max(), 2))
    ax.imshow(np.where(m > 0, m, np.nan), cmap=P.SEQ_CMAP, norm=norm)
    for i in range(3):
        for j in range(3):
            dark = m[i, j] > 0 and norm(m[i, j]) > 0.5
            ax.text(j, i, f"{int(m[i,j]):,}\n{m[i,j]/m.sum():.1%} spills\n{mp[i,j]/max(mp.sum(),1):.1%} POT",
                    ha="center", va="center", fontsize=9, color="#ffffff" if dark else INK)
    ax.set_xticks(range(3), P.CATS)
    ax.set_yticks(range(3), P.CATS)
    ax.set_xlabel("new (fixed)")
    ax.set_ylabel("old (CAF)")
    ax.grid(False)
    ax.set_title(f"{label}: category migration old → new")
    fig.tight_layout()
    figs.append(("common", "migration", fig))

    # 3 status bits
    names = list(BITS)
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.9))
    for ax, g, col, title in ((axes[0], "recovered", NEW, "recovered (old: no FOM → new: FOM)"),
                              (axes[1], "still", OLD, "still without a FOM")):
        cnt = a["bits"][g]
        y = np.arange(len(names))
        ax.barh(y, cnt, height=0.62, color=col)
        for yy, c in zip(y, cnt):
            if c:
                ax.text(c, yy, f" {int(c):,}", va="center", fontsize=9, color=INK2)
        ax.set_yticks(y, [P.STATUS_TEXT[k] for k in names])
        ax.invert_yaxis()
        ax.set_xlim(0, max(cnt.max(), 1) * 1.25)
        ax.set_xlabel("spills (a spill can have several flags)")
        ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
        ax.set_title(f"{title}  n={a['nset'][g]:,}", fontsize=10)
        ax.grid(axis="y", visible=False)
    fig.suptitle(f"{label}: why spills changed", fontweight="bold")
    fig.tight_layout()
    figs.append(("common", "status", fig))

    for grp, nom in (("measured_only", False), ("offset_removed", True)):
        ov, nv = ("old_nom", "new_nom") if nom else ("old", "new")
        tag2 = "+100 removed (nominal width included)" if nom else "measured width only"
        # 4 distributions
        fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
        for ax, key, e, xl in ((axes[0], "h_lin", E_LIN, "FOM"), (axes[1], "h_zoom", E_ZOOM, "FOM (zoom)"),
                               (axes[2], "h_log", E_LOG, "log10(1 − FOM)")):
            for v, c, nm in ((ov, OLD, "old (CAF)"), (nv, NEW, "new (fixed)")):
                ax.stairs(a[key][v], e, color=c, lw=2, label=f"{nm}: {a['h_lin'][v].sum():,.0f} spills")
            ax.set_yscale("log")
            ax.set_xlabel(xl)
            ax.set_ylabel("spills / bin")
        axes[0].set_title(f"{label}: FOM, {tag2}", fontsize=10)
        axes[1].set_title("zoom 0.9–1 (edges hold overflow)", fontsize=10)
        axes[2].set_title("loss fraction (1e-8 = FOM rounding to 1)", fontsize=10)
        axes[0].legend(loc="upper left", fontsize=8)
        fig.tight_layout()
        figs.append((grp, "distributions", fig))

        # 5 per-spill change
        k = "nom" if nom else "meas"
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
        H = a["h2"][k]
        e2 = E_LOG[::2]
        im = axes[0].pcolormesh(e2, e2, np.where(H.T > 0, H.T, np.nan), cmap=P.SEQ_CMAP, norm=LogNorm())
        fig.colorbar(im, ax=axes[0], label="spills")
        axes[0].plot([-8, 0], [-8, 0], color=INK2, lw=1, ls="--")
        axes[0].set_xlabel("old: log10(1 − FOM)")
        axes[0].set_ylabel("new: log10(1 − FOM)")
        axes[0].set_title(f"{label}: spills with a FOM in both ({H.sum():,.0f})", fontsize=10)
        axes[0].grid(False)
        c = ctr(E_D)
        axes[1].stairs(a["dpos"][k], E_D, color=NEW, lw=2, label=f"new > old: {a['dpos'][k].sum():,.0f}")
        axes[1].stairs(a["dneg"][k], E_D, color=OLD, lw=2, label=f"new < old: {a['dneg'][k].sum():,.0f}")
        axes[1].set_yscale("log")
        axes[1].set_xlabel("log10 |new − old|")
        axes[1].set_ylabel("spills / bin")
        axes[1].set_title(f"size of the change ({a['dzero'][k]:,} identical)", fontsize=10)
        axes[1].legend(fontsize=8)
        fig.tight_layout()
        figs.append((grp, "value_change", fig))

        # 6 per run / per day
        fig, axes = plt.subplots(2, 2, figsize=(14, 6.4), gridspec_kw={"height_ratios": [2, 1]})
        for col, (tot, pas, xs, xl) in enumerate((
                (a["run_pot"], a["run_pass"], np.arange(RMAX), "run"),
                (a["day_pot"], a["day_pass"], (np.arange(NDAY) + DAY0).astype("datetime64[D]"), "date"))):
            ok = tot > 0
            fo, fn = pas[ov][ok] / tot[ok], pas[nv][ok] / tot[ok]
            axes[0, col].plot(xs[ok], fo, "o", ms=3, color=OLD, label=f"old: {pas[ov].sum()/tot.sum():.1%} of POT")
            axes[0, col].plot(xs[ok], fn, "o", ms=3, color=NEW, label=f"new: {pas[nv].sum()/tot.sum():.1%} of POT")
            axes[0, col].set_ylim(-0.02, 1.02)
            axes[0, col].set_ylabel(f"POT fraction FOM > {cut}")
            axes[0, col].legend(loc="lower left", fontsize=8)
            axes[1, col].plot(xs[ok], fn - fo, "o", ms=3, color=INK2)
            axes[1, col].axhline(0, color=INK2, lw=1)
            axes[1, col].set_ylabel("new − old")
            axes[1, col].set_xlabel(xl)
        axes[0, 0].set_title(f"{label}: per run, {tag2}", fontsize=10)
        axes[0, 1].set_title("per day", fontsize=10)
        fig.autofmt_xdate()
        fig.tight_layout()
        figs.append((grp, "pot_per_run_day", fig))

    # 7 recovered vs always-good beam
    fig, axes = plt.subplots(2, 4, figsize=(16, 7))
    colors = [P.INK2, NEW, OLD, "#5ba300"]
    for ax, (k, (lo, hi, xl)) in zip(axes.flat, BEAM.items()):
        e = np.linspace(lo, hi, 101)
        for g, c in zip(GROUPS, colors):
            h = a["beam"][g][k]
            if a["gsize"][g] >= 1000:
                ax.stairs(h / h.sum(), e, color=c, lw=1.8, label=f"{g} ({a['gsize'][g]:,})")
        ax.set_yscale("log")
        ax.set_xlabel(xl)
    axes[0, 0].set_ylabel("fraction / bin")
    axes[1, 0].set_ylabel("fraction / bin")
    axes[0, 0].legend(fontsize=7, loc="upper left")
    fig.suptitle(f"{label}: do recovered spills look like spills that were always good? (edges hold overflow)",
                 fontweight="bold")
    fig.tight_layout()
    figs.append(("validation", "recovered_vs_good", fig))

    # 8 nominal vs measured
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.stairs(a["nomfit_pos"], E_D, color=NEW, lw=2, label=f"nominal > measured: {a['nomfit_pos'].sum():,.0f}")
    ax.stairs(a["nomfit_neg"], E_D, color=OLD, lw=2, label=f"nominal < measured: {a['nomfit_neg'].sum():,.0f}")
    ax.set_yscale("log")
    ax.set_xlabel("log10 |FOM_nomw − FOM(measured width)|")
    ax.set_ylabel("spills / bin")
    ax.set_title(f"{label}: nominal-width vs measured-width FOM, same spill ({a['nomfit_n']:,})", fontsize=10)
    ax.legend(fontsize=8)
    fig.tight_layout()
    figs.append(("validation", "nominal_vs_measured", fig))

    with PdfPages(os.path.join(outdir, f"fom_validation_{tag}.pdf")) as pdf:
        for grp, name, fig in figs:
            os.makedirs(os.path.join(outdir, tag, grp), exist_ok=True)
            fig.savefig(os.path.join(outdir, tag, grp, f"{tag}_{name}.png"), dpi=150)
            pdf.savefig(fig)
            plt.close(fig)


def dup_check(hashdir):
    hs = sorted(glob.glob(os.path.join(hashdir, "h*.npy")))
    h = np.concatenate([np.load(f) for f in hs])
    p = np.concatenate([np.load(os.path.join(os.path.dirname(f), "p" + os.path.basename(f)[1:])) for f in hs])
    o = np.argsort(h, kind="stable")
    h, p = h[o], p[o]
    rep = np.zeros(len(h), bool)
    rep[1:] = h[1:] == h[:-1]
    return int(rep.sum()), float(p[rep].sum())


def report(a, info, dups, label, cut):
    P = a["pot"]
    L = [f"## {label}", "",
         f"{a['n']:,} subrun spills (bnbinfo), {P:.4g} POT; spill tree: {info['n_records']:,} records, "
         f"{info['n_with_spill']:,} with spill info.", "", "### Validation", "",
         "| check | result | |", "|---|---|---|"]
    ok = lambda c: "PASS" if c else "**CHECK**"
    cr = a["clos_match"] / max(a["clos_real"], 1)
    L.append(f"| closure, CAF FOM real value (measured or 100+nominal) | {a['clos_match']:,}/{a['clos_real']:,} = {cr:.5%} | {ok(cr > 0.9999)} |")
    c9 = a["clos_999_match"] / max(a["clos_999"], 1)
    L.append(f"| closure, CAF FOM = -999 → legacy -999 | {a['clos_999_match']:,}/{a['clos_999']:,} = {c9:.4%} | {ok(c9 > 0.999 or a['clos_999'] == 0)} |")
    if a["clos_other"]:
        L.append("| other CAF codes (older FillExposure) | " +
                 "; ".join(f"{k}: {v:,}" for k, v in sorted(a["clos_other"].items()) if v) + " | info |")
    sc = info["closure_match"] / max(info["closure_real"], 1)
    L.append(f"| closure, spill tree (real CAF values) | {info['closure_match']:,}/{info['closure_real']:,} = {sc:.5%} | {ok(sc > 0.9999)} |")
    s = a["sanity"]
    L.append(f"| FOM_best outside [0,1] | {s['best_outside_01']:,} | {ok(s['best_outside_01'] == 0)} |")
    L.append(f"| fom_source valid but FOM_best NaN / none but value | {s['valid_src_but_nan']:,} / {s['none_src_but_value']:,} | "
             f"{ok(s['valid_src_but_nan'] == 0 and s['none_src_but_value'] == 0)} |")
    L.append(f"| new FOM on a spill without beam (POT ≤ 0) | {s['fom_without_beam']:,} (of {s['pot_le0']:,} beam-off spills) | {ok(s['fom_without_beam'] == 0)} |")
    L.append(f"| duplicate (run, subrun, spill time) entries | {dups[0]:,} spills, {dups[1]:.3g} POT ({dups[1]/P:.3%}) | {ok(dups[0] / a['n'] < 1e-4)} |")
    tf = a["trig_found"] / max(a["trig_n"], 1)
    te = a["trig_eq_new"] / max(a["trig_found"], 1)
    L.append(f"| triggered bnbinfo spills found in spill tree | {a['trig_found']:,}/{a['trig_n']:,} = {tf:.4%} | {ok(tf > 0.999)} |")
    L.append(f"| same FOM_best in spill and bnbinfo trees | {te:.5%} (CAF FOM same: {a['trig_eq_old']/max(a['trig_found'],1):.5%}) | {ok(te > 0.9999)} |")
    if "fb_v" in a:
        L.append(f"| FOM from the second-BPM fallback (HPTG2/VPTG2; offsets not calibrated, biased position) | "
                 f"H {a['fb_h']:,}, V {a['fb_v']:,} spills ({(a['fb_h_pot'] + a['fb_v_pot'])/P:.3%} POT) | "
                 f"{'info' if (a['fb_h_pot'] + a['fb_v_pot'])/P < 1e-3 else '**CHECK**'} |")
    if a["clos_bad"]:
        L += ["", f"First closure mismatches (run, subrun, sec, nsec, CAF, legacy): {a['clos_bad'][:5]}"]

    L += ["", "### Old vs new", "", "| category | old spills | new spills | old POT | new POT |", "|---|---:|---:|---:|---:|"]
    m, mp = a["mig_n"], a["mig_pot"]
    for c, name in enumerate(["measured width", "nominal width", "no FOM"]):
        L.append(f"| {name} | {m[c].sum()/m.sum():.2%} | {m[:, c].sum()/m.sum():.2%} | {mp[c].sum()/P:.2%} | {mp[:, c].sum()/P:.2%} |")
    L += ["", "POT fraction passing the cut:", "",
          "| cut | measured only: old | measured only: new | +100 removed: old | +100 removed: new |", "|---|---:|---:|---:|---:|"]
    for i, c in enumerate(CUTS):
        L.append(f"| FOM > {c} | " + " | ".join(f"{a['cut_pot'][v][i]/P:.2%}" for v in VARIANTS) + " |")
    src_names = ["none", "fit_MMBTBB", "fit_M876BB", "fit_M875BB", "prefit_M876", "prefit_M875", "no_multiwire"]
    L += ["", "New FOM source: " + ", ".join(f"{n} {a['src_pot'][i]/P:.2%} POT" for i, n in enumerate(src_names) if a["src_n"][i]), ""]
    rec = a["nset"]["recovered"]
    L.append(f"Recovered spills (old no FOM → new FOM): {rec:,}; still without FOM: {a['nset']['still']:,}. "
             "Flags on recovered spills: " + ", ".join(f"{k} {int(c):,}" for k, c in zip(BITS, a["bits"]["recovered"]) if c) + ".")
    L.append("")
    L.append("POT still without a FOM, by flag (a spill can carry several): " +
             ", ".join(f"{k} {c/P:.2%}" for k, c in zip(BITS, a["bits_pot"]["still"]) if c) + ".")
    L.append("")
    return "\n".join(L)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help="bnb_fom(_fast).py output")
    ap.add_argument("--label", required=True)
    ap.add_argument("-o", "--outdir", default="fom_validation")
    ap.add_argument("--cut", type=float, default=0.95)
    ap.add_argument("--step", type=int, default=2_000_000, help="spills per batch")
    ap.add_argument("--max-seconds", type=float, default=0, help="stop (checkpointed) after this long")
    ap.add_argument("--restart", action="store_true")
    a = ap.parse_args(argv)
    tag = a.label.lower().replace(" ", "_").replace("(", "").replace(")", "")
    wd = os.path.join(a.outdir, tag)
    hashdir = os.path.join(wd, "hashes")
    os.makedirs(hashdir, exist_ok=True)
    sp = os.path.join(wd, "state.pkl")
    spf = os.path.join(wd, "spill.pkl")
    if a.restart:
        for f in glob.glob(os.path.join(hashdir, "*.npy")) + [sp, spf]:
            if os.path.exists(f):
                os.remove(f)
    if os.path.exists(sp):
        with open(sp, "rb") as fh:
            acc = pickle.load(fh)
    else:
        acc = new_acc()
    if os.path.exists(spf):
        with open(spf, "rb") as fh:
            spill, info = pickle.load(fh)
    else:
        spill, info = load_spill(a.input)
        with open(spf, "wb") as fh:
            pickle.dump((spill, info), fh)
    t = uproot.open(a.input)["bnbinfo"]
    N = t.num_entries
    keys = set(t.keys())
    cols = [c for c in COLS if c in keys]
    t0 = time.time()
    nb = (N + a.step - 1) // a.step
    for b in range(acc["batches"], nb):
        d = t.arrays(cols, entry_start=b * a.step, entry_stop=min(N, (b + 1) * a.step), library="np")
        process(d, acc, a.cut, spill, hashdir, b)
        acc["batches"] = b + 1
        with open(sp + ".tmp", "wb") as fh:
            pickle.dump(acc, fh)
        os.replace(sp + ".tmp", sp)
        print(f"\r[{tag}] batch {b + 1}/{nb}  {time.time() - t0:.0f}s", end="", file=sys.stderr, flush=True)
        if a.max_seconds and time.time() - t0 > a.max_seconds and b + 1 < nb:
            print(f"\n[{tag}] paused after batch {b + 1}/{nb}; rerun the same command to continue", file=sys.stderr)
            return 3
    print(file=sys.stderr)
    dups = dup_check(hashdir)
    plots(acc, a.label, a.cut, a.outdir, tag)
    txt = report(acc, info, dups, a.label, a.cut)
    with open(os.path.join(wd, "report.md"), "w") as fh:
        fh.write(txt + "\n")
    print(txt)
    return 0


if __name__ == "__main__":
    sys.exit(main())
