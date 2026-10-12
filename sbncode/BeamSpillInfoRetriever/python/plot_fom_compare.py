#!/usr/bin/env python3
"""
plot_fom_compare.py -- compare the production (CAF) BNB FOM with the fixed FOM
from bnb_fom.py.

  ./plot_fom_compare.py beamqual_run2_fom.root beamqual_run4_fom.root \
        --labels "Run 2" "Run 4" -o fom_plots --cut 0.95

Writes one PNG per figure per input plus a combined PDF into the output folder,
and prints a summary table.

Definitions
  old  = rec.hdr.bnbinfo.FOM as stored in the CAF (production code)
  new  = FOM_best from bnb_fom.py (fixed code)
  Old FOM categories (FillExposure.cxx convention):
     measured width : FOM in [0,1]  (fitted or database multiwire width)
     nominal width  : FOM > 100     (100 + FOM with nominal width)
     no FOM         : FOM = -999
  New categories from fom_source: fit_* / prefit_* = measured width,
     no_multiwire = nominal width, none = no FOM.
  POT = TOR860 (TOR875 if TOR860 missing), only spills with positive intensity.
  FOM-cut plots are made twice: measured-width FOMs only (nominal-width spills
  fail the cut, as with a plain `FOM > cut` on the CAF value), and with the
  nominal-width FOMs included (100 subtracted from the CAF value).
"""
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.ticker
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap, LogNorm
import uproot

# palette (validated reference palette, light surface)
OLD, NEW = "#2a78d6", "#eb6834"                       # categorical slots 1, 2
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
SEQ = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
SEQ_CMAP = LinearSegmentedColormap.from_list("seqblue", SEQ)
STATUS_BITS = {"NO_TOR": 1, "NO_HBPM": 2, "NO_VBPM": 4, "OFFSET_FILLED": 8,
               "LOW_TOR": 16, "MW_TIME_CUT": 32, "NO_MW_WIDTH": 64}
STATUS_TEXT = {"NO_TOR": "no intensity (TOR ≤ 0 / missing)", "NO_HBPM": "horizontal BPM missing",
               "NO_VBPM": "vertical BPM missing", "OFFSET_FILLED": "BPM offset filled",
               "LOW_TOR": "low intensity (< 1e11)", "MW_TIME_CUT": "multiwire profile > 33 ms away",
               "NO_MW_WIDTH": "no measured width (nominal used)"}
CATS = ["measured width", "nominal width", "no FOM"]

plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False, "font.size": 10,
    "axes.titlesize": 11, "axes.titleweight": "bold", "legend.frameon": False,
    "lines.linewidth": 2,
})


# beam-quality selection (--beam-cuts): intensity and horn current; the FOM window
# 0.98 < FOM <= 1 is the --cut (FOM values used here are never > 1).
# LM875A/B/C > 1e-2 is part of the official selection but the loss monitors are not in the dumps.
BEAM_CUTS = dict(tor_min=1e11, thcurr=(173.0, 175.0))


def beam_mask(d):
    t860, t875, hc = d["TOR860"], d["TOR875"], d["THCURR"]
    with np.errstate(invalid="ignore"):
        return ((t860 > BEAM_CUTS["tor_min"]) & (t875 > BEAM_CUTS["tor_min"])
                & (hc >= BEAM_CUTS["thcurr"][0]) & (hc <= BEAM_CUTS["thcurr"][1]))


def load(path, beam_cuts=False):
    t = uproot.open(path)["bnbinfo"]
    want = ["run", "FOM", "FOM_legacy", "FOM_best", "fom_source", "fom_status", "TOR860", "TOR875", "THCURR"]
    d = t.arrays([w for w in want if w in t.keys()], library="np")
    if beam_cuts:
        m = beam_mask(d)
        d = {k: v[m] for k, v in d.items()}
    old = d["FOM"].astype(float)
    new = d["FOM_best"].astype(float)
    src = d["fom_source"].astype(int)
    # production never stores a genuine 0 (FillExposure needs FOM > 0); older productions
    # (e.g. SBND run1dev) store failure codes 0 (BPM missing) and 4 (no multiwire) instead of -999
    d["old_cat"] = np.where((old > 0) & (old <= 1), 0, np.where(old > 100, 1, 2))
    d["new_cat"] = np.where((src >= 0) & (src <= 4), 0, np.where(src == 5, 1, 2))
    tor = np.where(d["TOR860"] > 0, d["TOR860"], np.where(d.get("TOR875", d["TOR860"]) > 0, d.get("TOR875", d["TOR860"]), 0.0))
    d["pot"] = np.where(np.isfinite(tor), tor, 0.0)
    d["old_val"] = np.where(d["old_cat"] == 0, old, np.nan)          # usable old FOM
    d["new_val"] = np.where(d["new_cat"] == 0, new, np.nan)          # usable new FOM
    # same, also counting nominal-width FOMs (CAF stores them as 100 + FOM; subtract the 100)
    d["old_val_nom"] = np.where(d["old_cat"] == 0, old, np.where(d["old_cat"] == 1, old - 100., np.nan))
    d["new_val_nom"] = np.where(d["new_cat"] <= 1, new, np.nan)      # FOM_best has no +100
    return d


def fmt_pot(x):
    return f"{x:.3g}"


def fig_distributions(d, label, nominal=False):
    ov, nv = ("old_val_nom", "new_val_nom") if nominal else ("old_val", "new_val")
    tag = "+100 offset removed (nominal width included)" if nominal else "measured width only"
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, (lo, hi, nb) in zip(axes, [(0, 1, 100), (0.9, 1, 100)]):
        bins = np.linspace(lo, hi, nb + 1)
        for v, c, name in ((d[ov], OLD, "old (CAF)"), (d[nv], NEW, "new (fixed)")):
            v = v[np.isfinite(v)]
            ax.hist(np.clip(v, lo, hi), bins=bins, histtype="step", color=c, lw=2,
                    label=f"{name}: {len(v):,} spills")
        ax.set_yscale("log")
        ax.set_xlabel("FOM")
        ax.set_ylabel("spills / bin")
        ax.set_xlim(lo, hi)
    axes[0].set_title(f"{label}: FOM, {tag}", fontsize=10)
    axes[1].set_title("zoom 0.9 – 1 (overflow clipped into edge bins)")
    axes[0].legend(loc="upper left")
    fig.tight_layout()
    return fig


def fig_accounting(d, label):
    """old vs new spill and POT counts in the three categories."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8))
    for ax, weight, unit in ((axes[0], None, "spills"), (axes[1], d["pot"], "POT (TOR860)")):
        oc = [np.sum((d["old_cat"] == k) * (1 if weight is None else weight)) for k in range(3)]
        nc = [np.sum((d["new_cat"] == k) * (1 if weight is None else weight)) for k in range(3)]
        y = np.arange(3)
        h = 0.36
        ax.barh(y - h / 2 - 0.01, oc, height=h, color=OLD, label="old (CAF)")
        ax.barh(y + h / 2 + 0.01, nc, height=h, color=NEW, label="new (fixed)")
        tot = sum(oc)
        for yy, o, n in zip(y, oc, nc):
            ax.text(max(o, n) * 1.01, yy, f"{o/tot:.1%} → {n/tot:.1%}", va="center", color=INK2, fontsize=9)
        ax.set_yticks(y, CATS)
        ax.invert_yaxis()
        ax.set_xlabel(unit)
        ax.set_xlim(0, max(max(oc), max(nc)) * 1.32)
        ax.grid(axis="y", visible=False)
    axes[0].set_title(f"{label}: where spills land")
    axes[1].set_title("same, weighted by POT")
    axes[0].legend(loc="lower right")
    fig.tight_layout()
    return fig


def fig_migration(d, label):
    m = np.zeros((3, 3))
    mp = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            sel = (d["old_cat"] == i) & (d["new_cat"] == j)
            m[i, j] = sel.sum()
            mp[i, j] = d["pot"][sel].sum()
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    norm = LogNorm(vmin=max(1, m[m > 0].min()), vmax=max(m.max(), 2))
    ax.imshow(np.where(m > 0, m, np.nan), cmap=SEQ_CMAP, norm=norm)
    tot, totp = m.sum(), mp.sum()
    for i in range(3):
        for j in range(3):
            frac = m[i, j] / tot
            dark = m[i, j] > 0 and norm(m[i, j]) > 0.5      # white text only on the darker half of the ramp
            ax.text(j, i, f"{int(m[i,j]):,}\n{frac:.1%} spills\n{mp[i,j]/max(totp,1):.1%} POT",
                    ha="center", va="center", fontsize=9, color="#ffffff" if dark else INK)
    ax.set_xticks(range(3), CATS)
    ax.set_yticks(range(3), CATS)
    ax.set_xlabel("new (fixed)")
    ax.set_ylabel("old (CAF)")
    ax.grid(False)
    ax.set_title(f"{label}: category migration old → new")
    fig.tight_layout()
    return fig


def fig_value_change(d, label, nominal=False):
    ov, nv = ("old_val_nom", "new_val_nom") if nominal else ("old_val", "new_val")
    both = np.isfinite(d[ov]) & np.isfinite(d[nv])
    o, n = d[ov][both], d[nv][both]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    lo = max(0.0, np.percentile(o, 0.5)) if len(o) else 0
    h = axes[0].hist2d(o, n, bins=120, range=[[lo, 1], [lo, 1]], cmap=SEQ_CMAP, norm=LogNorm(), cmin=1)
    fig.colorbar(h[3], ax=axes[0], label="spills")
    axes[0].plot([lo, 1], [lo, 1], color=INK2, lw=1, ls="--")
    axes[0].set_xlabel("old FOM (CAF)")
    axes[0].set_ylabel("new FOM (fixed)")
    axes[0].set_title(f"{label}: spills with a FOM in both ({len(o):,}), "
                      + ("+100 removed" if nominal else "measured width only"), fontsize=10)
    axes[0].grid(False)
    dlt = n - o
    lim = max(1e-6, np.percentile(np.abs(dlt), 99.5)) if len(dlt) else 1e-3
    axes[1].hist(np.clip(dlt, -lim, lim), bins=120, range=(-lim, lim), color=NEW, histtype="stepfilled", alpha=0.85)
    axes[1].set_yscale("log")
    axes[1].set_xlabel(f"new − old (values beyond ±{lim:.1e} collected in the edge bins)")
    axes[1].set_ylabel("spills / bin")
    if len(dlt):
        axes[1].set_title(f"FOM change: median {np.median(dlt):+.1e}, 99% within ±{np.percentile(np.abs(dlt), 99):.1e}")
    fig.tight_layout()
    return fig


def fig_per_run(d, label, cut, nominal=False):
    ov, nv = ("old_val_nom", "new_val_nom") if nominal else ("old_val", "new_val")
    runs = np.unique(d["run"])
    pot_tot, pot_old, pot_new = [], [], []
    for r in runs:
        s = d["run"] == r
        p = d["pot"][s]
        pot_tot.append(p.sum())
        pot_old.append(p[np.nan_to_num(d[ov][s], nan=-1) > cut].sum())
        pot_new.append(p[np.nan_to_num(d[nv][s], nan=-1) > cut].sum())
    pot_tot, pot_old, pot_new = map(np.asarray, (pot_tot, pot_old, pot_new))
    ok = pot_tot > 0
    fig, axes = plt.subplots(2, 1, figsize=(11, 6.2), sharex=True, gridspec_kw={"height_ratios": [2, 1]})
    axes[0].plot(runs[ok], pot_old[ok] / pot_tot[ok], "o", ms=3.5, color=OLD, label=f"old (CAF): {pot_old.sum()/pot_tot.sum():.1%} of POT")
    axes[0].plot(runs[ok], pot_new[ok] / pot_tot[ok], "o", ms=3.5, color=NEW, label=f"new (fixed): {pot_new.sum()/pot_tot.sum():.1%} of POT")
    axes[0].set_ylabel(f"POT fraction with FOM > {cut}")
    axes[0].set_ylim(-0.02, 1.02)
    axes[0].set_title(f"{label}: POT passing FOM > {cut}, per run "
                      + ("(nominal-width FOMs included, 100 subtracted)" if nominal else "(measured-width FOMs only)"))
    axes[0].legend(loc="lower left")
    axes[1].plot(runs[ok], (pot_new[ok] - pot_old[ok]) / pot_tot[ok], "o", ms=3.5, color=INK2)
    axes[1].axhline(0, color=INK2, lw=1)
    axes[1].set_ylabel("new − old")
    axes[1].set_xlabel("run")
    fig.tight_layout()
    return fig, (pot_tot.sum(), pot_old.sum(), pot_new.sum())


def fig_status(d, label):
    st = d["fom_status"].astype(int)
    oldbad = d["old_cat"] == 2
    rec = oldbad & (d["new_cat"] != 2)
    still = d["new_cat"] == 2
    names = list(STATUS_BITS)
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.9))
    for ax, sel, title in ((axes[0], rec, f"recovered spills (old: no FOM → new: FOM)  n={rec.sum():,}"),
                           (axes[1], still, f"spills still without a FOM  n={still.sum():,}")):
        cnt = [np.sum(sel & ((st & STATUS_BITS[k]) != 0)) for k in names]
        y = np.arange(len(names))
        ax.barh(y, cnt, height=0.62, color=NEW if ax is axes[0] else OLD)
        for yy, c in zip(y, cnt):
            if c:
                ax.text(c, yy, f" {c:,}", va="center", fontsize=9, color=INK2)
        ax.set_yticks(y, [STATUS_TEXT[k] for k in names])
        ax.invert_yaxis()
        ax.set_xlim(0, max(max(cnt), 1) * 1.2)
        ax.set_xlabel("spills (a spill can have several flags)")
        ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
        ax.set_title(title, fontsize=10)
        ax.grid(axis="y", visible=False)
    fig.suptitle(f"{label}: why spills changed", fontweight="bold")
    fig.tight_layout()
    return fig


def pot_pass(d, key, cut):
    v = np.nan_to_num(d[key], nan=-1.0)
    return d["pot"][v > cut].sum() / d["pot"].sum()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="bnb_fom.py outputs (.root)")
    ap.add_argument("--labels", nargs="+", help="one label per input")
    ap.add_argument("-o", "--outdir", default="fom_plots")
    ap.add_argument("--cut", type=float, default=0.95, help="FOM cut for the per-run POT plots")
    ap.add_argument("--beam-cuts", action="store_true",
                    help="keep only spills with TOR860, TOR875 > 1e11 and 173 <= THCURR <= 175 kA "
                         "(denominator = POT of those spills)")
    ap.add_argument("--table-cuts", type=float, nargs="+", default=[0.90, 0.95, 0.98, 0.99])
    a = ap.parse_args(argv)
    labels = a.labels or [os.path.basename(p).replace("_fom.root", "") for p in a.inputs]
    os.makedirs(a.outdir, exist_ok=True)
    lines = ["# BNB FOM: old (CAF) vs new (fixed)", "",
             "POT = TOR860 (TOR875 if TOR860 missing), positive-intensity spills only.",
             "*Measured width only*: FOM in [0,1]; CAF values > 100 (nominal-width FOM) fail every cut.",
             "*+100 offset removed*: CAF values > 100 are counted as FOM - 100 (new FOM has no +100).", ""]
    if a.beam_cuts:
        lines[2:2] = ["Beam cuts applied: TOR860 > 1e11 and TOR875 > 1e11, 173 <= THCURR <= 175 kA "
                      "(LM875A/B/C not in the dumps, not applied). Fractions are of the POT passing these.", ""]
    for path, label in zip(a.inputs, labels):
        d = load(path, a.beam_cuts)
        tag = label.lower().replace(" ", "")
        groups = {
            "common": [("accounting", fig_accounting(d, label)), ("migration", fig_migration(d, label)),
                       ("status", fig_status(d, label))],
            "measured_only": [("distributions", fig_distributions(d, label)),
                              ("value_change", fig_value_change(d, label)),
                              ("pot_per_run", fig_per_run(d, label, a.cut)[0])],
            "offset_removed": [("distributions", fig_distributions(d, label, nominal=True)),
                               ("value_change", fig_value_change(d, label, nominal=True)),
                               ("pot_per_run", fig_per_run(d, label, a.cut, nominal=True)[0])],
        }
        with PdfPages(os.path.join(a.outdir, f"fom_comparison_{tag}.pdf")) as pdf:
            for grp, figs in groups.items():
                os.makedirs(os.path.join(a.outdir, tag, grp), exist_ok=True)
                for name, fig in figs:
                    fig.savefig(os.path.join(a.outdir, tag, grp, f"{tag}_{name}.png"), dpi=150)
                    pdf.savefig(fig)
                    plt.close(fig)

        n, P = len(d["FOM"]), d["pot"].sum()
        closure = np.mean(np.isclose(d["FOM"], d["FOM_legacy"], atol=1e-4)) if "FOM_legacy" in d else float("nan")
        cat = lambda k, c: (np.mean(d[k] == c), d["pot"][d[k] == c].sum() / P)
        lines += [f"## {label}", "",
                  f"{n:,} subrun spills, {P:.4g} POT. Legacy closure (recomputed old FOM = CAF FOM): {closure:.4%}.", "",
                  "| category | old spills | new spills | old POT | new POT |", "|---|---:|---:|---:|---:|"]
        for c, name in enumerate(CATS):
            os_, op = cat("old_cat", c)
            ns_, np_ = cat("new_cat", c)
            lines.append(f"| {name} | {os_:.2%} | {ns_:.2%} | {op:.2%} | {np_:.2%} |")
        lines += ["", "POT fraction passing the cut:", "",
                  "| cut | measured only: old | measured only: new | +100 removed: old | +100 removed: new |",
                  "|---|---:|---:|---:|---:|"]
        for c in a.table_cuts:
            lines.append(f"| FOM > {c} | {pot_pass(d, 'old_val', c):.2%} | {pot_pass(d, 'new_val', c):.2%} | "
                         f"{pot_pass(d, 'old_val_nom', c):.2%} | {pot_pass(d, 'new_val_nom', c):.2%} |")
        lines.append("")
    txt = "\n".join(lines)
    with open(os.path.join(a.outdir, "summary.md"), "w") as fh:
        fh.write(txt + "\n")
    print(txt)
    print(f"\nplots in {a.outdir}/<run>/{{common,measured_only,offset_removed}}/, one PDF per run, summary.md")


if __name__ == "__main__":
    main()
