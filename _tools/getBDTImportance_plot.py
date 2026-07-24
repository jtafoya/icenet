#!/usr/bin/env python
# ===========================================================================
#  getBDTImportance_plot.py  --  readable XGBoost feature-importance plots.
#
#  Remakes the icenet training plot (figs/.../train/xgboost-importance/), which
#  puts all 532 features on one page: the labels are illegible and ~90% of the
#  canvas is the flat tail (191 features are never split on at all).
#
#  Here: a ranked horizontal bar chart of the top N only, coloured by collection,
#  with the cumulative share stated so the truncation is honest. Two pages:
#    1. top N individual features (e.g. muonSV_dxy_0)
#    2. top N base variables, summed over the 8 leading-object slots _0.._7
#       -- this is what you want if the question is "which variable matters",
#       since one variable is spread over up to 8 features.
#
#  Importance is read from the trained model, so no retraining is needed.
#
#  Usage:  python getBDTImportance_plot.py
# ===========================================================================
import os
import re
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Patch
import xgboost as xgb

MODEL_FILE = "_tools/models/2024_2026-06-13_fourmuonSV_charge0/XGB_NOJETS_2024_A_combined.model"
TAG        = "XGB_NOJETS_2024_A_combined"
OUT_DIR    = "_tools/output_getSignalEff_noMET"
TOPN       = 25

# 'gain'       = average gain per split (what the icenet plot shows)
# 'total_gain' = gain summed over splits = the variable's total contribution
METRICS = ["gain", "total_gain"]

# Categorical slots 1-5, fixed order, from the dataviz reference palette
# (pre-validated; not cycled). Colour = collection, i.e. the entity, not rank.
CAT_COL = {
    "Muon":       "#2a78d6",   # blue
    "muonSV":     "#008300",   # green
    "fourmuonSV": "#e87ba4",   # magenta
    "SV":         "#eda100",   # yellow
    "Event":      "#1baf7a",   # aqua
}
INK, INK_MUTED, GRID = "#0b0b0b", "#52514e", "#d8d7d2"


def category(feat):
    """Collection a feature belongs to. Order matters: fourmuonSV_ before
    muonSV_, and Muon_ (capital M) is a different collection from muonSV_."""
    for p, c in (("fourmuonSV_", "fourmuonSV"), ("muonSV_", "muonSV"),
                 ("Muon_", "Muon"), ("SV_", "SV")):
        if feat.startswith(p):
            return c
    return "Event"


def barpage(pdf, items, total, title, subtitle, xlabel):
    """One ranked horizontal bar page. items = [(label, value), ...] desc."""
    fig, ax = plt.subplots(figsize=(8.5, 0.26 * len(items) + 2.2))
    labels = [k for k, _ in items][::-1]          # matplotlib y grows upward
    vals   = [v for _, v in items][::-1]
    cols   = [CAT_COL[category(k)] for k in labels]

    ax.barh(range(len(vals)), vals, color=cols, height=0.62)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=7.5, fontfamily="monospace", color=INK)
    ax.set_xlabel(xlabel, fontsize=9, color=INK)
    ax.tick_params(axis="x", labelsize=8, colors=INK_MUTED)
    ax.tick_params(axis="y", length=0)

    # share-of-total label per bar, in ink not the series colour
    span = max(vals) if vals else 1.0
    for i, v in enumerate(vals):
        ax.text(v + 0.012 * span, i, f"{100*v/total:.1f}%",
                va="center", fontsize=6.5, color=INK_MUTED)
    ax.set_xlim(0, span * 1.12)

    # recessive grid/axes
    ax.xaxis.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(GRID)

    # legend: only the categories actually present, in fixed palette order
    present = [c for c in CAT_COL if any(category(k) == c for k in labels)]
    ax.legend(handles=[Patch(facecolor=CAT_COL[c], label=c) for c in present],
              fontsize=7.5, frameon=False, loc="lower right", ncol=1)

    # subtitle may be multi-line; the title pad has to clear it
    nsub = subtitle.count("\n") + 1
    ax.set_title(title, fontsize=11, color=INK, loc="left", pad=8 + 10 * nsub)
    ax.text(0, 1.005, subtitle, transform=ax.transAxes, fontsize=7.5,
            color=INK_MUTED, va="bottom", linespacing=1.5)
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def main():
    booster = xgb.Booster()
    booster.load_model(MODEL_FILE)
    nfeat = len(booster.feature_names)

    os.makedirs(OUT_DIR, exist_ok=True)
    for metric in METRICS:
        score = booster.get_score(importance_type=metric)
        total = sum(score.values())
        used  = len(score)

        indiv = sorted(score.items(), key=lambda x: -x[1])[:TOPN]
        cum_i = 100 * sum(v for _, v in indiv) / total

        agg = defaultdict(float)
        for k, v in score.items():
            agg[re.sub(r"_\d+$", "", k)] += v
        base = sorted(agg.items(), key=lambda x: -x[1])[:TOPN]
        cum_b = 100 * sum(v for _, v in base) / total

        pretty = {"gain": "gain (average per split)",
                  "total_gain": "total gain (summed over splits)"}[metric]
        head = "Scenario A, 2024  |  Mu10 OR DoubleMu  |  XGB-NOJETS"
        out = f"{OUT_DIR}/bdt_importance_{metric}_{TAG}.pdf"
        with PdfPages(out) as pdf:
            barpage(pdf, indiv, total,
                    f"BDT feature importance — {pretty}",
                    f"{head}\ntop {TOPN} of {used} features split on "
                    f"({nfeat} inputs, {nfeat-used} never used) "
                    f"= {cum_i:.0f}% of total importance",
                    pretty)
            barpage(pdf, base, total,
                    f"BDT importance by variable — {pretty}",
                    f"{head}\nsummed over the 8 leading-object slots; "
                    f"top {TOPN} of {len(agg)} variables "
                    f"= {cum_b:.0f}% of total importance",
                    pretty)
        print(f"[imp] {metric:11s}: {used}/{nfeat} features used, "
              f"top{TOPN} = {cum_i:.0f}% (indiv) / {cum_b:.0f}% (by var) -> {out}")


if __name__ == "__main__":
    main()
