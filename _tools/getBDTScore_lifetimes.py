#!/usr/bin/env python
# ===========================================================================
#  getBDTScore_lifetimes.py  --  BDT score, four lifetimes of a signal point
#  overlaid on one panel (+ QCD reference), for all 6 latest BDTs.
#
#  One multi-page PDF per BDT (scenario x trigger); one page per signal mass
#  point, each page overlaying the four lifetimes (ctau = 0.1/1/10/100 mm) of
#  that point plus the full-QCD reference. Score = the BDT's own XGB-NOJETS
#  output (the deployed model).
#
#  QCD reference uses ALL background: icenet assigns GEN_ctau to background at
#  random, and QCD has no generator lifetime anyway.
#
#  Signal points:
#    scenarioA : (4,0.40), (4,1.33), (10,1.00), (10,3.33)   [as requested]
#    scenarioB1: (1,0.33), (2,0.40), (2,0.67), (4,1.33)      [B1's own grid --
#      the requested A-points (4/0.40, 10/1.00, 10/3.33) do NOT exist for B1,
#      whose grid tops out at mpi=4; only 4/1.33 overlaps]
#
#  Usage:  python getBDTScore_lifetimes.py
# ===========================================================================
import os
import pickle

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import AutoMinorLocator
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

BASE = "figs/dqcd"
# run label -> (scenario, trigger, eval-pickle relative path)
RUNS = [
    ("A_combined",  "A",  "combined",
     # "config__tune0_2024_new.yml/inputmap__mc_map__scenarioA_all_2024.yml--modeltag__scenarioA_all_no_DA_old_BDT_with_dRmSV_2024/2026-07-24_10-22-33_lxfw00"),
     "config__tune0_2024_new.yml/inputmap__mc_map__scenarioA_all_2024.yml--modeltag__scenarioA_all_no_DA_old_BDT_2024/2026-07-24_10-22-33_lxfw00"),
    ("A_Mu10",      "A",  "Mu10",
     # "config__tune0_2024_Mu10_new.yml/inputmap__mc_map__scenarioA_all_2024.yml--modeltag__scenarioA_all_no_DA_old_BDT_with_dRmSV_2024_Mu10/2026-07-23_21-20-15_lxbgpu01"),
     "config__tune0_2024_Mu10_new.yml/inputmap__mc_map__scenarioA_all_2024.yml--modeltag__scenarioA_all_no_DA_old_BDT_2024_Mu10/2026-07-23_21-20-15_lxbgpu01"),
    ("A_DoubleMu",  "A",  "DoubleMu",
     # "config__tune0_2024_DoubleMu_new.yml/inputmap__mc_map__scenarioA_all_2024.yml--modeltag__scenarioA_all_no_DA_old_BDT_with_dRmSV_2024_DoubleMu/2026-07-23_21-49-14_lxbgpu00"),
     "config__tune0_2024_DoubleMu_new.yml/inputmap__mc_map__scenarioA_all_2024.yml--modeltag__scenarioA_all_no_DA_old_BDT_2024_DoubleMu/2026-07-23_21-49-14_lxbgpu00"),
    ("B1_combined", "B1", "combined",
     # "config__tune0_2024_new.yml/inputmap__mc_map__scenarioB1_all_2024.yml--modeltag__scenarioB1_all_no_DA_old_BDT_with_dRmSV_2024/2026-07-24_08-15-48_lxfw00"),
     "config__tune0_2024_new.yml/inputmap__mc_map__scenarioB1_all_2024.yml--modeltag__scenarioB1_all_no_DA_old_BDT_2024/2026-07-24_08-15-48_lxfw00"),
    ("B1_Mu10",     "B1", "Mu10",
     # "config__tune0_2024_Mu10_new.yml/inputmap__mc_map__scenarioB1_all_2024.yml--modeltag__scenarioB1_all_no_DA_old_BDT_with_dRmSV_2024_Mu10/2026-07-23_21-33-35_lxcgpu03"),
     "config__tune0_2024_Mu10_new.yml/inputmap__mc_map__scenarioB1_all_2024.yml--modeltag__scenarioB1_all_no_DA_old_BDT_2024_Mu10/2026-07-23_21-33-35_lxcgpu03"),
    ("B1_DoubleMu", "B1", "DoubleMu",
     # "config__tune0_2024_DoubleMu_new.yml/inputmap__mc_map__scenarioB1_all_2024.yml--modeltag__scenarioB1_all_no_DA_old_BDT_with_dRmSV_2024_DoubleMu/2026-07-23_22-26-45_lxbgpu01"),
     "config__tune0_2024_DoubleMu_new.yml/inputmap__mc_map__scenarioB1_all_2024.yml--modeltag__scenarioB1_all_no_DA_old_BDT_2024_DoubleMu/2026-07-23_22-26-45_lxbgpu01"),
]

# signal (mpi, mA) points per scenario
POINTS = {
    "A":  [(4, 0.40), (4, 1.33), (10, 1.00), (10, 3.33)],
    "B1": [(1, 0.33), (2, 0.40), (2, 0.67), (4, 1.33)],
}
SCEN_TEX = {"A": "Scenario A", "B1": "Scenario B1"}

PRED_IDX = 1                 # 0 = XGB, 1 = XGB-NOJETS (deployed)
OUT_DIR  = "_tools/output_getSignalEff_noMET"

CTAUS   = [0.1, 1.0, 10.0, 100.0]
CCOL    = ["#d62a2a", "#eb6834", "#008300", "#2a78d6"]   # red, orange, green, blue
QCD_COL = "#8c8c8c"
NBINS   = 50
LUMI_TEX = r"($13.6\,\mathrm{TeV}$, 2024)"               # energy/year in parens
INK, INK_MUTED, GRID = "#0b0b0b", "#52514e", "#d8d7d2"


def make_page(score, y, w, mpi, mA, ctau, scen, trig, mpi_t, mA_t):
    """One page: the four lifetimes of (mpi_t, mA_t) overlaid, plus full QCD."""
    bins = np.linspace(0., 1., NBINS + 1)
    fig, ax = plt.subplots(figsize=(7.6, 5.2))

    def norm(mask):
        ww = w[mask]
        return ww / ww.sum() if ww.sum() > 0 else ww

    b = (y == 0)
    ax.hist(score[b], bins=bins, weights=norm(b), histtype="stepfilled",
            color=QCD_COL, alpha=0.30)
    ax.hist(score[b], bins=bins, weights=norm(b), histtype="step",
            color=QCD_COL, lw=1.4)

    handles, labels = [], []
    for ct, col in zip(CTAUS, CCOL):
        m = (y == 1) & np.isclose(mpi, mpi_t) & np.isclose(mA, mA_t) & np.isclose(ctau, ct)
        n = int(m.sum())
        if n > 0:
            ax.hist(score[m], bins=bins, weights=norm(m), histtype="step",
                    color=col, lw=1.8)
        handles.append(Line2D([], [], color=col, lw=1.8))
        lab = (rf"{SCEN_TEX[scen]}, $m_{{\bar{{\pi}}}} = {mpi_t:g}$, "
               rf"$m_{{A'}} = {mA_t:g}$, $c\tau = {ct:g}\,\mathrm{{mm}}$")
        if n == 0:
            lab += " (no events)"
        labels.append(lab)
        print(f"[score] {scen}_{trig}  mpi={mpi_t:g} mA={mA_t:g} ctau={ct:6g}: {n} events")

    ax.set_xlabel("BDT score", fontsize=10, color=INK, loc="right")
    ax.set_ylabel(r"Normalised entries (a.u.)", fontsize=10, color=INK, loc="top")
    ax.set_xlim(0., 1.)
    for s in ("top", "right", "bottom", "left"):
        ax.spines[s].set_visible(True)
        ax.spines[s].set_color(INK)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which="both", direction="in", top=True, right=True,
                   left=True, bottom=True, colors=INK)
    ax.tick_params(which="major", length=6, labelsize=9)
    ax.tick_params(which="minor", length=3)

    handles.append(Patch(facecolor=QCD_COL, alpha=0.30, edgecolor=QCD_COL))
    labels.append("QCD")
    leg = ax.legend(handles, labels, fontsize=8.5, frameon=False, loc="upper left",
                    handlelength=1.6, title=f"{SCEN_TEX[scen]} BDT ({trig})")
    leg.get_title().set_fontsize(9.5)
    leg._legend_box.align = "left"

    ax.text(0, 1.01, "Preliminary", transform=ax.transAxes, fontsize=11,
            color=INK, style="italic", va="bottom", ha="left")
    ax.text(1, 1.01, LUMI_TEX, transform=ax.transAxes, fontsize=11,
            color=INK, va="bottom", ha="right")
    fig.tight_layout()
    return fig


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    for run, scen, trig, rel in RUNS:
        path = f"{BASE}/{rel}/eval/eval_results.pkl"
        print(f"\n=== {run} : {path}")
        with open(path, "rb") as f:
            d = pickle.load(f)
        dat = d["data"]
        X, y, w = dat.x, dat.y, dat.w
        ids = list(dat.ids)
        score = np.asarray(d["y_preds"][PRED_IDX], dtype=float)
        mpi  = X[:, ids.index("GEN_mpi")].astype(float)
        mA   = X[:, ids.index("GEN_mA")].astype(float)
        ctau = X[:, ids.index("GEN_ctau")].astype(float)

        out = f"{OUT_DIR}/bdt_score_signalpoints_{run}_XGB-NOJETS.pdf"
        with PdfPages(out) as pdf:
            for mpi_t, mA_t in POINTS[scen]:
                fig = make_page(score, y, w, mpi, mA, ctau, scen, trig, mpi_t, mA_t)
                pdf.savefig(fig)
                plt.close(fig)
        print(f"[score] -> {out} ({len(POINTS[scen])} pages)")


if __name__ == "__main__":
    main()
