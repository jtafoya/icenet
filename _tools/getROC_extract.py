#!/usr/bin/env python
# ===========================================================================
#  getROC_extract.py  --  per-lxy-bin ROC curves from the BDT eval output.
#
#  Reads the icenet eval pickle (test set: per-event BDT score, class label,
#  weight, and the input features) and writes two CSVs consumed by
#  getROC_render.C:
#
#    roc_curve_<MODEL>.csv   lxy_bin, fpr, tpr, thr, n_bkg_raw
#    roc_points_<MODEL>.csv  lxy_bin, ctau, fpr, tpr, thr, n_bkg_raw, n_eff
#
#  Scenario A, SingleORDouble trigger (filter_standard_2024 = OR of
#  HLT_Mu10_Barrel_L1HP11_IP6 and HLT_DoubleMu4_3_LowMass), model XGB-NOJETS
#  (= the xgb01_NOJETS branch used by the deployment study).
#
#  Weights: icenet balances the two classes (sum w_sig == sum w_bkg), but
#  WITHIN a class the weights are physical -- each QCD bin carries xs*eff_acc,
#  each signal sample is normalised to 1. TPR and FPR are both within-class
#  ratios, so the balancing cancels: the FPR is xs-weighted QCD efficiency and
#  the pooled signal is the equal-per-sample mix the BDT was trained on.
#
#  Usage:  python getROC_extract.py
# ===========================================================================
import os
import pickle

import numpy as np

# --- inputs -----------------------------------------------------------------
# Scenario A / SingleORDouble training, same modeltag as the deployed model
# (_tools/models/2024_2026-06-13_fourmuonSV_charge0).
EVAL_PKL = (
    "figs/dqcd/config__tune0_2024_new.yml/"
    "inputmap__mc_map__scenarioA_all_2024.yml"
    # "--modeltag__scenarioA_all_no_DA_old_BDT_with_dRmSV_2024/"
    "--modeltag__scenarioA_all_no_DA_old_BDT_2024/"
    "2026-07-24_10-22-33_lxfw00/eval/eval_results.pkl"
)
MODEL     = "Mu10ORDoubleMu"
PRED_IDX  = 1            # y_preds index: 0 = XGB, 1 = XGB-NOJETS
OUT_DIR   = "_tools/output_getSignalEff_noMET"

# lxy (= muonSV_dxy_0) category boundaries [cm]; same as getSignalEff_render.C.
# A 4th "inclusive" bin (any reconstructed muonSV, dxy>=0) is appended so the
# overlay plot can show the three displacement bins against the whole sample.
LXY_LO  = [0.0, 1.0, 10.0]
LXY_HI  = [1.0, 10.0, 100.0]
LXY_TAG = ["lxy[0,1]cm", "lxy[1,10]cm", "lxy[10,100]cm", "inclusive"]

# marker mass point + its four lifetimes
MPI, MA  = 10.0, 1.00
CTAUS    = [0.1, 1.0, 10.0, 100.0]
FPR_MARK = 1e-4          # analysis working point

# ROC sampling grid (log-spaced in FPR)
FPR_GRID = np.logspace(-5, 0, 500)


def roc_from(scores, weights, grid):
    """Weighted reverse-cumulative background efficiency at each grid FPR.

    Scans the score from high to low and returns, per grid point,
      thr   : the cut yielding that background efficiency
      eff   : the efficiency actually attained (>= the request; a single
              high-xs event can carry more weight than the requested FPR,
              in which case the grid point is unreachable and eff jumps)
      n_raw : raw MC events above thr
      n_eff : effective events above thr, (sum w)^2 / sum w^2. The QCD weights
              span ~7 orders of magnitude (xs from 2.8e6 to 1.1 pb), so n_raw
              badly overstates the statistical power of the tail -- n_eff is
              what decides whether a point on the curve means anything.
    """
    order = np.argsort(-scores)
    s, w = scores[order], weights[order]
    cw, cw2 = np.cumsum(w), np.cumsum(w ** 2)
    cum = cw / w.sum()
    idx = np.clip(np.searchsorted(cum, grid), 0, len(s) - 1)
    n_eff = cw[idx] ** 2 / cw2[idx]
    return s[idx], cum[idx], idx + 1, n_eff


def main():
    print(f"[roc] reading {EVAL_PKL}")
    with open(EVAL_PKL, "rb") as f:
        d = pickle.load(f)

    dat = d["data"]
    X, y, w = dat.x, dat.y, dat.w
    ids = list(dat.ids)
    score = np.asarray(d["y_preds"][PRED_IDX], dtype=float)

    dxy  = X[:, ids.index("muonSV_dxy_0")].astype(float)
    mpi  = X[:, ids.index("GEN_mpi")].astype(float)
    mA   = X[:, ids.index("GEN_mA")].astype(float)
    ctau = X[:, ids.index("GEN_ctau")].astype(float)

    # events without a reconstructed muonSV carry the -999 sentinel; they fall
    # in no lxy bin and are dropped by the bin masks below.
    n_nosv = int(np.sum(dxy < 0))
    print(f"[roc] {len(y)} events, {n_nosv} without a muonSV (dxy<0, dropped)")

    os.makedirs(OUT_DIR, exist_ok=True)
    f_curve = open(f"{OUT_DIR}/roc_curve_{MODEL}.csv", "w")
    f_point = open(f"{OUT_DIR}/roc_points_{MODEL}.csv", "w")
    # lxy tags contain a comma -> always quoted, as in the other CSVs here
    f_curve.write("lxy_bin,fpr,tpr,thr,n_bkg_raw,n_bkg_eff\n")
    f_point.write("lxy_bin,ctau_mm,fpr,tpr,thr,n_bkg_raw,n_bkg_eff,n_sig_raw\n")

    for c in range(len(LXY_TAG)):
        if LXY_TAG[c] == "inclusive":
            inb = dxy >= 0.0                       # any reconstructed muonSV
        else:
            inb = (dxy >= LXY_LO[c]) & (dxy < LXY_HI[c])
        b, s = inb & (y == 0), inb & (y == 1)

        # ROC curve: pooled scenario A signal vs xs-weighted QCD, in this bin
        thr, fpr, nraw, neff = roc_from(score[b], w[b], FPR_GRID)
        ws, ss = w[s], score[s]
        tpr = np.array([ws[ss >= t].sum() / ws.sum() for t in thr])
        for i in range(len(FPR_GRID)):
            f_curve.write(f'"{LXY_TAG[c]}",{fpr[i]:.6e},{tpr[i]:.6e},'
                          f"{thr[i]:.6f},{nraw[i]:d},{neff[i]:.3f}\n")

        # working-point threshold in this bin, and the four lifetime markers
        t_wp, f_wp, nraw_wp, neff_wp = roc_from(score[b], w[b],
                                                np.array([FPR_MARK]))
        t_wp, f_wp = float(t_wp[0]), float(f_wp[0])

        for ct in CTAUS:
            # GEN_* are float32: exact == against 0.1 never matches
            m = (inb & (y == 1) & np.isclose(mpi, MPI) & np.isclose(mA, MA)
                 & np.isclose(ctau, ct))
            wm, sm = w[m], score[m]
            tpr_m = wm[sm >= t_wp].sum() / wm.sum() if wm.sum() > 0 else 0.0
            f_point.write(f'"{LXY_TAG[c]}",{ct:g},{f_wp:.6e},{tpr_m:.6e},'
                          f"{t_wp:.6f},{int(nraw_wp[0]):d},{neff_wp[0]:.2f},"
                          f"{int(m.sum()):d}\n")

        flag = "" if abs(f_wp - FPR_MARK) < 0.5 * FPR_MARK else "  <-- 1e-4 UNREACHABLE"
        print(f"[roc] {LXY_TAG[c]:14s} bkg_raw={int(b.sum()):7d} "
              f"sig_raw={int(s.sum()):6d} | thr={t_wp:.5f} FPR={f_wp:.2e} "
              f"raw_above={int(nraw_wp[0]):3d} N_eff={neff_wp[0]:.1f}{flag}")

    f_curve.close()
    f_point.close()
    print(f"[roc] -> {OUT_DIR}/roc_curve_{MODEL}.csv")
    print(f"[roc] -> {OUT_DIR}/roc_points_{MODEL}.csv")


if __name__ == "__main__":
    main()
