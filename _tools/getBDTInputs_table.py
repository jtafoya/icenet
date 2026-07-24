#!/usr/bin/env python
# ===========================================================================
#  getBDTInputs_table.py  --  compact table of the BDT input variables.
#
#  The variable list is read from the TRAINED MODEL's feature names, not from
#  mvavars_*.py, so the table always reflects what the BDT actually eats after
#  the `exclude_MVA_vars` regexes in models.yml have been applied. Any feature
#  not covered by VARMAP below is reported and dumped into an "unmapped" row --
#  the table can never silently omit an input.
#
#  Writes bdt_inputs_<TAG>.tex and compiles it with pdflatex.
#
#  Usage:  python getBDTInputs_table.py
# ===========================================================================
import os
import re
import subprocess
import sys

import xgboost as xgb

# --- inputs -----------------------------------------------------------------
MODEL_FILE = "_tools/models/2024_2026-06-13_fourmuonSV_charge0/XGB_NOJETS_2024_A_combined.model"
TAG        = "XGB_NOJETS_2024_A_combined"
OUT_DIR    = "_tools/output_getSignalEff_noMET"

# Column order (as requested)
COLS = ["Kinematic", "Position", "Quality/ID"]

# Column widths [cm]: Category, then one per COLS. Narrow on purpose -- cells
# wrap onto several lines rather than the table running the full page width.
W_CAT   = 2.6
W_COL   = {"Kinematic": 3.0, "Position": 2.5, "Quality/ID": 5.4}
TABCOLSEP_CM = 4 * 0.03515                    # 4pt in cm, matches \tabcolsep below
# span of the three variable columns (for the unmapped row) and of the whole
# table (for the footnote), including the inter-column separation
W_SPAN3 = sum(W_COL.values()) + 4 * TABCOLSEP_CM
W_TABLE = W_CAT + W_SPAN3 + 2 * TABCOLSEP_CM

# Category order in the table
CATS = ["Muon", "muonSV", "fourmuonSV", "SV", "Event"]

# Pretty category labels
CAT_TEX = {
    "Muon":       r"Muon",
    "muonSV":     r"muonSV $(2\mu)$",
    "fourmuonSV": r"fourmuonSV $(4\mu)$",
    "SV":         r"SV",
    "Event":      r"Event",
}

# base variable -> (category, column, latex label).
# mu1/mu2/mu3/mu4 variants map to the SAME label on purpose: they are deduped
# per cell, so the four-muon kinematics collapse to one "per muon" entry.
_MU_KIN = {
    "pt":  ("Kinematic", r"$p_{\mathrm{T}}$"),
    "eta": ("Kinematic", r"$\eta$"),
    "phi": ("Kinematic", r"$\phi$"),
}
_VTX = {
    # SV carries its own kinematics; muonSV/fourmuonSV carry per-muon ones only
    # (entries for fields a collection does not have are simply never hit).
    "pt":      ("Kinematic",  r"$p_{\mathrm{T}}$"),
    "eta":     ("Kinematic",  r"$\eta$"),
    "phi":     ("Kinematic",  r"$\phi$"),
    "chi2":    ("Quality/ID", r"$\chi^{2}/\mathrm{dof}$"),
    "pAngle":  ("Quality/ID", r"pAngle"),
    "dlen":    ("Position",   r"$L_{3\mathrm{D}}$"),
    "dlenSig": ("Quality/ID", r"$L_{3\mathrm{D}}/\sigma$"),
    "dxy":     ("Position",   r"$L_{xy}$"),
    "dxySig":  ("Quality/ID", r"$L_{xy}/\sigma$"),
    "x":       ("Position",   r"$x$"),
    "y":       ("Position",   r"$y$"),
    "z":       ("Position",   r"$z$"),
    "mass":    ("Kinematic",  r"$m$"),
    "ndof":    ("Quality/ID", r"$n_{\mathrm{dof}}$"),
}

VARMAP = {
    # --- event multiplicities ---
    "nMuon":       ("Event", "Kinematic", r"$n_{\mu}$"),
    "nSV":         ("Event", "Kinematic", r"$n_{\mathrm{SV}}$"),
    "nmuonSV":     ("Event", "Kinematic", r"$n_{\mu\mathrm{SV}}$"),
    "nfourmuonSV": ("Event", "Kinematic", r"$n_{4\mu\mathrm{SV}}$"),
    # --- muons ---
    "Muon_pt":               ("Muon", "Kinematic",  r"$p_{\mathrm{T}}$"),
    "Muon_eta":              ("Muon", "Kinematic",  r"$\eta$"),
    "Muon_phi":              ("Muon", "Kinematic",  r"$\phi$"),
    "Muon_dxy":              ("Muon", "Position",   r"$d_{xy}$"),
    "Muon_dz":               ("Muon", "Position",   r"$d_{z}$"),
    "Muon_ip3d":             ("Muon", "Position",   r"$\mathrm{IP}_{3\mathrm{D}}$"),
    "Muon_ptErr":            ("Muon", "Quality/ID", r"$\sigma(p_{\mathrm{T}})$"),
    "Muon_dxyErr":           ("Muon", "Quality/ID", r"$\sigma(d_{xy})$"),
    "Muon_dzErr":            ("Muon", "Quality/ID", r"$\sigma(d_{z})$"),
    "Muon_sip3d":            ("Muon", "Quality/ID", r"$\mathrm{IP}_{3\mathrm{D}}/\sigma$"),
    "Muon_charge":           ("Muon", "Quality/ID", r"$q$"),
    "Muon_tightId":          ("Muon", "Quality/ID", r"tightId"),
    "Muon_softMva":          ("Muon", "Quality/ID", r"softMVA"),
    "Muon_pfRelIso03_all":   ("Muon", "Quality/ID", r"$I_{\mathrm{PF}}^{\Delta R<0.3}$"),
    "Muon_miniPFRelIso_all": ("Muon", "Quality/ID", r"$I_{\mathrm{mini}}$"),
    "Muon_jetIdx":           ("Muon", "Quality/ID", r"jet idx"),
}
# vertex collections: <coll>_<field> and <coll>_mu<N><kin>
for _coll in ("muonSV", "fourmuonSV", "SV"):
    for _f, (_c, _l) in _VTX.items():
        VARMAP[f"{_coll}_{_f}"] = (_coll, _c, _l)
    for _n in (1, 2, 3, 4):
        for _k, (_c, _l) in _MU_KIN.items():
            VARMAP[f"{_coll}_mu{_n}{_k}"] = (_coll, _c, _l)


def latex_escape(s):
    return s.replace("_", r"\_")


def main():
    booster = xgb.Booster()
    booster.load_model(MODEL_FILE)
    feats = booster.feature_names
    if not feats:
        sys.exit(f"[inputs] ERROR: no feature names stored in {MODEL_FILE}")

    # collapse the _0.._7 leading-object index suffix
    base, mult = [], {}
    for f in feats:
        b = re.sub(r"_\d+$", "", f)
        if b not in base:
            base.append(b)
            mult[b] = 0
        mult[b] += 1

    # cells[cat][col] = list of labels (deduped, order preserved)
    cells = {c: {k: [] for k in COLS} for c in CATS}
    unmapped = []
    for b in base:
        if b not in VARMAP:
            unmapped.append(b)
            continue
        cat, col, lab = VARMAP[b]
        if lab not in cells[cat][col]:
            cells[cat][col].append(lab)

    if unmapped:
        print(f"[inputs] WARNING: {len(unmapped)} feature(s) not in VARMAP: {unmapped}")

    print(f"[inputs] {len(feats)} BDT features, {len(base)} distinct variables")
    for c in CATS:
        n = sum(len(cells[c][k]) for k in COLS)
        print(f"[inputs]   {c:12s} {n:2d} entries")

    # ---- LaTeX (preamble matches the other tables in this directory) --------
    os.makedirs(OUT_DIR, exist_ok=True)
    tex = f"{OUT_DIR}/bdt_inputs_{TAG}.tex"
    with open(tex, "w") as f:
        f.write(r"""\documentclass[10pt,a4paper]{article}
\usepackage[a4paper,top=1.5cm,bottom=1.5cm,left=1.5cm,right=1.5cm]{geometry}
\usepackage{booktabs}
\usepackage{amsmath,amssymb}
\usepackage[T1]{fontenc}
\usepackage{helvet}
\renewcommand{\familydefault}{\sfdefault}
\setlength{\parindent}{0pt}
\pagestyle{empty}
\setlength{\tabcolsep}{4pt}
\renewcommand{\arraystretch}{1.25}
\begin{document}

{\large BDT input variables}

\vspace{0.3em}
{\footnotesize Scenario A, 2024 $|$ \texttt{Mu10 OR DoubleMu} $|$ model \texttt{XGB-NOJETS}}

\vspace{0.8em}

\footnotesize
""" + (r"\begin{tabular}{@{}p{%.2fcm} p{%.2fcm} p{%.2fcm} p{%.2fcm}@{}}"
       % (W_CAT, W_COL["Kinematic"], W_COL["Position"], W_COL["Quality/ID"])) + r"""
\toprule
\textbf{Category} & \textbf{Kinematic} & \textbf{Position} & \textbf{Quality/ID} \\
\midrule
""")
        for i, c in enumerate(CATS):
            row = [CAT_TEX[c]]
            for k in COLS:
                row.append(", ".join(cells[c][k]) if cells[c][k] else r"--")
            # thin full-width rule between rows (\bottomrule closes the last one)
            if i:
                f.write(r"\midrule[\cmidrulewidth]" + "\n")
            f.write(" & ".join(row) + r" \\" + "\n")
        if unmapped:
            f.write(r"\midrule" + "\n")
            f.write(r"\textit{unmapped} & \multicolumn{3}{@{}p{%.2fcm}@{}}{" % W_SPAN3
                    + ", ".join(r"\texttt{%s}" % latex_escape(u) for u in unmapped)
                    + r"} \\" + "\n")
        f.write(r"""\bottomrule
\end{tabular}

\vspace{0.8em}
""" + (r"\begin{minipage}{%.2fcm}" % W_TABLE) + r"""
{\footnotesize
Jagged collections are stored for the 8 leading objects (suffix \texttt{\_0}\,--\,\texttt{\_7}):
""" + f"{len(base)} distinct variables give {len(feats)} BDT inputs." + r"""
Per-muon kinematics of a vertex ($\mu_{1..4}$) are listed once.
Variables removed by \texttt{exclude\_MVA\_vars}: all \texttt{Jet\_*} (NOJETS),
the SV masses \texttt{muonSV\_mass}/\texttt{SV\_mass} (mass-sculpting),
\texttt{PFMET\_*}, the $\Delta R$ variables and the SV charges.
}
\end{minipage}

\end{document}
""")

    print(f"[inputs] -> {tex}")
    r = subprocess.run(["pdflatex", "-interaction=nonstopmode",
                        "-output-directory", OUT_DIR, tex],
                       capture_output=True, text=True)
    pdf = f"{OUT_DIR}/bdt_inputs_{TAG}.pdf"
    if not os.path.exists(pdf):
        print(r.stdout[-2500:])
        sys.exit("[inputs] ERROR: pdflatex failed")
    for ext in ("aux", "log"):
        p = f"{OUT_DIR}/bdt_inputs_{TAG}.{ext}"
        if os.path.exists(p):
            os.remove(p)
    print(f"[inputs] -> {pdf}")


if __name__ == "__main__":
    main()
