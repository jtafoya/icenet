// ===========================================================================
//  getSignalEff_compare.C  --  overlay THIS analysis' Asimov-significance (Z)
//  results against Javier's parallel scouting study.
//
//  * THIS analysis  -> solid lines + solid (filled) markers.
//    Read from the flat CSV produced by getSignalEff_render.C:
//        output_getSignalEff_noMET/[test_]full_by_lxy_Mu10ORDoubleMu.csv
//    (columns: signal_point, mA_GeV, ctau_mm, lxy_bin, FPR, ..., Z, ...)
//
//  * Javier (scouting) -> dashed lines + open (hollow) markers.
//    Hardcoded below from "Z values Scouting.docx" (4 mass points x 4 lxy
//    categories x 4 FPR rows x 4 ctau columns).
//
//  Colour  = FPR working point (10^-1..10^-4).
//  Style   = author (solid/filled = us, dashed/open = Javier).
//
//  Layout: one PDF, 4 pages (one per mass point), each a 2x2 grid of the four
//  lxy categories; per-page a title band + a shared horizontal legend.
//
//  Usage:
//    root -l -b -q 'getSignalEff_compare.C("test")'   # read test_ CSV (default off)
//    root -l -b -q 'getSignalEff_compare.C()'         # read production CSV
// ===========================================================================
#include "TCanvas.h"
#include "TPad.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TColor.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>

// ---- axis grids ------------------------------------------------------------
static const int NM = 4;   // mass points
static const int NL = 4;   // lxy categories
static const int NF = 4;   // FPR working points
static const int NC = 4;   // ctau columns

static const double CTAU[NC] = { 0.1, 1.0, 10.0, 100.0 };      // mm
static const double FPRV[NF] = { 1e-1, 1e-2, 1e-3, 1e-4 };

// mass points (order matches Javier's docx)
struct MassPt { int mpi; double mA; const char* tex; };
static const MassPt MASS[NM] = {
    { 4,  0.40, "m_{#bar{#pi}}=4 GeV, m_{A'}=0.40 GeV" },
    { 1,  0.33, "m_{#bar{#pi}}=1 GeV, m_{A'}=0.33 GeV" },
    { 4,  1.33, "m_{#bar{#pi}}=4 GeV, m_{A'}=1.33 GeV" },
    { 10, 1.00, "m_{#bar{#pi}}=10 GeV, m_{A'}=1.00 GeV" },
};

// Dimuon mass window used for the Z values in the CSV: MWIN_FRAC_LO*m_A to
// MWIN_FRAC_HI*m_A. Mirrored from getSignalEff_render.C (USE_MASS_WINDOW,
// MWIN_FRAC_LO/HI there), which is the source of truth -- this script only
// reads the CSV, so the values must be kept in sync by hand.
static const double MWIN_FRAC_LO = 0.9;
static const double MWIN_FRAC_HI = 1.1;

// lxy categories: our CSV tag vs pretty label (order matches Javier's docx)
static const char* LXY_CSV[NL] = { "lxy[0,1]cm", "lxy[1,10]cm", "lxy[10,100]cm", "inclusive" };
static const char* LXY_TEX[NL] = { "l_{xy}#in[0,1] cm", "l_{xy}#in[1,10] cm",
                                   "l_{xy}#in[10,100] cm", "all l_{xy}" };

// ---- Javier's Z values [mass][lxy][fpr][ctau] (from the docx) --------------
static const double JAV[NM][NL][NF][NC] = {
  // ---- mpi=4, mA=0.40 ----
  {
    {{1.9,0.6,0.1,0.0},{3.9,1.4,0.2,0.0},{6.0,2.6,0.3,0.0},{8.2,4.2,0.5,0.0}},   // 0-1
    {{1.3,4.4,1.1,0.1},{2.4,8.0,2.0,0.2},{3.7,13.7,3.4,0.4},{4.8,20.1,4.9,0.6}}, // 1-10
    {{0.0,0.5,0.9,0.2},{0.0,0.9,1.7,0.3},{0.0,1.7,3.2,0.5},{0.0,13.0,22.6,4.8}}, // 10-100
    {{2.2,1.6,0.4,0.1},{4.4,4.0,1.1,0.1},{6.8,9.0,2.6,0.3},{9.2,18.3,6.2,0.9}},  // all
  },
  // ---- mpi=1, mA=0.33 ----
  {
    {{0.9,0.3,0.0,0.0},{1.5,0.6,0.1,0.0},{1.9,1.1,0.2,0.0},{1.8,1.4,0.2,0.0}},   // 0-1
    {{0.6,2.1,0.5,0.1},{1.1,3.8,0.9,0.1},{1.6,6.5,1.6,0.2},{0.9,6.0,1.7,0.2}},   // 1-10
    {{0.0,0.2,0.4,0.1},{0.0,0.4,0.8,0.1},{0.0,1.5,2.7,0.5},{0.0,2.0,3.9,0.7}},   // 10-100
    {{1.1,0.8,0.2,0.0},{1.8,1.9,0.5,0.1},{2.4,4.3,1.2,0.2},{1.8,5.9,2.2,0.3}},   // all
  },
  // ---- mpi=4, mA=1.33 ----
  {
    {{0.9,0.7,0.1,0.0},{0.8,1.3,0.2,0.0},{0.7,1.9,0.4,0.0},{0.7,1.7,0.3,0.0}},   // 0-1
    {{0.0,2.2,1.9,0.3},{0.1,4.0,3.4,0.5},{0.1,6.6,6.3,0.9},{0.1,7.4,7.5,1.0}},   // 1-10
    {{0.0,0.0,0.8,0.4},{0.0,0.1,2.0,0.9},{0.0,0.1,5.6,2.5},{0.0,0.3,11.9,5.2}},  // 10-100
    {{0.8,1.2,0.7,0.1},{0.7,2.8,2.0,0.4},{0.6,5.4,5.4,1.2},{0.5,6.2,9.3,2.4}},   // all
  },
  // ---- mpi=10, mA=1.00 ----
  {
    {{1.4,0.6,0.1,0.0},{3.0,1.7,0.2,0.0},{5.1,3.4,0.4,0.0},{10.0,6.9,0.8,0.1}},  // 0-1
    {{0.3,3.6,1.5,0.2},{0.6,7.1,3.0,0.4},{1.2,14.2,6.1,0.7},{2.8,30.5,12.6,1.5}},// 1-10
    {{0.0,0.1,1.1,0.3},{0.0,0.2,2.0,0.5},{0.0,0.6,4.9,1.3},{0.0,1.9,14.5,4.0}},  // 10-100
    {{1.4,1.4,0.6,0.1},{2.9,4.3,1.9,0.3},{4.4,10.8,5.5,0.9},{7.0,26.0,17.0,3.3}},// all
  },
};

// ---- small helpers ---------------------------------------------------------
static std::vector<std::string> split_csv(const std::string& s)
{
    std::vector<std::string> out; std::string cur; bool q = false;
    for (char ch : s) {
        if (ch == '"')            { q = !q; continue; }
        else if (ch == ',' && !q) { out.push_back(cur); cur.clear(); }
        else                       cur += ch;
    }
    out.push_back(cur);
    return out;
}

static int idx_ctau(double t) {
    for (int c = 0; c < NC; ++c) if (std::fabs(t - CTAU[c]) < 1e-6 * std::max(1.0, CTAU[c])) return c;
    // tolerant match
    int best = 0; double bd = 1e30;
    for (int c = 0; c < NC; ++c) { double d = std::fabs(t - CTAU[c]); if (d < bd) { bd = d; best = c; } }
    return best;
}
static int idx_fpr(double f) {
    int k = (int)std::lround(-std::log10(f)) - 1;   // 1e-1 -> 0 ... 1e-4 -> 3
    return (k >= 0 && k < NF) ? k : -1;
}

// ===========================================================================
void getSignalEff_compare(const char* opt = "")
{
    TString O(opt); O.ToLower();
    const bool RUN_TEST = O.Contains("test");
    const char* MODEL   = "Mu10ORDoubleMu";

    gStyle->SetOptStat(0);
    gStyle->SetTextFont(42);

    TString script_dir = gSystem->DirName(__FILE__);
    TString out_dir    = script_dir + "/output_getSignalEff_noMET";
    TString pfx        = RUN_TEST ? "test_" : "";
    TString csv        = out_dir + "/" + pfx + "full_by_lxy_" + MODEL + ".csv";

    // ---- load OUR Z values from the flat CSV -------------------------------
    // MINE[m][l][f][c]; NaN = missing
    double MINE[NM][NL][NF][NC];
    for (int m=0;m<NM;++m) for (int l=0;l<NL;++l) for (int f=0;f<NF;++f) for (int c=0;c<NC;++c)
        MINE[m][l][f][c] = std::nan("");

    std::ifstream in(csv.Data());
    if (!in) { printf("[compare] ERROR: cannot open %s\n", csv.Data()); return; }

    std::string line;
    std::getline(in, line);                       // header
    auto hdr = split_csv(line);
    int i_sp=-1, i_mA=-1, i_ct=-1, i_lxy=-1, i_fpr=-1, i_Z=-1;
    for (size_t i=0;i<hdr.size();++i) {
        std::string h = hdr[i];
        if (h=="signal_point") i_sp=i;
        else if (h=="mA_GeV")  i_mA=i;
        else if (h=="ctau_mm") i_ct=i;
        else if (h=="lxy_bin") i_lxy=i;
        else if (h=="FPR")     i_fpr=i;
        else if (h=="Z")       i_Z=i;
    }
    if (i_sp<0||i_mA<0||i_ct<0||i_lxy<0||i_fpr<0||i_Z<0) {
        printf("[compare] ERROR: unexpected CSV header in %s\n", csv.Data()); return;
    }

    int n_rows=0, n_used=0;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        ++n_rows;
        auto v = split_csv(line);
        if ((int)v.size() <= i_Z) continue;

        // mpi from signal_point ("mpi=4, mA=1.33, ctau=0.1mm")
        int mpi = -1;
        { const std::string& sp = v[i_sp]; size_t p = sp.find("mpi=");
          if (p != std::string::npos) mpi = std::atoi(sp.c_str()+p+4); }
        double mA = std::atof(v[i_mA].c_str());
        double ct = std::atof(v[i_ct].c_str());
        double fp = std::atof(v[i_fpr].c_str());
        double Z  = std::atof(v[i_Z].c_str());

        // match mass point
        int m = -1;
        for (int mm=0; mm<NM; ++mm)
            if (MASS[mm].mpi==mpi && std::fabs(MASS[mm].mA-mA) < 0.05) { m=mm; break; }
        if (m<0) continue;
        // match lxy
        int l = -1;
        for (int ll=0; ll<NL; ++ll) if (v[i_lxy] == LXY_CSV[ll]) { l=ll; break; }
        if (l<0) continue;
        int f = idx_fpr(fp); if (f<0) continue;
        int c = idx_ctau(ct);

        if (Z >= 0.0) { MINE[m][l][f][c] = Z; ++n_used; }
    }
    in.close();
    printf("[compare] read %s : %d data rows, %d matched cells\n", csv.Data(), n_rows, n_used);

    // ---- styling (matches getSignalEff_render.C significance plots) --------
    const int COL[NF] = { kBlue+1, kRed+1, kGreen+2, kMagenta+1 };
    const int MKF[NF] = { 20, 21, 22, 33 };   // filled  markers (this analysis)
    const int MKO[NF] = { 24, 25, 26, 27 };   // hollow  markers (Javier)
    const double LUMI_FB = 109.95;

    // B(A'->mumu) per mass point (MASS order = set1,set3,set0,set2 -> 0.440,0.464,0.317,0.307)
    const double BR_MASS[NM] = { 0.440, 0.464, 0.317, 0.307 };
    const int L_ALL = 3;   // "all lxy" index in the lxy arrays

    std::vector<TObject*> keep;   // keep drawn objects alive until Print

    // helper: overlay the 8 curves (4 FPR x {mine solid/filled, Javier dashed/open})
    // for mass m, lxy l, on the current pad. lw = line width, ms = marker size.
    // Returns per-FPR "mine" graphs (for the FPR legend). ymax updated by ref.
    auto draw_overlay = [&](int m, int l, int lw, double ms,
                            std::vector<TGraph*>& mine_g, double& ymax) {
        for (int f=0; f<NF; ++f) {
            TGraph* gm = new TGraph(); keep.push_back(gm);
            for (int cc=0; cc<NC; ++cc)
                if (!std::isnan(MINE[m][l][f][cc])) {
                    gm->SetPoint(gm->GetN(), CTAU[cc], MINE[m][l][f][cc]);
                    ymax = std::max(ymax, MINE[m][l][f][cc]);
                }
            gm->SetLineColor(COL[f]); gm->SetMarkerColor(COL[f]);
            gm->SetLineWidth(lw); gm->SetLineStyle(1);
            gm->SetMarkerStyle(MKF[f]); gm->SetMarkerSize(ms);
            mine_g.push_back(gm);

            TGraph* gj = new TGraph(); keep.push_back(gj);
            for (int cc=0; cc<NC; ++cc) {
                gj->SetPoint(cc, CTAU[cc], JAV[m][l][f][cc]);
                ymax = std::max(ymax, JAV[m][l][f][cc]);
            }
            gj->SetLineColor(COL[f]); gj->SetMarkerColor(COL[f]);
            gj->SetLineWidth(lw); gj->SetLineStyle(2);
            gj->SetMarkerStyle(MKO[f]); gj->SetMarkerSize(ms);
        }
    };
    // black dummy graphs for the "author" style key (solid=us, dashed=Javier)
    auto make_key = [&](int style, int marker) {
        TGraph* g = new TGraph(1); keep.push_back(g);
        g->SetLineColor(kBlack); g->SetMarkerColor(kBlack);
        g->SetLineWidth(2); g->SetLineStyle(style); g->SetMarkerStyle(marker);
        return g;
    };

    // Draw ONE full-format page (single plot) for mass m, lxy category l, on the
    // current pad/canvas. lxy_label != nullptr adds an "l_xy" line to the info
    // legend (used for the partitioned pages). This is the format the user liked
    // for "all lxy"; the partitioned PDF reuses it verbatim (one plot per page).
    auto draw_full_page = [&](int m, int l, const char* lxy_label) {
        gPad->SetLogx();
        gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.12); gPad->SetRightMargin(0.04);

        std::vector<TGraph*> mine_g; double zmax = -1e30;
        draw_overlay(m, l, /*lw=*/2, /*ms=*/1.4, mine_g, zmax);
        if (zmax <= 0.) zmax = 1.;
        // frame top set so the tallest curve reaches 70% of the frame height
        TH1F* fr = gPad->DrawFrame(0.05, 0., 200., zmax/0.70);
        fr->GetXaxis()->SetTitle("c#kern[0.25]{#tau} [mm]");
        fr->GetYaxis()->SetTitle("Asymptotic significance  #it{Z}");
        fr->GetXaxis()->SetTitleOffset(1.1); fr->GetYaxis()->SetTitleOffset(1.5);
        // the 2*NF graphs just created by draw_overlay are the last keep entries
        for (size_t i=keep.size()-2*NF; i<keep.size(); ++i)
            if (auto* g = dynamic_cast<TGraph*>(keep[i])) g->Draw("LP");

        // FPR working-point legend (top-right)
        TLegend* lg2 = new TLegend(0.71,0.69,0.94,0.90); keep.push_back(lg2);
        lg2->SetFillStyle(0); lg2->SetBorderSize(0); lg2->SetTextSize(0.026);
        lg2->SetHeader("BDT working point:");
        for (int f=0; f<NF; ++f) lg2->AddEntry(mine_g[f], Form("FPR=10^{-%d}",f+1), "lp");
        lg2->Draw();

        // signal info (+ lxy) + author-style key (top-left); grow box if lxy shown
        double y0 = lxy_label ? 0.57 : 0.60;
        TLegend* lgs = new TLegend(0.17,y0,0.50,0.90); keep.push_back(lgs);
        lgs->SetFillStyle(0); lgs->SetBorderSize(0); lgs->SetTextSize(0.026);
        lgs->AddEntry((TObject*)nullptr, "Scenario A", "");
        lgs->AddEntry((TObject*)nullptr, Form("m_{#bar{#pi}} = %d GeV", MASS[m].mpi), "");
        lgs->AddEntry((TObject*)nullptr, Form("m_{A'} = %.2f GeV", MASS[m].mA), "");
        lgs->AddEntry((TObject*)nullptr, Form("#it{B}(A'#rightarrow#mu#mu) = %.3g", BR_MASS[m]), "");
        lgs->AddEntry((TObject*)nullptr,
                      Form("Mass window: [%.3g, %.3g] GeV",
                           MWIN_FRAC_LO*MASS[m].mA, MWIN_FRAC_HI*MASS[m].mA), "");
        if (lxy_label) lgs->AddEntry((TObject*)nullptr, Form("#font[62]{%s}", lxy_label), "");
        lgs->AddEntry(make_key(1,20), "Parking", "lp");
        lgs->AddEntry(make_key(2,24), "Scouting", "lp");
        lgs->Draw();

        // title band, right-aligned to the frame's right edge (1 - right margin)
        TLatex tg; tg.SetNDC(); tg.SetTextSize(0.040); tg.SetTextAlign(32);
        tg.DrawLatex(0.96,0.955, Form("%.2f fb#kern[0.30]{^{-1}} (13.6 TeV, 2024)", LUMI_FB));
    };

    // =====================================================================
    //  PDF 1: "all lxy" -- one page per mass point (single plot per page)
    // =====================================================================
    {
        TString pdf = out_dir + "/" + pfx + "compare_Z_Javier_alllxy_" + MODEL + ".pdf";
        for (int m=0; m<NM; ++m) {
            TCanvas cz(Form("cz%d",m), "", 1000, 800);
            cz.cd();
            draw_full_page(m, L_ALL, /*lxy_label=*/nullptr);
            TString pg = (m==0 && NM==1) ? pdf : (m==0) ? (pdf+"(") : (m==NM-1) ? (pdf+")") : pdf;
            cz.Print(pg);
        }
        printf("[compare] -> %s\n", pdf.Data());
    }

    // =====================================================================
    //  PDF 2: lxy partitions -- SAME single-plot format as "all lxy", one plot
    //  per page: NM mass points x 3 lxy categories (0-1,1-10,10-100 cm), all in
    //  one PDF. (3-per-page couldn't reproduce the full-legend style cleanly.)
    // =====================================================================
    {
        TString pdf = out_dir + "/" + pfx + "compare_Z_Javier_by_lxy_" + MODEL + ".pdf";
        const int total = NM * 3;
        int idx = 0;
        for (int m=0; m<NM; ++m) {
            for (int l=0; l<3; ++l, ++idx) {
                TCanvas cv(Form("cv%d_%d",m,l), "", 1000, 800);
                cv.cd();
                draw_full_page(m, l, LXY_TEX[l]);
                TString pg = (total==1) ? pdf : (idx==0) ? (pdf+"(")
                           : (idx==total-1) ? (pdf+")") : pdf;
                cv.Print(pg);
            }
        }
        printf("[compare] -> %s (%d pages)\n", pdf.Data(), total);
    }
}
