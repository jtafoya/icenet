// ===========================================================================
//  getROC_render.C  --  BDT ROC curves per lxy (displacement) bin.
//
//  One PDF, three pages (one per displacement bin: 0-1, 1-10, 10-100 cm).
//  Per page:
//    * the BDT ROC curve  -- pooled scenario A signal (the equal-per-sample
//      training mix) vs xs-weighted QCD, both restricted to that lxy bin.
//      Drawn SOLID where the weighted QCD tail has N_eff >= NEFF_MIN effective
//      events and DASHED below that, where the curve is MC-stat limited and
//      should not be read quantitatively.
//    * the four scenario A lifetimes of mpi=10, mA=1.00 marked at the
//      FPR=1e-4 working point (colour = ctau).
//
//  Input CSVs are produced by getROC_extract.py from the icenet eval output
//  (Scenario A, SingleORDouble trigger, model XGB-NOJETS).
//
//  Usage:  root -l -b -q 'getROC_render.C()'        # log FPR axis
//          root -l -b -q 'getROC_render.C("linx")'  # linear FPR axis
//                                                   # (-> ..._linx.pdf)
// ===========================================================================
#include "TCanvas.h"
#include "TPad.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>

static const int NL  = 3;   // lxy displacement bins (one per-bin page each)
static const int NLL = 4;   // + inclusive, for the overlay page
static const int NT  = 4;   // lifetimes

static const char* LXY_CSV[NLL] = { "lxy[0,1]cm", "lxy[1,10]cm", "lxy[10,100]cm",
                                    "inclusive" };
static const char* LXY_TEX[NLL] = { "l_{xy}#in[0,1] cm", "l_{xy}#in[1,10] cm",
                                    "l_{xy}#in[10,100] cm", "inclusive" };
// overlay colours per lxy bin (inclusive = black)
static const int   LXY_COL[NLL] = { kBlue+1, kGreen+2, kRed+1, kBlack };
static const double CTAU[NT]    = { 0.1, 1.0, 10.0, 100.0 };

// Below this many EFFECTIVE background events the weighted QCD tail is carried
// by one or two high-xs events -> curve drawn dashed (see getROC_extract.py).
static const double NEFF_MIN = 10.0;
static const double FPR_WP   = 1e-4;   // analysis working point

// Split a CSV line into fields, honouring "..." quoting (the lxy tags contain
// a comma, e.g. "lxy[0,1]cm").
static std::vector<std::string> csv_split(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    bool q = false;
    for (char ch : line) {
        if (ch == '"')       q = !q;
        else if (ch == ',' && !q) { out.push_back(cur); cur.clear(); }
        else                 cur += ch;
    }
    out.push_back(cur);
    return out;
}

void getROC_render(const char* opt = "")
{
    gStyle->SetOptStat(0);
    TString pfx = TString(opt).Contains("test") ? "test_" : "";
    const Bool_t LOGX = !TString(opt).Contains("linx");
    const char* MODEL = "Mu10ORDoubleMu";

    TString script_dir = gSystem->DirName(__FILE__);
    TString out_dir    = script_dir + "/output_getSignalEff_noMET";

    // ---- read the ROC curves (3 displacement bins + inclusive) -------------
    std::vector<double> fpr[NLL], tpr[NLL], neff[NLL];
    {
        TString f = out_dir + "/" + pfx + "roc_curve_" + MODEL + ".csv";
        std::ifstream in(f.Data());
        if (!in) { printf("[roc] ERROR: cannot open %s\n", f.Data()); return; }
        std::string line; std::getline(in, line);           // header
        int nrow = 0;
        while (std::getline(in, line)) {
            if (line.empty()) continue;
            auto v = csv_split(line);
            for (int c = 0; c < NLL; ++c)
                if (v[0] == LXY_CSV[c]) {
                    fpr [c].push_back(std::stod(v[1]));
                    tpr [c].push_back(std::stod(v[2]));
                    neff[c].push_back(std::stod(v[5]));
                    ++nrow;
                }
        }
        printf("[roc] read %s : %d curve points\n", f.Data(), nrow);
    }

    // ---- read the lifetime working points -----------------------------------
    double p_fpr[NL][NT] = {}, p_tpr[NL][NT] = {}, p_thr[NL] = {};
    double p_neff[NL] = {};  int p_nsig[NL][NT] = {}, p_nraw[NL] = {};
    {
        TString f = out_dir + "/" + pfx + "roc_points_" + MODEL + ".csv";
        std::ifstream in(f.Data());
        if (!in) { printf("[roc] ERROR: cannot open %s\n", f.Data()); return; }
        std::string line; std::getline(in, line);           // header
        while (std::getline(in, line)) {
            if (line.empty()) continue;
            auto v = csv_split(line);
            for (int c = 0; c < NL; ++c) {
                if (v[0] != LXY_CSV[c]) continue;
                for (int t = 0; t < NT; ++t) {
                    if (std::fabs(std::stod(v[1]) - CTAU[t]) > 1e-6 * CTAU[t]) continue;
                    p_fpr [c][t] = std::stod(v[2]);
                    p_tpr [c][t] = std::stod(v[3]);
                    p_thr [c]    = std::stod(v[4]);
                    p_nraw[c]    = std::stoi(v[5]);
                    p_neff[c]    = std::stod(v[6]);
                    p_nsig[c][t] = std::stoi(v[7]);
                }
            }
        }
    }

    // ---- styling (colour = lifetime; matches the compare/significance plots)
    const int COL[NT] = { kBlue+1, kRed+1, kGreen+2, kMagenta+1 };
    const int MKR[NT] = { 20, 21, 22, 33 };

    std::vector<TObject*> keep;
    TString pdf = out_dir + "/" + pfx + "roc_by_lxy_" + MODEL
                + (LOGX ? "" : "_linx") + ".pdf";

    for (int c = 0; c < NL; ++c) {
        TCanvas cv(Form("cr%d", c), "", 1000, 800);
        cv.cd();
        if (LOGX) gPad->SetLogx();
        gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.12);
        gPad->SetRightMargin(0.04);

        TH1F* fr = gPad->DrawFrame(LOGX ? 1e-5 : 0., 0., 1., 1.0);
        fr->GetXaxis()->SetTitle("False positive rate (QCD efficiency)");
        fr->GetYaxis()->SetTitle("True positive rate (signal efficiency)");
        fr->GetXaxis()->SetTitleOffset(1.1); fr->GetYaxis()->SetTitleOffset(1.5);

        // ROC split into a reliable (solid) and an MC-stat-limited (dashed) part
        TGraph* g_ok = new TGraph(); keep.push_back(g_ok);
        TGraph* g_lo = new TGraph(); keep.push_back(g_lo);
        for (size_t i = 0; i < fpr[c].size(); ++i) {
            TGraph* g = (neff[c][i] >= NEFF_MIN) ? g_ok : g_lo;
            g->SetPoint(g->GetN(), fpr[c][i], tpr[c][i]);
        }
        // bridge the gap so the two segments join visually
        if (g_lo->GetN() && g_ok->GetN()) {
            double x, y; g_ok->GetPoint(0, x, y);
            g_lo->SetPoint(g_lo->GetN(), x, y);
        }
        g_lo->SetLineColor(kGray+2); g_lo->SetLineWidth(2); g_lo->SetLineStyle(2);
        g_ok->SetLineColor(kBlack);  g_ok->SetLineWidth(2); g_ok->SetLineStyle(1);
        g_lo->Draw("L"); g_ok->Draw("L");

        // working-point guide line
        TLine* l = new TLine(FPR_WP, 0., FPR_WP, 1.0); keep.push_back(l);
        l->SetLineStyle(3); l->SetLineColor(kGray+1); l->Draw();

        // the four lifetimes at the working point
        TGraph* gp[NT];
        for (int t = 0; t < NT; ++t) {
            gp[t] = new TGraph(1); keep.push_back(gp[t]);
            gp[t]->SetPoint(0, p_fpr[c][t], p_tpr[c][t]);
            gp[t]->SetMarkerColor(COL[t]); gp[t]->SetLineColor(COL[t]);
            gp[t]->SetMarkerStyle(MKR[t]); gp[t]->SetMarkerSize(1.6);
            gp[t]->Draw("P");
        }

        // lifetime key: the ROC runs through the top-right on a log axis and
        // hugs the top-left on a linear one, so the bottom-right is free in both
        TLegend* lg = LOGX ? new TLegend(0.68, 0.25, 0.93, 0.47)
                           : new TLegend(0.72, 0.22, 0.94, 0.44);
        keep.push_back(lg);
        lg->SetFillStyle(0); lg->SetBorderSize(0); lg->SetTextSize(0.026);
        lg->SetHeader(Form("at FPR = 10^{-4}:"));
        for (int t = 0; t < NT; ++t)
            lg->AddEntry(gp[t], Form("c#tau = %g mm", CTAU[t]), "p");
        lg->Draw();

        // signal info + curve key (top-left on log; pushed under the curve on linear)
        TLegend* ls = LOGX ? new TLegend(0.17, 0.53, 0.50, 0.90)
                           : new TLegend(0.38, 0.50, 0.71, 0.87);
        keep.push_back(ls);
        ls->SetFillStyle(0); ls->SetBorderSize(0); ls->SetTextSize(0.026);
        ls->AddEntry((TObject*)nullptr, "Scenario A", "");
        ls->AddEntry((TObject*)nullptr, "m_{#bar{#pi}} = 10 GeV, m_{A'} = 1.00 GeV", "");
        ls->AddEntry((TObject*)nullptr, Form("#font[62]{%s}", LXY_TEX[c]), "");
        ls->AddEntry(g_ok, "BDT (all scenario A vs QCD)", "l");
        ls->AddEntry(g_lo, Form("same, N_{eff}^{QCD} < %g (stat. limited)", NEFF_MIN), "l");
        ls->Draw();

        // caveats that change how the markers should be read. Kept inside the
        // frame (bottom margin is 0.12) and below the lifetime key.
        TLatex tn; tn.SetNDC(); tn.SetTextSize(0.023); tn.SetTextAlign(12);
        double yn = 0.205 + 0.035;
        if (p_fpr[c][0] > 1.5 * FPR_WP)
            tn.DrawLatex(0.17, yn -= 0.035, Form(
                "#color[2]{FPR = 10^{-4} unreachable: highest-score QCD event alone gives %.1e}",
                p_fpr[c][0]));
        if (p_nsig[c][0] < 10)
            tn.DrawLatex(0.17, yn -= 0.035, Form(
                "#color[2]{c#tau = 0.1 mm has only %d signal events in this bin}", p_nsig[c][0]));

        TLatex tg; tg.SetNDC(); tg.SetTextSize(0.040); tg.SetTextAlign(32);
        tg.DrawLatex(0.96, 0.955, "13.6 TeV, 2024");

        // per-bin pages open the multipage PDF; the overlay page (below) closes it
        TString pg = (c == 0) ? (pdf + "(") : pdf;
        cv.Print(pg);
    }

    // ---- overlay page: all lxy bins (+ inclusive) on one axes ---------------
    {
        TCanvas cv("cr_all", "", 1000, 800);
        cv.cd();
        if (LOGX) gPad->SetLogx();
        gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.12);
        gPad->SetRightMargin(0.04);

        TH1F* fr = gPad->DrawFrame(LOGX ? 1e-5 : 0., 0., 1., 1.0);
        fr->GetXaxis()->SetTitle("False positive rate (QCD efficiency)");
        fr->GetYaxis()->SetTitle("True positive rate (signal efficiency)");
        fr->GetXaxis()->SetTitleOffset(1.1); fr->GetYaxis()->SetTitleOffset(1.5);

        TLine* l = new TLine(FPR_WP, 0., FPR_WP, 1.0); keep.push_back(l);
        l->SetLineStyle(3); l->SetLineColor(kGray+1); l->Draw();

        TLegend* lg = new TLegend(0.62, 0.16, 0.93, 0.40); keep.push_back(lg);
        lg->SetFillStyle(0); lg->SetBorderSize(0); lg->SetTextSize(0.028);
        lg->SetHeader("Displacement bin:");

        for (int c = 0; c < NLL; ++c) {
            // split each curve into a reliable (solid) and stat-limited (dashed) part
            TGraph* g_ok = new TGraph(); keep.push_back(g_ok);
            TGraph* g_lo = new TGraph(); keep.push_back(g_lo);
            for (size_t i = 0; i < fpr[c].size(); ++i)
                ((neff[c][i] >= NEFF_MIN) ? g_ok : g_lo)
                    ->SetPoint(((neff[c][i] >= NEFF_MIN) ? g_ok : g_lo)->GetN(),
                               fpr[c][i], tpr[c][i]);
            if (g_lo->GetN() && g_ok->GetN()) {
                double x, y; g_ok->GetPoint(0, x, y);
                g_lo->SetPoint(g_lo->GetN(), x, y);   // bridge the gap
            }
            g_lo->SetLineColor(LXY_COL[c]); g_lo->SetLineWidth(2); g_lo->SetLineStyle(2);
            g_ok->SetLineColor(LXY_COL[c]); g_ok->SetLineWidth(3); g_ok->SetLineStyle(1);
            g_lo->Draw("L"); g_ok->Draw("L");
            lg->AddEntry(g_ok, LXY_TEX[c], "l");
        }
        lg->Draw();

        TLegend* ls = new TLegend(0.17, 0.74, 0.50, 0.90); keep.push_back(ls);
        ls->SetFillStyle(0); ls->SetBorderSize(0); ls->SetTextSize(0.028);
        ls->AddEntry((TObject*)nullptr, "Scenario A (all points vs QCD)", "");
        ls->AddEntry((TObject*)nullptr, "#font[62]{by displacement bin}", "");
        ls->AddEntry((TObject*)nullptr, Form("solid: N_{eff}^{QCD} #geq %g; dashed: below", NEFF_MIN), "");
        ls->Draw();

        TLatex tg; tg.SetNDC(); tg.SetTextSize(0.040); tg.SetTextAlign(32);
        tg.DrawLatex(0.96, 0.955, "13.6 TeV, 2024");

        cv.Print(pdf + ")");   // close the multipage PDF
    }
    printf("[roc] -> %s (%d pages)\n", pdf.Data(), NL + 1);
}
