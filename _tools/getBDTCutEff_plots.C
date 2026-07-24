// ===========================================================================
//  getBDTCutEff_plots.C  --  "BDT cut efficiency" plots from deployment output,
//  for the 6 latest BDTs (2024_2026-06-13_fourmuonSV_charge0).
//
//  For each of the 6 configs it builds the summed, xs-weighted QCD BDT-score
//  histogram (12 pT bins), normalises it, and takes the reverse cumulative =
//  BACKGROUND EFFICIENCY (FPR) vs BDT score cut. That curve is what the analysis
//  working point is read off. Only the QCD (background) curve is drawn; a signal
//  curve would need signal deployment output, which exists for scenarioA but NOT
//  scenarioB1, so the 6 plots are kept uniform.
//
//  CACHING: the reverse-cumulative histograms are the expensive part (reads all
//  deployment output). They are cached to BDTcutEff_cache.root, so styling edits
//  re-plot in seconds:
//     root -l -b -q 'getBDTCutEff_plots.C'         # read deployment, write cache
//     root -l -b -q 'getBDTCutEff_plots.C(true)'   # re-plot from cache (fast)
// ===========================================================================
#include "ROOT/RDataFrame.hxx"
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <array>

void getBDTCutEff_plots(bool use_cache=false) {
   ROOT::EnableImplicitMT();
   gStyle->SetOptStat(0);

   const TString BASE = "/home/hep/jtafoyav/vols/parking/bdt/icenet/output/dqcd/deploy";
   const TString MID  = "gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/tafoyava/samples/bParking/2024";
   const TString ODIR = "/vols/cms/jtafoyav/parking/bdt/icenet/_tools/models/2024_2026-06-13_fourmuonSV_charge0";
   const TString CACHE = ODIR + "/BDTcutEff_cache.root";

   std::vector<std::string> QCD = {
     "QCD_Bin-PT-15to20","QCD_Bin-PT-20to30","QCD_Bin-PT-30to50","QCD_Bin-PT-50to80",
     "QCD_Bin-PT-80to120","QCD_Bin-PT-120to170","QCD_Bin-PT-170to300","QCD_Bin-PT-300to470",
     "QCD_Bin-PT-470to600","QCD_Bin-PT-600to800","QCD_Bin-PT-800to1000","QCD_Bin-PT-1000"};
   std::vector<double> XS = {2799000.0,2526000.0,1362000.0,376600.0,88930.0,21230.0,7055.0,619.3,59.24,18.21,3.275,1.078};
   const TString SUF = "_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8";
   const double LUMI = 33.6*1000.0;
   const int NB = 10000;
   const double FPR_WP[4] = {1e-1, 1e-2, 1e-3, 1e-4};

   std::vector<std::pair<std::string,std::string>> MT = {
     // {"A_combined", "scenarioA_all_no_DA_old_BDT_with_dRmSV_2024"},
     {"A_combined", "scenarioA_all_no_DA_old_BDT_2024"},
     // {"A_DoubleMu", "scenarioA_all_no_DA_old_BDT_with_dRmSV_2024_DoubleMu"},
     {"A_DoubleMu", "scenarioA_all_no_DA_old_BDT_2024_DoubleMu"},
     // {"A_Mu10",     "scenarioA_all_no_DA_old_BDT_with_dRmSV_2024_Mu10"},
     {"A_Mu10",     "scenarioA_all_no_DA_old_BDT_2024_Mu10"},
     // {"B1_combined","scenarioB1_all_no_DA_old_BDT_with_dRmSV_2024"},
     {"B1_combined","scenarioB1_all_no_DA_old_BDT_2024"},
     // {"B1_DoubleMu","scenarioB1_all_no_DA_old_BDT_with_dRmSV_2024_DoubleMu"},
     {"B1_DoubleMu","scenarioB1_all_no_DA_old_BDT_2024_DoubleMu"},
     // {"B1_Mu10",    "scenarioB1_all_no_DA_old_BDT_with_dRmSV_2024_Mu10"}};
     {"B1_Mu10",    "scenarioB1_all_no_DA_old_BDT_2024_Mu10"}};
   std::vector<std::string> BR = {"xgb01","xgb01_NOJETS"};

   TFile* fcache = use_cache ? TFile::Open(CACHE, "READ") : TFile::Open(CACHE, "RECREATE");
   if (!fcache || fcache->IsZombie()) { printf("[cuteff] ERROR: cache %s\n", CACHE.Data()); return; }

   std::ofstream csv((ODIR + "/thresholds_by_working_point_deploy.csv").Data());
   csv << "run,model,thr_FPR_1e-1,thr_FPR_1e-2,thr_FPR_1e-3,thr_FPR_1e-4\n";

   for (auto& mtp : MT) {
      std::string lbl = mtp.first, mt = mtp.second;
      std::map<std::string,TH1D*> rc;   // reverse-cumulative FPR-vs-score per model

      for (auto& b : BR) {
         TString hname = Form("rc_%s_%s", lbl.c_str(), b.c_str());
         if (use_cache) {
            TH1D* c = (TH1D*)fcache->Get(hname);
            if (c) c->SetDirectory(0);
            rc[b] = c;
            continue;
         }
         // compute from deployment output (xs-weighted, summed over pT bins)
         TH1D* tot = nullptr;
         for (size_t k=0;k<QCD.size();++k) {
            TChain ch;
            ch.Add(TString::Format("%s/modeltag__%s/%s/%s%s/*/*/*icenet.root/Events",
                     BASE.Data(), mt.c_str(), MID.Data(), QCD[k].c_str(), SUF.Data()));
            ROOT::RDataFrame df(ch);
            auto cnt = df.Count();
            auto h = df.Histo1D({Form("h_%s_%zu_%s",lbl.c_str(),k,b.c_str()),"",NB,0.0,1.0}, b.c_str());
            double n = (double)*cnt;
            if (n<=0) { printf("[warn] %s %s %s: 0 events\n", lbl.c_str(), QCD[k].c_str(), b.c_str()); continue; }
            TH1D* hc = (TH1D*)h->Clone(); hc->SetDirectory(0); hc->Scale(LUMI*XS[k]/n);
            if (!tot) { tot=(TH1D*)hc->Clone(Form("tot_%s_%s",lbl.c_str(),b.c_str())); tot->SetDirectory(0); }
            else tot->Add(hc);
            delete hc;
         }
         if (!tot || tot->Integral()<=0) { printf("[warn] %s %s: no QCD\n", lbl.c_str(), b.c_str()); rc[b]=nullptr; continue; }
         tot->Scale(1.0/tot->Integral());
         TH1D* c = (TH1D*)tot->GetCumulative(kFALSE);   // reverse cumulative = FPR(>score)
         c->SetName(hname); c->SetDirectory(0); rc[b]=c; delete tot;
         fcache->cd(); c->Write(hname);                 // cache for fast re-plotting
      }

      // thresholds at each FPR working point (per model)
      std::map<std::string, std::array<double,4>> thr;
      for (auto& b : BR) {
         thr[b] = {-1.0,-1.0,-1.0,-1.0};
         if (!rc[b]) continue;
         for (int w=0;w<4;++w)
            for (int i=1;i<=NB;i++)
               if (rc[b]->GetBinContent(i) <= FPR_WP[w]) { thr[b][w]=rc[b]->GetBinLowEdge(i); break; }
         csv << lbl << "," << b;
         for (int w=0;w<4;++w) csv << "," << Form("%.4f", thr[b][w]);
         csv << "\n";
      }

      // ---- plot: FPR vs BDT score, NOJETS (deployed) model only ----
      std::string scen = (lbl.rfind("B1",0)==0) ? "Scenario B1" : "Scenario A";
      std::string trig = lbl.substr(lbl.find('_')+1);
      TString title = TString::Format("%s BDT (%s)", scen.c_str(), trig.c_str());

      const std::string BM = "xgb01_NOJETS";
      TH1D* cN = rc[BM];
      if (cN) {
         TCanvas cv("cv","",900,700);
         cv.SetLogy();                                  // linear x, log y
         cv.SetLeftMargin(0.13); cv.SetRightMargin(0.04);
         cN->SetLineColor(kBlue+1); cN->SetLineWidth(2); cN->SetTitle("");
         cN->GetXaxis()->SetTitle("BDT score");
         cN->GetYaxis()->SetTitle("QCD efficiency (FPR)");
         cN->GetYaxis()->SetTitleOffset(1.3);
         cN->GetXaxis()->SetRangeUser(0.0, 1.0);
         cN->GetYaxis()->SetRangeUser(3e-5, 9.0);
         cN->Draw("hist");

         // dashed threshold markers: for each working point, a vertical line
         // from the x-axis and a horizontal line from the y-axis, both meeting
         // the curve at (thr[w], FPR_WP[w]) -- the curve value at thr[w] is
         // exactly FPR_WP[w], so neither line runs past the curve.
         const double YMIN = 3e-5;
         for (int w=0;w<4;++w) {
            double x = thr[BM][w]; if (x<=0) continue;
            TLine* vl = new TLine(x, YMIN, x, FPR_WP[w]);      // x-axis -> curve
            vl->SetLineStyle(2); vl->SetLineColor(kGray+2); vl->SetLineWidth(2); vl->Draw();
            TLine* hl = new TLine(0.0, FPR_WP[w], x, FPR_WP[w]); // y-axis -> curve
            hl->SetLineStyle(2); hl->SetLineColor(kGray+2); hl->SetLineWidth(2); hl->Draw();
         }

         // header as a legend with a semi-transparent white background (drawn
         // last so it sits over the dashed lines but the lines show through);
         // top-right, in the empty region above the curve (the y-max=9 headroom)
         TLegend* lg = new TLegend(0.53, 0.71, 0.95, 0.88);
         lg->SetFillColorAlpha(kWhite, 0.6); lg->SetFillStyle(1001);
         lg->SetBorderSize(0); lg->SetTextAlign(12); lg->SetMargin(0.03);
         lg->SetTextSize(0.028);
         lg->SetHeader(title);
         lg->AddEntry((TObject*)nullptr, "thr @ FPR 10^{-1}/10^{-2}/10^{-3}/10^{-4}:", "");
         lg->AddEntry((TObject*)nullptr, Form("%.4f / %.4f / %.4f / %.4f",
             thr[BM][0], thr[BM][1], thr[BM][2], thr[BM][3]), "");
         lg->Draw();

         // energy/year top-right, non-bold (font 42 = regular helvetica)
         TLatex tr; tr.SetNDC(); tr.SetTextFont(42); tr.SetTextSize(0.034); tr.SetTextAlign(32);
         tr.DrawLatex(0.96,0.945, "(13.6 TeV, 2024)");

         TString out = ODIR + "/BDTcutEff_" + lbl.c_str() + "_xgb01_NOJETS.pdf";
         cv.SaveAs(out);
         printf("[cuteff] -> %s\n", out.Data());
      }
      if (use_cache) for (auto& b : BR) if (rc[b]) delete rc[b];
   }
   csv.close();
   if (!use_cache) { fcache->Close(); printf("[cuteff] wrote cache %s\n", CACHE.Data()); }
   printf("[cuteff] wrote thresholds CSV\n");
}
