#include "ROOT/RDataFrame.hxx"
#include "TChain.h"
#include "TH1D.h"
#include <fstream>
#include <vector>
#include <string>
#include <map>

void threshold_from_deploy_all() {
   ROOT::EnableImplicitMT();

   const TString BASE = "/home/hep/jtafoyav/vols/parking/bdt/icenet/output/dqcd/deploy";
   const TString MID  = "gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/tafoyava/samples/bParking/2024";
   const char* OUT = "/vols/cms/jtafoyav/parking/bdt/icenet/_tools/models/2024_2026-06-13_fourmuonSV_charge0/getBDTScoreThreshold_fromDeploymentOutput.csv";

   std::vector<std::string> QCD = {
     "QCD_Bin-PT-15to20","QCD_Bin-PT-20to30","QCD_Bin-PT-30to50","QCD_Bin-PT-50to80",
     "QCD_Bin-PT-80to120","QCD_Bin-PT-120to170","QCD_Bin-PT-170to300","QCD_Bin-PT-300to470",
     "QCD_Bin-PT-470to600","QCD_Bin-PT-600to800","QCD_Bin-PT-800to1000","QCD_Bin-PT-1000"};
   std::vector<double> XS = {2799000.0,2526000.0,1362000.0,376600.0,88930.0,21230.0,7055.0,619.3,59.24,18.21,3.275,1.078};
   const TString SUF = "_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8";
   const double LUMI = 33.6*1000.0;
   const int NB = 10000; const double TARGET = 1e-4;

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

   std::ofstream csv(OUT);
   csv << "run,model,threshold,bkg_eff_at_threshold\n";

   for (auto& mtp : MT) {
      std::string lbl = mtp.first, mt = mtp.second;
      std::map<std::string,TH1D*> tot; for (auto& b : BR) tot[b]=nullptr;

      for (size_t k=0;k<QCD.size();++k) {
         TChain* ch = new TChain();
         TString patt = TString::Format("%s/modeltag__%s/%s/%s%s/*/*/*icenet.root/Events",
                          BASE.Data(), mt.c_str(), MID.Data(), QCD[k].c_str(), SUF.Data());
         ch->Add(patt);
         ROOT::RDataFrame df(*ch);
         auto cnt = df.Count();
         std::map<std::string, ROOT::RDF::RResultPtr<TH1D>> hh;
         for (auto& b : BR)
            hh[b] = df.Histo1D({Form("h_%s_%zu_%s",lbl.c_str(),k,b.c_str()),"",NB,0.0,1.0}, b.c_str());
         double n = (double)*cnt;
         if (n<=0) { printf("[warn] %s %s: 0 events\n", lbl.c_str(), QCD[k].c_str()); delete ch; continue; }
         for (auto& b : BR) {
            TH1D* h = (TH1D*)hh[b]->Clone(Form("c_%s_%zu_%s",lbl.c_str(),k,b.c_str()));
            h->SetDirectory(0); h->Scale(LUMI*XS[k]/n);
            if (!tot[b]) { tot[b]=(TH1D*)h->Clone(Form("tot_%s_%s",lbl.c_str(),b.c_str())); tot[b]->SetDirectory(0); }
            else tot[b]->Add(h);
            delete h;
         }
         delete ch;
      }

      for (auto& b : BR) {
         TH1D* t = tot[b]; double thr=-1.0, eff=-1.0;
         if (t && t->Integral()>0) {
            t->Scale(1.0/t->Integral());
            TH1D* rc = (TH1D*)t->GetCumulative(kFALSE);
            for (int i=1;i<=NB;i++){ if (rc->GetBinContent(i)<=TARGET){ thr=rc->GetBinLowEdge(i); eff=rc->GetBinContent(i); break; } }
         }
         printf("%-12s %-14s threshold=%.4f  bkg_eff=%.2e\n", lbl.c_str(), b.c_str(), thr, eff);
         csv << lbl << "," << b << "," << Form("%.4f",thr) << "," << Form("%.2e",eff) << "\n";
      }
   }
   csv.close();
   printf("WROTE %s\n", OUT);
}
