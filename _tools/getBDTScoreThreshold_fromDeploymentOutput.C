#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TStyle.h"
#include <iostream>
#include <cmath>

// Based on a copy of Kai's nanotools_analysis_code_bdt_background_new_models.cpp

using namespace ROOT::VecOps;

void hist_stack(const char* stack_name, ROOT::RDF::RResultPtr<::TH1D > h1 , ROOT::RDF::RResultPtr<::TH1D > h, std::string name, bool normalize) {
     TCanvas* c1 = new TCanvas("", "", 800, 700);
     h->SetTitle("");
     if (normalize == true) {
        h->SetMaximum(1);
     }
     else {
        h->SetMaximum(1E11);
        h->SetMinimum(1E-1);
     }
     h1->SetLineWidth(2);

     h1->SetLineColor(kRed);

     h->SetLineColor(kBlue+1);
     h->SetLineWidth(2);

     h->DrawClone("hist");
     h1->DrawClone("same hist");
     TLegend* out_legend = new TLegend(0.62, 0.70, 0.82, 0.88);
     out_legend->SetFillColor(0);
     out_legend->SetBorderSize(0);
     out_legend->SetTextSize(0.025);
     out_legend->AddEntry(h1->GetName(), "QCD Background", "f");
     out_legend->AddEntry(h->GetName(), "#splitline{Signal: Scenario A,}{#splitline{m_{#pi} = 4 GeV, m_{A}=1.33 GeV, c #tau = 1 mm,}{#xi_{#omega} = 1, #xi_{#Lambda} = 1}}", "f");
     out_legend->Draw("Same");
     gPad->SetLogy();

     TLatex cms_label;
     cms_label.SetTextSize(0.04);
     cms_label.DrawLatexNDC(0.10, 0.92, "#bf{CMS} #it{Work in progress}");

     TLatex header;
     header.SetTextSize(0.03);
     header.DrawLatexNDC(0.63, 0.92, "#sqrt{s} = 13 TeV, L_{int} = 41.6 fb^{-1}");

     std::string hname = h->GetName();
     std::string fname = hname + "_" + name + ".pdf";
     c1->SaveAs(fname.c_str());
     }
void hist_stack_old(const char* stack_name, std::vector<TH1D*> hist_vector) {

	
     TCanvas* c1 = new TCanvas();
     const char* stack_title = hist_vector.at(0)->GetTitle();
     THStack* out_stack = new THStack(stack_name, stack_title);
     out_stack->Add(hist_vector.at(0));
     out_stack->Add(hist_vector.at(1),"S");
     out_stack->Draw("hist nostack");
     out_stack->GetXaxis()->SetTitle(hist_vector.at(0)->GetXaxis()->GetTitle());
     out_stack->GetYaxis()->SetTitle(hist_vector.at(0)->GetYaxis()->GetTitle());
     TLegend* out_legend = new TLegend(0.78, 0.695, 0.98, 0.775);
     out_legend->AddEntry(hist_vector.at(0)->GetName(), "Background", "l");
     out_legend->AddEntry(hist_vector.at(1)->GetName(), "Signal", "l");
     out_legend->Draw("Same");
     gPad->SetLogy();
    }

double Z(Double_t s, Double_t b){
    Double_t Z_score;
    if(b>0.0 and 2*((s+b)*log(1+s/b)-s) > 0.0){
        Z_score = sqrt(2*((s+b)*log(1+s/b)-s));
    }
    else{
        Z_score = 0.0;
    }

    return Z_score;
}

double Z_error(Double_t s, Double_t b){
    Double_t error7;

    Double_t Z2 = 2*((s+b)*log(1+s/b)-s);
    if(Z2>0.0 and b>0.0){
        Double_t Z = sqrt(Z2);
        error7 = sqrt(pow(log(1+s/b)*sqrt(s)/Z, 2)+pow((b*log(1+s/b)-s)*sqrt(b)/(b*Z),2));
    }
    else{
        error7 = 0.0;
    }


    return error7;
}

TH1D* make_hist(const char* hist_name, const char* hist_title, Int_t nbins, Double_t bin_start, Double_t bin_end, const char* xaxis, const char* yaxis, Color_t lcolor, Width_t lwidth) {
        TH1D* out_hist = new TH1D(hist_name, hist_title, nbins, bin_start, bin_end);
        out_hist->SetLineColor(lcolor);
        out_hist->SetLineWidth(lwidth);
        out_hist->GetXaxis()->SetTitle(xaxis);
        out_hist->GetYaxis()->SetTitle(yaxis);
        return out_hist;

}

double ssqrtb(Double_t s, Double_t b){
    Double_t ssqrtb_score;
    if(b>0.0){
        ssqrtb_score = s/sqrt(b);
    }
    else{
        ssqrtb_score = 0.0;
    }
    return ssqrtb_score;
}

double ssqrtb_error(Double_t s, Double_t b){
    Double_t error;
    if(b>0.0){
        error = sqrt(pow(1.0/sqrt(b)*sqrt(s), 2) + pow(s/(2*pow(b, 1.5))*sqrt(b), 2));
    }
    else{
        error = 0.0;
    }
    return error;
}

void hist_draw_3(std::vector<TH1D*> hist_vector) {
   TH1D *h1 = hist_vector.at(0);
   TH1D *h2 = hist_vector.at(1);

   TCanvas *c = new TCanvas();
   c->SetCanvasSize(800, 800);

   TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
   pad1->SetBottomMargin(0);
   pad1->SetGridx();
   pad1->SetLogy();
   pad1->Draw();
   pad1->cd();


   h1->SetLineColor(kBlue+1);
   h1->SetLineWidth(2);

   h1->GetYaxis()->SetTitleSize(20);
   h1->GetYaxis()->SetTitleFont(43);
   h1->GetYaxis()->SetTitleOffset(1.55);

   h2->SetLineColor(kRed);
   h2->SetLineWidth(2);
   h1->Draw("hist");
   h2->Draw("hist same");



   h1->SetMaximum(1E9);

   TLegend* out_legend = new TLegend(0.78, 0.695, 0.98, 0.775);
   out_legend->AddEntry(h2, "Background", "l");
   out_legend->AddEntry(h1, "Signal", "l");
   out_legend->Draw("Same");

   c->cd();
   TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
   pad2->SetTopMargin(0);
   pad2->SetBottomMargin(0.2);
   pad2->SetGridx();
   pad2->SetLogy();
   pad2->Draw();
   pad2->cd();
   pad2->SetLogy();

   TH1D *h3 = hist_vector.at(2);
   TH1D *h4 = hist_vector.at(3);

   h3->SetLineColor(kBlack);
   h3->SetMarkerColor(kBlack);
   h4->SetLineColor(kGreen);
   h4->SetMarkerColor(kGreen);

   h3->Sumw2();
   h4->Sumw2();
   h3->SetStats(0);
   h4->SetStats(0);
   h3->SetMarkerStyle(20);
   h4->SetMarkerStyle(20);
   h3->SetMarkerSize(1.0);
   h4->SetMarkerSize(1.0);
   h3->Draw("ep");
   h4->Draw("ep same");
   TLegend* out_legend2 = new TLegend(0.78, 0.695, 0.98, 0.775);
   out_legend2->AddEntry((TObject*)0, "", "");
   out_legend2->AddEntry(h3, "#sqrt{q_{0,A}}", "ep");
   out_legend2->AddEntry((TObject*)0, "", "");
   out_legend2->AddEntry(h4, "s/#sqrt{b}", "ep");
   out_legend2->Draw("Same");
   h3->SetTitle("");

   h3->SetMaximum(0.0);
   h3->GetYaxis()->SetNdivisions(-505);
   h3->GetYaxis()->SetTitleSize(20);
   h3->GetYaxis()->SetTitleFont(43);
   h3->GetYaxis()->SetTitleOffset(1.55);
   h3->GetYaxis()->SetLabelFont(43);
   h3->GetYaxis()->SetLabelSize(15);

   h3->GetXaxis()->SetTitleSize(20);
   h3->GetXaxis()->SetTitleFont(43);
   h3->GetXaxis()->SetTitleOffset(4.);
   h3->GetXaxis()->SetLabelFont(43);
   h3->GetXaxis()->SetLabelSize(15);
}


void nano_analysis(){
   // Selector to run over a subset of the sample
   Bool_t runTest = kFALSE;	// run over all events (production)
   //Bool_t runTest = kTRUE;	// process only 1% of events for fast debugging

   if (!runTest){
      printf("Running over entirety of the sample\n");
   } else {
      printf("Debugging mode: running over 1%% of the sample\n");
   }

   TChain *signal_original = new TChain();
   TChain *signal = new TChain();
   TChain *background1 = new TChain();
   TChain *background2 = new TChain();
   TChain *background3 = new TChain();
   TChain *background4 = new TChain();
   TChain *background5 = new TChain();
   TChain *background6 = new TChain();
   TChain *background7 = new TChain();
   TChain *background8 = new TChain();
   TChain *background9 = new TChain();
   TChain *background10 = new TChain();
   TChain *background11 = new TChain();
   TChain *background12 = new TChain();

   TChain *background5_ext = new TChain();
   TChain *background6_ext = new TChain();
   TChain *background8_ext = new TChain();
   TChain *background9_ext = new TChain();
   TChain *background11_ext = new TChain();   

   TChain *signal_bdt = new TChain();
   TChain *signal_bdt2 = new TChain();
   TChain *background1_bdt = new TChain();
   TChain *background2_bdt = new TChain();
   TChain *background3_bdt = new TChain();
   TChain *background4_bdt = new TChain();
   TChain *background5_bdt = new TChain();
   TChain *background6_bdt = new TChain();
   TChain *background7_bdt = new TChain();
   TChain *background8_bdt = new TChain();
   TChain *background9_bdt = new TChain();
   TChain *background10_bdt = new TChain();
   TChain *background11_bdt = new TChain();
   TChain *background12_bdt = new TChain();

   TChain *background5_bdt_ext = new TChain();
   TChain *background6_bdt_ext = new TChain();
   TChain *background8_bdt_ext = new TChain();
   TChain *background9_bdt_ext = new TChain();
   TChain *background11_bdt_ext = new TChain();

   /*
   //N.B. need to uncomment each individial signal/background line (commented out with // to avoid warnings due to /*

   ///// 2018
   TString path_bdt_output_base="/home/hep/jtafoyav/vols/parking/bdt/icenet/output/dqcd/deploy";
   TString modeltag="modeltag__scenarioA_all_no_DA_old_BDT_with_dRmSV_2018_firstTry";
   TString path_bdt_qcd="gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mcitron/darkshowersamples/bparkProductionAll_V1p3/tmp";
   TString path_bdt_signal="vols/cms/khl216/bparkProductionAll_V1p3";
   TString BDTmodel="xgb01_NOJETS";

   //signal_bdt->Add( TString::Format("%s/%s/%s/scenarioA_mpi_4_mA_1p33_ctau_1p0/*.root/Events", path_bdt_output_base.Data(), modeltag.Data(), path_bdt_signal.Data() ) );
   //signal_bdt2->Add( TString::Format("%s/%s/%s/scenarioA_mpi_4_mA_1p33_ctau_10/*.root/Events", path_bdt_output_base.Data(), modeltag.Data(), path_bdt_signal.Data() ) );
   //background1_bdt->Add( TString::Format("%s/%s/%s/QCD_Pt-15To20_MuEnrichedPt5_TuneCP5_13TeV-pythia8_RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2_MINIAODSIM_v1p1_generationSync/*.root/Events", path_bdt_output_base.Data(), modeltag.Data(), path_bdt_qcd.Data() ) );
   //background2_bdt->Add( TString::Format("%s/%s/%s/QCD_Pt-20To30_MuEnrichedPt5_TuneCP5_13TeV-pythia8_RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2_MINIAODSIM_v1p1_generationSync/*.root/Events", path_bdt_output_base.Data(), modeltag.Data(), path_bdt_qcd.Data() ) );
   //background3_bdt->Add( TString::Format("%s/%s/%s/QCD_Pt-30To50_MuEnrichedPt5_TuneCP5_13TeV-pythia8_RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2_MINIAODSIM_v1p1_generationSync/*.root/Events", path_bdt_output_base.Data(), modeltag.Data(), path_bdt_qcd.Data() ) );
   //background4_bdt->Add( TString::Format("%s/%s/%s/QCD_Pt-50To80_MuEnrichedPt5_TuneCP5_13TeV-pythia8_RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2_MINIAODSIM_v1p1_generationSync/*.root/Events", path_bdt_output_base.Data(), modeltag.Data(), path_bdt_qcd.Data() ) );
   //background5_bdt->Add( TString::Format("%s/%s/%s/QCD_Pt-80To120_MuEnrichedPt5_TuneCP5_13TeV-pythia8_RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2_MINIAODSIM_v1p1_generationSync/*.root/Events", path_bdt_output_base.Data(), modeltag.Data(), path_bdt_qcd.Data() ) );
   //background6_bdt->Add( TString::Format("%s/%s/%s/QCD_Pt-120To170_MuEnrichedPt5_TuneCP5_13TeV-pythia8_RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2_MINIAODSIM_v1p1_generationSync/*.root/Events", path_bdt_output_base.Data(), modeltag.Data(), path_bdt_qcd.Data() ) );

   //background7_bdt->Add( TString::Format("%s/%s/%s/QCD_Pt-170To300_MuEnrichedPt5_TuneCP5_13TeV-pythia8_RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2_MINIAODSIM_v1p1_generationSync/*.root/Events", path_bdt_output_base.Data(), modeltag.Data(), path_bdt_qcd.Data() ) );
   //background8_bdt->Add( TString::Format("%s/%s/%s/QCD_Pt-300To470_MuEnrichedPt5_TuneCP5_13TeV-pythia8_RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2_MINIAODSIM_v1p1_generationSync/*.root/Events", path_bdt_output_base.Data(), modeltag.Data(), path_bdt_qcd.Data() ) );
   //background9_bdt->Add( TString::Format("%s/%s/%s/QCD_Pt-470To600_MuEnrichedPt5_TuneCP5_13TeV-pythia8_RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2_MINIAODSIM_v1p1_generationSync/*.root/Events", path_bdt_output_base.Data(), modeltag.Data(), path_bdt_qcd.Data() ) );
   //background10_bdt->Add( TString::Format("%s/%s/%s/QCD_Pt-600To800_MuEnrichedPt5_TuneCP5_13TeV-pythia8_RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2_MINIAODSIM_v1p1_generationSync/*.root/Events", path_bdt_output_base.Data(), modeltag.Data(), path_bdt_qcd.Data() ) );
   //background11_bdt->Add( TString::Format("%s/%s/%s/QCD_Pt-800To1000_MuEnrichedPt5_TuneCP5_13TeV-pythia8_RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2_MINIAODSIM_v1p1_generationSync/*.root/Events", path_bdt_output_base.Data(), modeltag.Data(), path_bdt_qcd.Data() ) );
   //background12_bdt->Add( TString::Format("%s/%s/%s/QCD_Pt-1000_MuEnrichedPt5_TuneCP5_13TeV-pythia8_RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2_MINIAODSIM_v1p1_generationSync/*.root/Events", path_bdt_output_base.Data(), modeltag.Data(), path_bdt_qcd.Data() ) );

   //signal->Add( "/vols/cms/khl216/nano_out/scenario_A_no_conditional/bparkProductionAll_V1p3/scenarioA_mpi_4_mA_1p33_ctau_10/*.root/Friends" );
   //background1->Add( "/vols/cms/khl216/nano_out/scenario_A_no_conditional/bparkProductionAll_V1p3/QCD_Pt-15To20_MuEnrichedPt5_TuneCP5_13TeV-pythia8_RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2_MINIAODSIM_v1p1_generationSync/*.root/Friends" );
   //background2->Add( "/vols/cms/khl216/nano_out/scenario_A_no_conditional/bparkProductionAll_V1p3/QCD_Pt-20To30_MuEnrichedPt5_TuneCP5_13TeV-pythia8_RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2_MINIAODSIM_v1p1_generationSync/*.root/Friends" );
   //background3->Add( "/vols/cms/khl216/nano_out/scenario_A_no_conditional/bparkProductionAll_V1p3/QCD_Pt-30To50_MuEnrichedPt5_TuneCP5_13TeV-pythia8_RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2_MINIAODSIM_v1p1_generationSync/*.root/Friends" );
   //background4->Add( "/vols/cms/khl216/nano_out/scenario_A_no_conditional/bparkProductionAll_V1p3/QCD_Pt-50To80_MuEnrichedPt5_TuneCP5_13TeV-pythia8_RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2_MINIAODSIM_v1p1_generationSync/*.root/Friends" );
   //background5->Add( "/vols/cms/khl216/nano_out/scenario_A_no_conditional/bparkProductionAll_V1p3/QCD_Pt-80To120_MuEnrichedPt5_TuneCP5_13TeV-pythia8_RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2_MINIAODSIM_v1p1_generationSync/*.root/Friends" );
   //background6->Add( "/vols/cms/khl216/nano_out/scenario_A_no_conditional/bparkProductionAll_V1p3/QCD_Pt-120To170_MuEnrichedPt5_TuneCP5_13TeV-pythia8_RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2_MINIAODSIM_v1p1_generationSync/*.root/Friends" );

   //background7->Add( "/vols/cms/khl216/nano_out/scenario_A_no_conditional/bparkProductionAll_V1p3/QCD_Pt-170To300_MuEnrichedPt5_TuneCP5_13TeV-pythia8_RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2_MINIAODSIM_v1p1_generationSync/*.root/Friends" );
   //background8->Add( "/vols/cms/khl216/nano_out/scenario_A_no_conditional/bparkProductionAll_V1p3/QCD_Pt-300To470_MuEnrichedPt5_TuneCP5_13TeV-pythia8_RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2_MINIAODSIM_v1p1_generationSync/*.root/Friends" );
   //background9->Add( "/vols/cms/khl216/nano_out/scenario_A_no_conditional/bparkProductionAll_V1p3/QCD_Pt-470To600_MuEnrichedPt5_TuneCP5_13TeV-pythia8_RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2_MINIAODSIM_v1p1_generationSync/*.root/Friends" );
   //background10->Add( "/vols/cms/khl216/nano_out/scenario_A_no_conditional/bparkProductionAll_V1p3/QCD_Pt-600To800_MuEnrichedPt5_TuneCP5_13TeV-pythia8_RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2_MINIAODSIM_v1p1_generationSync/*.root/Friends" );
   //background11->Add( "/vols/cms/khl216/nano_out/scenario_A_no_conditional/bparkProductionAll_V1p3/QCD_Pt-800To1000_MuEnrichedPt5_TuneCP5_13TeV-pythia8_RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2_MINIAODSIM_v1p1_generationSync/*.root/Friends" );
   //background12->Add( "/vols/cms/khl216/nano_out/scenario_A_no_conditional/bparkProductionAll_V1p3/QCD_Pt-1000_MuEnrichedPt5_TuneCP5_13TeV-pythia8_RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2_MINIAODSIM_v1p1_generationSync/*.root/Friends" );
   */

   ///// 2024

   TString path_bdt_output_base="/home/hep/jtafoyav/vols/parking/bdt/icenet/output/dqcd/deploy";
   TString path_bdt_2024="gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/tafoyava/samples/bParking/2024";

   //// noMET

   // 2024, Mu10 OR DoubleMu, noMET deployment output
   //TString modeltag="modeltag__scenarioA_all_no_DA_old_BDT_with_dRmSV_2024_noMET";
   //TString BDTmodel="xgb01_NOJETS_NOMET";

   // 2024, Mu10, noMET deployment output
   //TString modeltag="modeltag__scenarioA_all_no_DA_old_BDT_with_dRmSV_2024_Mu10_noMET";
   //TString BDTmodel="xgb01_NOJETS_NOMET";

   // 2024, DoubleMu, noMET deployment output
   //TString modeltag="modeltag__scenarioA_all_no_DA_old_BDT_with_dRmSV_2024_DoubleMu_noMET";
   //TString BDTmodel="xgb01_NOJETS_NOMET";

   //// withMET v1 (with MET_pt and MET_phi)

   // 2024, Mu10 OR DoubleMu, withMET deployment output
   //TString modeltag="modeltag__scenarioA_all_no_DA_old_BDT_with_dRmSV_2024_withMET";
   //TString BDTmodel="xgb01_NOJETS";

   // 2024, Mu10, withMET deployment output
   TString modeltag="modeltag__scenarioA_all_no_DA_old_BDT_with_dRmSV_2024_Mu10_withMET";
   TString BDTmodel="xgb01_NOJETS";

   // 2024, DoubleMu, withMET deployment output
   //TString modeltag="modeltag__scenarioA_all_no_DA_old_BDT_with_dRmSV_2024_DoubleMu_withMET";
   //TString BDTmodel="xgb01_NOJETS";


   printf("Running on %s, for result %s\n", modeltag.Data(), BDTmodel.Data());

   signal_bdt->Add(  TString::Format("%s/%s/%s/GluGluHToDarkShowers-ScenarioA_Par-ctau-1p0-mA-1p33-mpi-4_TuneCP5_13p6TeV_powheg-pythia8/*/*/*icenet.root/Events",  path_bdt_output_base.Data(), modeltag.Data(), path_bdt_2024.Data()) );
   signal_bdt2->Add( TString::Format("%s/%s/%s/GluGluHToDarkShowers-ScenarioA_Par-ctau-10-mA-1p33-mpi-4_TuneCP5_13p6TeV_powheg-pythia8/*/*/*icenet.root/Events",   path_bdt_output_base.Data(), modeltag.Data(), path_bdt_2024.Data()) );
   background1_bdt->Add(  TString::Format("%s/%s/%s/QCD_Bin-PT-15to20_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8/*/*/*icenet.root/Events",   path_bdt_output_base.Data(), modeltag.Data(), path_bdt_2024.Data()) );
   background2_bdt->Add(  TString::Format("%s/%s/%s/QCD_Bin-PT-20to30_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8/*/*/*icenet.root/Events",   path_bdt_output_base.Data(), modeltag.Data(), path_bdt_2024.Data()) );
   background3_bdt->Add(  TString::Format("%s/%s/%s/QCD_Bin-PT-30to50_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8/*/*/*icenet.root/Events",   path_bdt_output_base.Data(), modeltag.Data(), path_bdt_2024.Data()) );
   background4_bdt->Add(  TString::Format("%s/%s/%s/QCD_Bin-PT-50to80_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8/*/*/*icenet.root/Events",   path_bdt_output_base.Data(), modeltag.Data(), path_bdt_2024.Data()) );
   background5_bdt->Add(  TString::Format("%s/%s/%s/QCD_Bin-PT-80to120_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8/*/*/*icenet.root/Events",  path_bdt_output_base.Data(), modeltag.Data(), path_bdt_2024.Data()) );
   background6_bdt->Add(  TString::Format("%s/%s/%s/QCD_Bin-PT-120to170_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8/*/*/*icenet.root/Events", path_bdt_output_base.Data(), modeltag.Data(), path_bdt_2024.Data()) );
   background7_bdt->Add(  TString::Format("%s/%s/%s/QCD_Bin-PT-170to300_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8/*/*/*icenet.root/Events", path_bdt_output_base.Data(), modeltag.Data(), path_bdt_2024.Data()) );
   background8_bdt->Add(  TString::Format("%s/%s/%s/QCD_Bin-PT-300to470_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8/*/*/*icenet.root/Events", path_bdt_output_base.Data(), modeltag.Data(), path_bdt_2024.Data()) );
   background9_bdt->Add(  TString::Format("%s/%s/%s/QCD_Bin-PT-470to600_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8/*/*/*icenet.root/Events", path_bdt_output_base.Data(), modeltag.Data(), path_bdt_2024.Data()) );
   background10_bdt->Add( TString::Format("%s/%s/%s/QCD_Bin-PT-600to800_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8/*/*/*icenet.root/Events", path_bdt_output_base.Data(), modeltag.Data(), path_bdt_2024.Data()) );
   background11_bdt->Add( TString::Format("%s/%s/%s/QCD_Bin-PT-800to1000_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8/*/*/*icenet.root/Events",path_bdt_output_base.Data(), modeltag.Data(), path_bdt_2024.Data()) );
   background12_bdt->Add( TString::Format("%s/%s/%s/QCD_Bin-PT-1000_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8/*/*/*icenet.root/Events",     path_bdt_output_base.Data(), modeltag.Data(), path_bdt_2024.Data()) );


   /*
   ROOT::RDataFrame df(*signal);
   ROOT::RDataFrame df1(*background1);
   ROOT::RDataFrame df2(*background2);
   ROOT::RDataFrame df3(*background3);
   ROOT::RDataFrame df4(*background4);
   ROOT::RDataFrame df5(*background5);
   ROOT::RDataFrame df6(*background6);
   ROOT::RDataFrame df7(*background7);
   ROOT::RDataFrame df8(*background8);
   ROOT::RDataFrame df9(*background9);
   ROOT::RDataFrame df10(*background10);
   ROOT::RDataFrame df11(*background11);
   ROOT::RDataFrame df12(*background12);
   */

   ROOT::RDataFrame df_bdt(*signal_bdt);
   ROOT::RDataFrame df_bdt2(*signal_bdt2);
   ROOT::RDataFrame df1_bdt(*background1_bdt);
   ROOT::RDataFrame df2_bdt(*background2_bdt);
   ROOT::RDataFrame df3_bdt(*background3_bdt);
   ROOT::RDataFrame df4_bdt(*background4_bdt);
   ROOT::RDataFrame df5_bdt(*background5_bdt);
   ROOT::RDataFrame df6_bdt(*background6_bdt);
   ROOT::RDataFrame df7_bdt(*background7_bdt);
   ROOT::RDataFrame df8_bdt(*background8_bdt);
   ROOT::RDataFrame df9_bdt(*background9_bdt);
   ROOT::RDataFrame df10_bdt(*background10_bdt);
   ROOT::RDataFrame df11_bdt(*background11_bdt);
   ROOT::RDataFrame df12_bdt(*background12_bdt);


   std::cout << "RDataFrame done" << std::endl;
   /*
   auto entries = df.Count();
   auto entries1 = df1.Count();
   auto entries2 = df2.Count();
   auto entries3 = df3.Count();
   auto entries4 = df4.Count();
   auto entries5 = df5.Count();
   auto entries6 = df6.Count();
   auto entries7 = df7.Count();
   auto entries8 = df8.Count();
   auto entries9 = df9.Count();
   auto entries10 = df10.Count();
   auto entries11 = df11.Count();
   auto entries12 = df12.Count();
   */

   auto entries_bdt = df_bdt.Count();
   auto entries_bdt2 = df_bdt2.Count();
   auto entries1_bdt = df1_bdt.Count();
   auto entries2_bdt = df2_bdt.Count();
   auto entries3_bdt = df3_bdt.Count();
   auto entries4_bdt = df4_bdt.Count();
   auto entries5_bdt = df5_bdt.Count();
   auto entries6_bdt = df6_bdt.Count();
   auto entries7_bdt = df7_bdt.Count();
   auto entries8_bdt = df8_bdt.Count();
   auto entries9_bdt = df9_bdt.Count();
   auto entries10_bdt = df10_bdt.Count();
   auto entries11_bdt = df11_bdt.Count();
   auto entries12_bdt = df12_bdt.Count();
  
   /*
   Double_t no_of_entries = *entries;
   Double_t no_of_entries1 = *entries1;
   Double_t no_of_entries2 = *entries2;
   Double_t no_of_entries3 = *entries3;
   Double_t no_of_entries4 = *entries4;
   Double_t no_of_entries5 = *entries5;
   Double_t no_of_entries6 = *entries6;
   Double_t no_of_entries7 = *entries7;
   Double_t no_of_entries8 = *entries8;
   Double_t no_of_entries9 = *entries9;
   Double_t no_of_entries10 = *entries10;
   Double_t no_of_entries11 = *entries11;
   Double_t no_of_entries12 = *entries12;
   */
   
   Double_t no_of_entries_original = *entries_bdt;
   Double_t no_of_entries_original2 = *entries_bdt2;
   Double_t no_of_entries1_original = *entries1_bdt; 
   Double_t no_of_entries2_original = *entries2_bdt;
   Double_t no_of_entries3_original = *entries3_bdt;
   Double_t no_of_entries4_original = *entries4_bdt;
   Double_t no_of_entries5_original = *entries5_bdt;
   Double_t no_of_entries6_original = *entries6_bdt;
   Double_t no_of_entries7_original = *entries7_bdt;
   Double_t no_of_entries8_original = *entries8_bdt;
   Double_t no_of_entries9_original = *entries9_bdt;
   Double_t no_of_entries10_original = *entries10_bdt; 
   Double_t no_of_entries11_original = *entries11_bdt;
   Double_t no_of_entries12_original = *entries12_bdt;

   
   /*
   Double_t no_of_entries1_original = 4576065;
   Double_t no_of_entries2_original = 30612338;
   Double_t no_of_entries3_original = 29884616;
   Double_t no_of_entries4_original = 20116013;
   Double_t no_of_entries5_original = 612919;
   Double_t no_of_entries6_original = 584368;
   Double_t no_of_entries7_original = 35187520;
   Double_t no_of_entries8_original = 492418;
   Double_t no_of_entries9_original = 492716;
   Double_t no_of_entries10_original = 16618977;
   Double_t no_of_entries11_original = 16749914;
   Double_t no_of_entries12_original = 10719790;
   */

   /*  
   std::cout << *entries << "entries in signal" << std::endl;
   std::cout << *entries1 << " entries in background1" << std::endl;
   std::cout << *entries2 << " entries in background2" << std::endl;
   std::cout << *entries3 << " entries in background3" << std::endl;
   std::cout << *entries4 << " entries in background4" << std::endl;
   std::cout << *entries5 << " entries in background5" << std::endl;
   std::cout << *entries6 << " entries in background6" << std::endl;
   std::cout << *entries7 << " entries in background7" << std::endl;
   std::cout << *entries8 << " entries in background8" << std::endl;
   std::cout << *entries9 << " entries in background9" << std::endl;
   std::cout << *entries10 << " entries in background10" << std::endl;
   std::cout << *entries11 << " entries in background11" << std::endl;
   std::cout << *entries12 << " entries in background12" << std::endl;
   */

   std::cout << no_of_entries_original << "entries in original signal" << std::endl;
   std::cout << no_of_entries1_original << " entries in original background1" << std::endl;
   std::cout << no_of_entries2_original << " entries in original background2" << std::endl;
   std::cout << no_of_entries3_original << " entries in original background3" << std::endl;
   std::cout << no_of_entries4_original << " entries in original background4" << std::endl;
   std::cout << no_of_entries5_original << " entries in original background5" << std::endl;
   std::cout << no_of_entries6_original << " entries in original background6" << std::endl;
   std::cout << no_of_entries7_original << " entries in original background7" << std::endl;
   std::cout << no_of_entries8_original << " entries in original background8" << std::endl;
   std::cout << no_of_entries9_original << " entries in original background9" << std::endl;
   std::cout << no_of_entries10_original << " entries in original background10" << std::endl;
   std::cout << no_of_entries11_original << " entries in original background11" << std::endl;
   std::cout << no_of_entries12_original << " entries in original background12" << std::endl;

   // When runTest=kTRUE: limit each dataframe to 1% of its events.
   // When runTest=kFALSE: pass the full dataframe through unchanged.
   // Scaling uses no_of_entries_*_original (full counts) in both cases.
   ROOT::RDF::RNode node_bdt   = runTest ? (ROOT::RDF::RNode)df_bdt.Range((ULong64_t)(no_of_entries_original   * 0.01)) : (ROOT::RDF::RNode)df_bdt;
   ROOT::RDF::RNode node_bdt2  = runTest ? (ROOT::RDF::RNode)df_bdt2.Range((ULong64_t)(no_of_entries_original2  * 0.01)) : (ROOT::RDF::RNode)df_bdt2;
   ROOT::RDF::RNode node1_bdt  = runTest ? (ROOT::RDF::RNode)df1_bdt.Range((ULong64_t)(no_of_entries1_original  * 0.01)) : (ROOT::RDF::RNode)df1_bdt;
   ROOT::RDF::RNode node2_bdt  = runTest ? (ROOT::RDF::RNode)df2_bdt.Range((ULong64_t)(no_of_entries2_original  * 0.01)) : (ROOT::RDF::RNode)df2_bdt;
   ROOT::RDF::RNode node3_bdt  = runTest ? (ROOT::RDF::RNode)df3_bdt.Range((ULong64_t)(no_of_entries3_original  * 0.01)) : (ROOT::RDF::RNode)df3_bdt;
   ROOT::RDF::RNode node4_bdt  = runTest ? (ROOT::RDF::RNode)df4_bdt.Range((ULong64_t)(no_of_entries4_original  * 0.01)) : (ROOT::RDF::RNode)df4_bdt;
   ROOT::RDF::RNode node5_bdt  = runTest ? (ROOT::RDF::RNode)df5_bdt.Range((ULong64_t)(no_of_entries5_original  * 0.01)) : (ROOT::RDF::RNode)df5_bdt;
   ROOT::RDF::RNode node6_bdt  = runTest ? (ROOT::RDF::RNode)df6_bdt.Range((ULong64_t)(no_of_entries6_original  * 0.01)) : (ROOT::RDF::RNode)df6_bdt;
   ROOT::RDF::RNode node7_bdt  = runTest ? (ROOT::RDF::RNode)df7_bdt.Range((ULong64_t)(no_of_entries7_original  * 0.01)) : (ROOT::RDF::RNode)df7_bdt;
   ROOT::RDF::RNode node8_bdt  = runTest ? (ROOT::RDF::RNode)df8_bdt.Range((ULong64_t)(no_of_entries8_original  * 0.01)) : (ROOT::RDF::RNode)df8_bdt;
   ROOT::RDF::RNode node9_bdt  = runTest ? (ROOT::RDF::RNode)df9_bdt.Range((ULong64_t)(no_of_entries9_original  * 0.01)) : (ROOT::RDF::RNode)df9_bdt;
   ROOT::RDF::RNode node10_bdt = runTest ? (ROOT::RDF::RNode)df10_bdt.Range((ULong64_t)(no_of_entries10_original * 0.01)) : (ROOT::RDF::RNode)df10_bdt;
   ROOT::RDF::RNode node11_bdt = runTest ? (ROOT::RDF::RNode)df11_bdt.Range((ULong64_t)(no_of_entries11_original * 0.01)) : (ROOT::RDF::RNode)df11_bdt;
   ROOT::RDF::RNode node12_bdt = runTest ? (ROOT::RDF::RNode)df12_bdt.Range((ULong64_t)(no_of_entries12_original * 0.01)) : (ROOT::RDF::RNode)df12_bdt;

   // Number of bins. More bins = better threshold estimation, but worse plot drawing
   Int_t nbins_histo = 10000;
   //Int_t nbins_histo = 1000;

   auto h_bdt = node_bdt.Histo1D({"hist_bdt_nano", "BDT score distribution; BDT score; Number of events", nbins_histo, 0.0, 1.0}, BDTmodel.Data());
   auto h_bdt_signal2 = node_bdt2.Histo1D({"hist_bdt_nano", "BDT score distribution; BDT score; Number of events", nbins_histo, 0.0, 1.0}, BDTmodel.Data());
   auto h_bdt1 = node1_bdt.Histo1D({"hist_bdt1_nano", "BDT cut efficiency; BDT score; Efficiency", nbins_histo, 0.0, 1.0}, BDTmodel.Data());
   auto h_bdt2 = node2_bdt.Histo1D({"hist_bdt2_nano", "BDT score distribution; BDT score; Number of events", nbins_histo, 0.0, 1.0}, BDTmodel.Data());
   auto h_bdt3 = node3_bdt.Histo1D({"hist_bdt3_nano", "BDT score distribution; BDT score; Number of events", nbins_histo, 0.0, 1.0}, BDTmodel.Data());
   auto h_bdt4 = node4_bdt.Histo1D({"hist_bdt4_nano", "BDT score distribution; BDT score; Number of events", nbins_histo, 0.0, 1.0}, BDTmodel.Data());
   auto h_bdt5 = node5_bdt.Histo1D({"hist_bdt5_nano", "BDT score distribution; BDT score; Number of events", nbins_histo, 0.0, 1.0}, BDTmodel.Data());
   auto h_bdt6 = node6_bdt.Histo1D({"hist_bdt6_nano", "BDT score distribution; BDT score; Number of events", nbins_histo, 0.0, 1.0}, BDTmodel.Data());
   auto h_bdt7 = node7_bdt.Histo1D({"hist_bdt7_nano", "BDT score distribution; BDT score; Number of events", nbins_histo, 0.0, 1.0}, BDTmodel.Data());
   auto h_bdt8 = node8_bdt.Histo1D({"hist_bdt8_nano", "BDT score distribution; BDT score; Number of events", nbins_histo, 0.0, 1.0}, BDTmodel.Data());
   auto h_bdt9 = node9_bdt.Histo1D({"hist_bdt9_nano", "BDT score distribution; BDT score; Number of events", nbins_histo, 0.0, 1.0}, BDTmodel.Data());
   auto h_bdt10 = node10_bdt.Histo1D({"hist_bdt10_nano", "BDT score distribution; BDT score; Number of events", nbins_histo, 0.0, 1.0}, BDTmodel.Data());
   auto h_bdt11 = node11_bdt.Histo1D({"hist_bdt11_nano", "BDT score distribution; BDT score; Number of events", nbins_histo, 0.0, 1.0}, BDTmodel.Data());
   auto h_bdt12 = node12_bdt.Histo1D({"hist_bdt12_nano", "BDT score distribution; BDT score; Number of events", nbins_histo, 0.0, 1.0}, BDTmodel.Data());

   Double_t branching_ratio = 0.01;
   Double_t xs = 43.9;  

   Double_t xs1 = 2799000.0;
   Double_t xs2 = 2526000.0;
   Double_t xs3 = 1362000.0;
   Double_t xs4 = 376600.0;
   Double_t xs5 = 88930.0;
   Double_t xs6 = 21230.0; 

   Double_t xs7 = 7055.0;
   Double_t xs8 = 619.3;
   Double_t xs9 = 59.24;
   Double_t xs10 = 18.21;
   Double_t xs11 = 3.275;
   Double_t xs12 = 1.078; 

   Double_t lumi = 33.6*1000;    

   h_bdt->Scale(lumi*xs*branching_ratio/no_of_entries_original);
   h_bdt_signal2->Scale(lumi*xs*branching_ratio/no_of_entries_original2);
   h_bdt1->Scale(lumi*xs1/no_of_entries1_original);
   h_bdt2->Scale(lumi*xs2/no_of_entries2_original);
   h_bdt3->Scale(lumi*xs3/no_of_entries3_original);
   h_bdt4->Scale(lumi*xs4/no_of_entries4_original);
   h_bdt5->Scale(lumi*xs5/no_of_entries5_original);
   h_bdt6->Scale(lumi*xs6/no_of_entries6_original);
   h_bdt7->Scale(lumi*xs7/no_of_entries7_original);
   h_bdt8->Scale(lumi*xs8/no_of_entries8_original);
   h_bdt9->Scale(lumi*xs9/no_of_entries9_original);
   h_bdt10->Scale(lumi*xs10/no_of_entries10_original);
   h_bdt11->Scale(lumi*xs11/no_of_entries11_original);
   h_bdt12->Scale(lumi*xs12/no_of_entries12_original);


   gStyle->SetOptStat(0); gStyle->SetTextFont(42);

   h_bdt1->Add(h_bdt2.GetPtr());
   h_bdt1->Add(h_bdt3.GetPtr());
   h_bdt1->Add(h_bdt4.GetPtr());
   h_bdt1->Add(h_bdt5.GetPtr());
   h_bdt1->Add(h_bdt6.GetPtr());
   h_bdt1->Add(h_bdt7.GetPtr());
   h_bdt1->Add(h_bdt8.GetPtr());
   h_bdt1->Add(h_bdt9.GetPtr());
   h_bdt1->Add(h_bdt10.GetPtr());
   h_bdt1->Add(h_bdt11.GetPtr());
   h_bdt1->Add(h_bdt12.GetPtr());

   
   TH1D* hist_bdt_expected = (TH1D*)h_bdt->Clone("hist_bdt_expected");
   TH1D* hist_bdt_signal2_expected = (TH1D*)h_bdt_signal2->Clone("hist_bdt_signal2_expected");
   TH1D* hist_bdt1_expected = (TH1D*)h_bdt1->Clone("hist_bdt1_expected");

   hist_bdt_expected->Scale(1.0/hist_bdt_expected->Integral());
   hist_bdt_signal2_expected->Scale(1.0/hist_bdt_signal2_expected->Integral());
   hist_bdt1_expected->Scale(1.0/hist_bdt1_expected->Integral());

   TCanvas* c1 = new TCanvas("", "", 800, 700);
   c1->SetLogy();
   hist_bdt1_expected->SetLineColor(kRed);
   hist_bdt1_expected->GetYaxis()->SetTitleOffset(1.40);
   hist_bdt_expected->SetLineColor(kBlue+1);
   hist_bdt_signal2_expected->SetLineColor(kMagenta);
   hist_bdt1_expected->SetTitle("");
   hist_bdt1_expected->GetYaxis()->SetTitle("Fraction of events");

   TH1D* hist_bdt1_expected_rebinned = (TH1D*)hist_bdt1_expected->Clone("hist_bdt1_expected_rebinned");
   TH1D* hist_bdt_expected_rebinned = (TH1D*)hist_bdt_expected->Clone("hist_bdt_expected_rebinned");
   TH1D* hist_bdt_signal2_expected_rebinned = (TH1D*)hist_bdt_signal2_expected->Clone("hist_bdt_signal2_expected_rebinned");
   hist_bdt1_expected_rebinned->Rebin(10);
   hist_bdt_expected_rebinned->Rebin(10);
   hist_bdt_signal2_expected_rebinned->Rebin(10);

   //hist_bdt1_expected->DrawClone("hist");
   //hist_bdt_expected->DrawClone("hist Same");
   //hist_bdt_signal2_expected->DrawClone("hist Same");
   hist_bdt1_expected_rebinned->DrawClone("hist");
   hist_bdt_expected_rebinned->DrawClone("hist Same");
   hist_bdt_signal2_expected_rebinned->DrawClone("hist Same");
   
   TLegend* out_legend1 = new TLegend(0.2, 0.12, 0.88, 0.22);
   out_legend1->AddEntry(hist_bdt1_expected->GetName(), "Background", "l");
   out_legend1->AddEntry(hist_bdt_expected->GetName(), "Scenario A, m_{#pi} = 4 GeV, m_{A}=1.33 GeV, c #tau = 1 mm", "l");
   out_legend1->AddEntry(hist_bdt_signal2_expected->GetName(), "Scenario A, m_{#pi} = 4 GeV, m_{A}=1.33 GeV, c #tau = 10 mm", "l");
   out_legend1->Draw("Same"); 



   TH1D* hist_bdt_c = (TH1D*)hist_bdt_expected ->GetCumulative(kFALSE);
   TH1D* hist_bdt_signal2_c = (TH1D*)hist_bdt_signal2_expected ->GetCumulative(kFALSE);
   TH1D* hist_bdt1_c = (TH1D*)hist_bdt1_expected ->GetCumulative(kFALSE);
   
   TH1D* ssqrtb_bdt = make_hist("ssqrtb_bdtscore", "BDT score distribution", 10000, 0., 1.0, "BDT score", "#splitline{Median}{discovery significance}", kBlack, 2);

   for(Int_t i = 1; i < 10001; i++){
       ssqrtb_bdt->SetBinContent(i, ssqrtb(hist_bdt_c->GetBinContent(i), hist_bdt1_c->GetBinContent(i)));
       ssqrtb_bdt->SetBinError(i, ssqrtb_error(hist_bdt_c->GetBinContent(i), hist_bdt1_c->GetBinContent(i)));
   }
   TH1D* Z_bdt = make_hist("Z_bdt", "BDT score distribution", 10000, 0., 1.0, "BDT score", "#splitline{Median}{discovery significance}", kBlack, 2);

   for(Int_t i = 1; i < 10001; i++){
       Z_bdt->SetBinContent(i, Z(hist_bdt_c->GetBinContent(i), hist_bdt1_c->GetBinContent(i)));
       Z_bdt->SetBinError(i, Z_error(hist_bdt_c->GetBinContent(i), hist_bdt1_c->GetBinContent(i)));
   } 


   TCanvas* c2 = new TCanvas("", "", 800, 700);
   c2->SetLogy();
   hist_bdt1_c->SetLineColor(kRed);
   hist_bdt1_c->DrawClone("hist");
  
   TH1D* hist_bdt1_c_inverted = make_hist("hist_bdt1_c_inverted", "BDT cut efficiency", nbins_histo, 0.0, 1.0, "1 - BDT score", "Efficiency", kRed, 2);
   for(Int_t i = 1; i <= nbins_histo; i++){
       hist_bdt1_c_inverted->SetBinContent(i, hist_bdt1_c->GetBinContent(nbins_histo+1-i));
   }
   
   TH1D* hist_bdt_c_inverted = make_hist("hist_bdt_c_inverted", "BDT cut efficiency", nbins_histo, 0.0, 1.0, "1 - BDT score", "Efficiency", kRed, 2);
   for(Int_t i = 1; i <= nbins_histo; i++){
       hist_bdt_c_inverted->SetBinContent(i, hist_bdt_c->GetBinContent(nbins_histo+1-i));
   }

   TH1D* hist_bdt_signal2_c_inverted = make_hist("hist_bdt_signal2_c_inverted", "BDT cut efficiency", nbins_histo, 0.0, 1.0, "1 - BDT score", "Efficiency", kRed, 2);
   for(Int_t i = 1; i <= nbins_histo; i++){
       hist_bdt_signal2_c_inverted->SetBinContent(i, hist_bdt_signal2_c->GetBinContent(nbins_histo+1-i));
   }
   
   TCanvas* c3 = new TCanvas("", "", 800, 700);
   c3->SetLogy();
   c3->SetLogx();
   hist_bdt1_c_inverted->SetLineColor(kRed);
   hist_bdt1_c_inverted->GetYaxis()->SetTitleOffset(1.35);
   //hist_bdt_c_inverted->SetLineColor(kBlue+1);
   hist_bdt1_c_inverted->DrawClone("hist");
   //hist_bdt_c_inverted->DrawClone("hist Same");  


   TCanvas* c4 = new TCanvas("", "", 800, 700);
   c4->SetLogy();
   c4->SetLogx();
   hist_bdt1_c_inverted->SetLineColor(kRed);
   hist_bdt_c_inverted->SetLineColor(kBlue+1);
   hist_bdt_signal2_c_inverted->SetLineColor(kMagenta);

   hist_bdt1_c_inverted->GetXaxis()->SetRangeUser(1E-4,1.0);
   hist_bdt1_c_inverted->SetTitle("");

   hist_bdt1_c_inverted->DrawClone("hist");
   hist_bdt_c_inverted->DrawClone("hist Same");
   hist_bdt_signal2_c_inverted->DrawClone("hist Same");
   TLegend* out_legend2 = new TLegend(0.78, 0.695, 0.98, 0.775);
   out_legend2->AddEntry(hist_bdt1_c_inverted->GetName(), "Background", "l");
   out_legend2->AddEntry(hist_bdt_c_inverted->GetName(), "Scenario A, m_{#pi} = 4 GeV, m_{A}=1.33 GeV, c #tau = 1 mm", "l");
   out_legend2->AddEntry(hist_bdt_signal2_c_inverted->GetName(), "Scenario A, m_{#pi} = 4 GeV, m_{A}=1.33 GeV, c #tau = 10 mm", "l");
   //out_legend2->Draw("Same");    


   //hist_draw_3({ hist_bdt_expected, hist_bdt1_expected, Z_bdt, ssqrtb_bdt });
   //hist_stack("h_svmass_nano_stack", h_svmass1, h_svmass, Ntuple_name, true);
  // hist_stack_old("h_bdt_stack", { hist_bdt1_c, hist_bdt_c });

   // Find BDT threshold for background rejection = 1e-4
   // hist_bdt1_c  = reverse cumulative of normalised background (eff for score >= threshold)
   // hist_bdt_c   = reverse cumulative of normalised signal 1 (ctau=1p0)
   // hist_bdt_signal2_c = reverse cumulative of normalised signal 2 (ctau=10)
   Double_t bdt_threshold = -1.0, bkg_eff_at_thr = -1.0;
   Double_t sig1_eff_at_thr = -1.0, sig2_eff_at_thr = -1.0;
   for(Int_t i = 1; i <= 10000; i++){
      if(hist_bdt1_c->GetBinContent(i) <= 1e-4){
         bdt_threshold    = hist_bdt1_c->GetBinLowEdge(i);
         bkg_eff_at_thr   = hist_bdt1_c->GetBinContent(i);
         sig1_eff_at_thr  = hist_bdt_c->GetBinContent(i);
         sig2_eff_at_thr  = hist_bdt_signal2_c->GetBinContent(i);
         printf("BDT threshold for bkg rejection 1e-4: score > %.4f  (bkg eff = %.2e)\n",
                bdt_threshold, bkg_eff_at_thr);
         printf("  Signal eff (ctau=1p0, mA=1p33, mpi=4): %.4f\n", sig1_eff_at_thr);
         printf("  Signal eff (ctau=10,  mA=1p33, mpi=4): %.4f\n", sig2_eff_at_thr);
         Double_t sig_eff_mean = (sig1_eff_at_thr + sig2_eff_at_thr)*0.5;
         printf("  Signal eff: %.4f\n", sig_eff_mean);
         break;
      }
   }

   // Redraw c4 legend with threshold and efficiency values
   c4->cd();
   TLegend* out_legend2b = new TLegend(0.35, 0.12, 0.90, 0.30);
   out_legend2b->SetFillColor(0);
   out_legend2b->SetBorderSize(0);
   out_legend2b->SetTextSize(0.025);
   out_legend2b->AddEntry((TObject*)0,
      TString::Format("BDT cut: score > %.4f", bdt_threshold), "");
   out_legend2b->AddEntry(hist_bdt1_c_inverted->GetName(),
      TString::Format("Background  (#varepsilon = %.2e)", bkg_eff_at_thr), "l");
   out_legend2b->AddEntry(hist_bdt_c_inverted->GetName(),
      TString::Format("scA, m_{#pi} = 4 GeV, m_{A}=1.33 GeV, c #tau = 1 mm  ( #varepsilon = %.3f)", sig1_eff_at_thr), "l");
   out_legend2b->AddEntry(hist_bdt_signal2_c_inverted->GetName(),
      TString::Format("scA, m_{#pi} = 4 GeV, m_{A}=1.33 GeV, c #tau = 10 mm  ( #varepsilon = %.3f)", sig2_eff_at_thr), "l");
   out_legend2b->Draw("Same");
   c4->Update();

   // Save all canvases to a single multi-page PDF named after the modeltag
   TString script_dir = gSystem->DirName(__FILE__);
   TString pdf_file   = script_dir + "/" + modeltag + ".pdf";
   c1->Print((pdf_file + "(").Data());
   c2->Print(pdf_file.Data());
   c3->Print(pdf_file.Data());
   c4->Print((pdf_file + ")").Data());
   printf("Plots saved to %s\n", pdf_file.Data());
}


int getBDTScoreThreshold_fromDeploymentOutput(){
   nano_analysis();
   return 0;
}
