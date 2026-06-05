// Signal efficiency study for 2024 noMET BDT models.
//
// For each of three noMET trigger strategies (Mu10ORDoubleMu, Mu10, DoubleMu):
//   1. Finds the xgb01_NOJETS_NOMET BDT score threshold that achieves
//      1e-4 background (QCD) fake-rate after trigger+kinematic pre-selection.
//   2. Evaluates trigger+kinematic and BDT cut efficiencies for 4 signal points
//      (mpi=4 GeV, mA=1.33 GeV, ctau in {0.1, 1, 10, 100} mm)
//      and prints a formatted cutflow table.
//   [Groups 1-4 (mA=0.40, mpi=10/mA=1.00, mpi=1/mA=0.33, mpi=10/mA=3.33)
//    are temporarily disabled via #if 0 blocks -- flip to #if 1 to re-enable.]
//   3. Saves efficiency-curve plots to a multi-page PDF.
//
// Deployment output convention (icenet dqcd_deploy.py):
//   Events failing the trigger+kinematic pre-selection are assigned score = -1
//   (via aux.unmask default_value=-1).  Valid BDT scores are in [0, 1].
//     score >= 0  -> passed trigger + kinematic pre-selection
//     score >  T  -> also passes BDT cut at threshold T
//   The [0,1] histogram range automatically excludes score=-1 events from the
//   bins, so all efficiency curves are relative to the trigger-selected sample.
//
// Kinematic cuts encoded per deployment (configs/dqcd/filter.py, cuts.py):
//   noMET (OR):     HLT_Mu10_Barrel_L1HP11_IP6 OR HLT_DoubleMu4_3_LowMass
//                   Mu10 leg:     (mu1pt>10 && |mu1eta|<1.2) || (mu2pt>10 && |mu2eta|<1.2)
//                   DoubleMu leg: max(mu1pt,mu2pt)>4 && min(mu1pt,mu2pt)>3
//   Mu10_noMET:     HLT_Mu10_Barrel_L1HP11_IP6 + Mu10 kinematic leg only
//   DoubleMu_noMET: HLT_DoubleMu4_3_LowMass    + DoubleMu kinematic leg only
//
// lxy displacement categories (mpi=4 GeV, mean lxy ~ ctau x <pT_pi/mpi> ~ ctau x 4-6):
//   LOW   ctau=0.1mm  -> mean lxy ~0.05 cm  >99.9% of dimuons below 1 cm
//   MED   ctau=1mm    -> mean lxy ~0.5 cm   transitional (most below 1 cm)
//   HIGH  ctau=10mm   -> mean lxy ~5 cm     majority above 1 cm
//   XHIGH ctau=100mm  -> mean lxy ~50 cm    extreme (near tracker edge)
//
// Output directory: _tools/output_getSignalEff_noMET/
//
// Efficiency-curve PDFs (one per HLT model, "tables" keyword not required):
//   getSignalEff_noMET_Mu10ORDoubleMu.pdf
//   getSignalEff_noMET_Mu10.pdf
//   getSignalEff_noMET_DoubleMu.pdf
// Each PDF:
//   Page 1 - BDT efficiency overview for that model
//   Page 2 - muonSV mass distributions (3 lxy cols x 4 ctau rows)
// NB: mass pages read NanoAOD via xrootd -> run `source setproxy.sh` first.
//
// Tables PDFs + CSVs ("tables" keyword required):
//   tables_Mu10ORDoubleMu.pdf / .csv
//   tables_Mu10.pdf           / .csv
//   tables_DoubleMu.pdf       / .csv
// Each PDF (compiled from LaTeX via pdflatex) has two pages per signal SET
// (= fixed mpi,mA with its 4 ctau values as columns), then one shared QCD page:
//   per signal set:
//     cutflow page - 4 ctau columns: absolute counts + relative eff w.r.t. N_total,
//                    plus a QCD Total (weighted) / Weighted-eps comparison column
//     lxy page     - 4 ctau x 3 lxy bins (landscape): counts + eff rel. to N in bin
//   QCD cutflow page - 12 pT bins + xs-weighted total/eff column (landscape)
//   => 5 sets give 5x2 + 1 = 11 pages; "setN" restricts to one set (3 pages).
//
// Usage:
//   root -l -b -q getSignalEff_noMET.C                   # full run (all pages)
//   root -l -b -q 'getSignalEff_noMET.C("test")'         # fast debug (1 file each)
//   root -l -b -q 'getSignalEff_noMET.C("skip_mass")'    # skip mass pages
//   root -l -b -q 'getSignalEff_noMET.C("compute_thr")'  # recompute BDT thresholds
//   root -l -b -q 'getSignalEff_noMET.C("model1")'       # run only one model (0/1/2)
//   root -l -b -q 'getSignalEff_noMET.C("set2")'         # run only one signal set (0-4)
//
// compute_thr weights each QCD event by its bin cross section sigma_q (no /N, no
// luminosity), so the FPR=1e-4 threshold is set on the cross-section-weighted QCD
// event count Sigma_q sigma_q N_q / Sigma_q sigma_q -- the same weighting as the
// "Total (weighted)" QCD column in the tables.

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TCanvas.h"
#include "TChain.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include "TChainElement.h"
#include <iostream>
#include <cmath>

using namespace ROOT::VecOps;

// --- constants ----------------------------------------------------------------

static Bool_t         RUN_TEST     = kFALSE; // set via argument: .C("test")
static Bool_t         COMPUTE_THR  = kFALSE; // set via argument: .C("compute_thr")
static Bool_t         SKIP_MASS    = kFALSE; // set via argument: .C("skip_mass") to omit mass pages
static Bool_t         MAKE_TABLES  = kFALSE; // set via argument: .C("tables") to produce tables PDF
static Int_t          ONLY_MODEL   = -1;     // set via "model0"/"model1"/"model2"; -1 = all 3 models
static Int_t          ONLY_SET     = -1;     // set via "set0".."set4"; -1 = all 5 signal mass sets
static const Long64_t TEST_EVTS    = 5000;   // events per sample in test mode

// True if model index m should be processed given the ONLY_MODEL selector.
static inline Bool_t run_model(Int_t m) { return ONLY_MODEL < 0 || m == ONLY_MODEL; }
// True if signal mass set g (group of 4 ctau) should be processed.
static inline Bool_t run_set(Int_t g) { return ONLY_SET < 0 || g == ONLY_SET; }
// True if signal point index s should be processed (its set = s/4 is selected).
static inline Bool_t run_sig(Int_t s) { return ONLY_SET < 0 || s / 4 == ONLY_SET; }

// NanoAOD access for mass plots (requires setproxy.sh before running).
// The deployment output mirrors the grid path locally; derive xrootd URLs from it.
static const char* XROOTD_HOST = "gfe02.grid.hep.ph.ic.ac.uk";

// lxy (= muonSV_dxy) category boundaries [cm]
static const Float_t  LXY_LO[3]  = {  0.f,  1.f, 10.f };
static const Float_t  LXY_HI[3]  = {  1.f, 10.f, 100.f };
static const char*    LXY_LABEL[3] = {
    "l_{xy}#in[0,1] cm", "l_{xy}#in[1,10] cm", "l_{xy}#in[10,100] cm" };

// Hardcoded BDT thresholds at FPR = 1e-4 (from full QCD background study).
// Order matches MODELS[]: [0] Mu10ORDoubleMu, [1] Mu10, [2] DoubleMu.
static const Double_t HARDCODED_THR[3] = { 0.9971, 0.9969, 0.9978 };

static const char*    BRANCH   = "xgb01_NOJETS_NOMET";
static const Int_t    NBINS    = 10000;
static const Double_t TARGET   = 1e-4;   // target background fake rate
static const Double_t LUMI     = 109.95e3; // pb^-1  (2024 total)
static const Double_t SIG_XS   = 43.9;  // pb  (ggH production)
static const Double_t SIG_BR   = 0.01;  // assumed signal branching ratio

static const char* DEPLOY_BASE =
    "/home/hep/jtafoyav/vols/parking/bdt/icenet/output/dqcd/deploy";
static const char* GRID_PATH =
    "gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/"
    "user/tafoyava/samples/bParking/2024";

// --- configuration structs ----------------------------------------------------

struct ModelCfg { const char* name; const char* modeltag; };
static const ModelCfg MODELS[3] = {
    { "Mu10ORDoubleMu",
      "modeltag__scenarioA_all_no_DA_old_BDT_with_dRmSV_2024_noMET" },
    { "Mu10",
      "modeltag__scenarioA_all_no_DA_old_BDT_with_dRmSV_2024_Mu10_noMET" },
    { "DoubleMu",
      "modeltag__scenarioA_all_no_DA_old_BDT_with_dRmSV_2024_DoubleMu_noMET" },
};
static const Int_t N_MODELS = 3;

// Signal points: 5 mass sets (groups) x 4 ctau values = 20 points.
// Each group of 4 consecutive entries is one signal "set" (fixed mpi, mA);
// the "setN" run argument selects a single set, and the tables PDF devotes one
// cutflow page + one lxy page to each set (4 ctau columns), exactly as before.
// lxy labels are approximate (mean lxy ~ ctau x pT_pi/mpi):
//   mpi=4:  boost ~5 -> ctau=1mm gives ~0.5cm (MED)
//   mpi=10: boost ~2 -> ctau=1mm gives ~0.2cm (LOW/MED boundary)
//   mpi=1:  boost ~20 -> ctau=0.1mm already gives ~0.2cm (MED)
struct SigCfg { const char* label; const char* dirname; const char* lxy; };
static const SigCfg SIG[20] = {
    // group 0 (set0): mpi=4, mA=1.33
    { "mpi=4, mA=1.33, ctau=0.1mm",
      "GluGluHToDarkShowers-ScenarioA_Par-ctau-0p1-mA-1p33-mpi-4_TuneCP5_13p6TeV_powheg-pythia8",
      "LOW"   },
    { "mpi=4, mA=1.33, ctau=1mm",
      "GluGluHToDarkShowers-ScenarioA_Par-ctau-1p0-mA-1p33-mpi-4_TuneCP5_13p6TeV_powheg-pythia8",
      "MED"   },
    { "mpi=4, mA=1.33, ctau=10mm",
      "GluGluHToDarkShowers-ScenarioA_Par-ctau-10-mA-1p33-mpi-4_TuneCP5_13p6TeV_powheg-pythia8",
      "HIGH"  },
    { "mpi=4, mA=1.33, ctau=100mm",
      "GluGluHToDarkShowers-ScenarioA_Par-ctau-100-mA-1p33-mpi-4_TuneCP5_13p6TeV_powheg-pythia8",
      "XHIGH" },
    // group 1 (set1): mpi=4, mA=0.40
    { "mpi=4, mA=0.40, ctau=0.1mm",
      "GluGluHToDarkShowers-ScenarioA_Par-ctau-0p1-mA-0p40-mpi-4_TuneCP5_13p6TeV_powheg-pythia8",
      "LOW"   },
    { "mpi=4, mA=0.40, ctau=1mm",
      "GluGluHToDarkShowers-ScenarioA_Par-ctau-1p0-mA-0p40-mpi-4_TuneCP5_13p6TeV_powheg-pythia8",
      "MED"   },
    { "mpi=4, mA=0.40, ctau=10mm",
      "GluGluHToDarkShowers-ScenarioA_Par-ctau-10-mA-0p40-mpi-4_TuneCP5_13p6TeV_powheg-pythia8",
      "HIGH"  },
    { "mpi=4, mA=0.40, ctau=100mm",
      "GluGluHToDarkShowers-ScenarioA_Par-ctau-100-mA-0p40-mpi-4_TuneCP5_13p6TeV_powheg-pythia8",
      "XHIGH" },
    // group 2: mpi=10, mA=1.00
    { "mpi=10, mA=1.00, ctau=0.1mm",
      "GluGluHToDarkShowers-ScenarioA_Par-ctau-0p1-mA-1p00-mpi-10_TuneCP5_13p6TeV_powheg-pythia8",
      "LOW"   },
    { "mpi=10, mA=1.00, ctau=1mm",
      "GluGluHToDarkShowers-ScenarioA_Par-ctau-1p0-mA-1p00-mpi-10_TuneCP5_13p6TeV_powheg-pythia8",
      "LOW"   },
    { "mpi=10, mA=1.00, ctau=10mm",
      "GluGluHToDarkShowers-ScenarioA_Par-ctau-10-mA-1p00-mpi-10_TuneCP5_13p6TeV_powheg-pythia8",
      "MED"   },
    { "mpi=10, mA=1.00, ctau=100mm",
      "GluGluHToDarkShowers-ScenarioA_Par-ctau-100-mA-1p00-mpi-10_TuneCP5_13p6TeV_powheg-pythia8",
      "HIGH"  },
    // group 3: mpi=1, mA=0.33
    { "mpi=1, mA=0.33, ctau=0.1mm",
      "GluGluHToDarkShowers-ScenarioA_Par-ctau-0p1-mA-0p33-mpi-1_TuneCP5_13p6TeV_powheg-pythia8",
      "MED"   },
    { "mpi=1, mA=0.33, ctau=1mm",
      "GluGluHToDarkShowers-ScenarioA_Par-ctau-1p0-mA-0p33-mpi-1_TuneCP5_13p6TeV_powheg-pythia8",
      "HIGH"  },
    { "mpi=1, mA=0.33, ctau=10mm",
      "GluGluHToDarkShowers-ScenarioA_Par-ctau-10-mA-0p33-mpi-1_TuneCP5_13p6TeV_powheg-pythia8",
      "XHIGH" },
    { "mpi=1, mA=0.33, ctau=100mm",
      "GluGluHToDarkShowers-ScenarioA_Par-ctau-100-mA-0p33-mpi-1_TuneCP5_13p6TeV_powheg-pythia8",
      "XHIGH" },
    // group 4: mpi=10, mA=3.33
    { "mpi=10, mA=3.33, ctau=0.1mm",
      "GluGluHToDarkShowers-ScenarioA_Par-ctau-0p1-mA-3p33-mpi-10_TuneCP5_13p6TeV_powheg-pythia8",
      "LOW"   },
    { "mpi=10, mA=3.33, ctau=1mm",
      "GluGluHToDarkShowers-ScenarioA_Par-ctau-1p0-mA-3p33-mpi-10_TuneCP5_13p6TeV_powheg-pythia8",
      "LOW"   },
    { "mpi=10, mA=3.33, ctau=10mm",
      "GluGluHToDarkShowers-ScenarioA_Par-ctau-10-mA-3p33-mpi-10_TuneCP5_13p6TeV_powheg-pythia8",
      "MED"   },
    { "mpi=10, mA=3.33, ctau=100mm",
      "GluGluHToDarkShowers-ScenarioA_Par-ctau-100-mA-3p33-mpi-10_TuneCP5_13p6TeV_powheg-pythia8",
      "HIGH"  },
};
static const Int_t N_SIG    = 20;          // 5 sets x 4 ctau
static const Int_t N_CTAU   = 4;           // ctau variations per set
static const Int_t N_SETS   = N_SIG / N_CTAU;  // = 5 signal mass sets

// 12 QCD background bins [name, cross section in pb]
struct QCDCfg { const char* dirname; Double_t xs; };
static const QCDCfg QCD[12] = {
    { "QCD_Bin-PT-15to20_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8",    3018000. },
    { "QCD_Bin-PT-20to30_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8",    2701000. },
    { "QCD_Bin-PT-30to50_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8",    1461000. },
    { "QCD_Bin-PT-50to80_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8",     407600. },
    { "QCD_Bin-PT-80to120_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8",    96070.  },
    { "QCD_Bin-PT-120to170_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8",   23140.  },
    { "QCD_Bin-PT-170to300_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8",    7754.  },
    { "QCD_Bin-PT-300to470_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8",     699.6 },
    { "QCD_Bin-PT-470to600_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8",      67.67},
    { "QCD_Bin-PT-600to800_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8",      21.27},
    { "QCD_Bin-PT-800to1000_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8",      3.89 },
    { "QCD_Bin-PT-1000_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8",           1.323},
};
static const Int_t N_QCD = 12;

// --- helpers ------------------------------------------------------------------

TChain* make_chain(const char* modeltag, const char* dataset, bool verbose = true)
{
    TChain* ch = new TChain("Events");
    TString pat = TString::Format("%s/%s/%s/%s/*/*/*icenet.root",
                                   DEPLOY_BASE, modeltag, GRID_PATH, dataset);
    ch->Add(pat);
    Int_t nfiles = ch->GetListOfFiles()->GetEntries();

    if (nfiles == 0) {
        if (verbose) printf("    WARNING: no files matched  %s\n", pat.Data());
        return ch;
    }

    if (RUN_TEST && nfiles > 1) {
        // In test mode keep only the first matched file to avoid globbing overhead
        // and to stop GetEntries() from opening every file.
        TString first = ((TChainElement*)ch->GetListOfFiles()->First())->GetTitle();
        delete ch;
        ch = new TChain("Events");
        ch->Add(first);
        nfiles = 1;
        if (verbose) printf("    [TEST] using 1 file: %s\n", first.Data());
    } else if (verbose) {
        printf("    chain: %d file(s)  |  %s\n", nfiles, pat.Data());
    }

    return ch;
}

// Return a range-limited RNode in test mode, or a pass-through node in full mode.
// Using RNode as a common type allows both branches to compile even though
// RDataFrame and the result of .Range() are different concrete types.
ROOT::RDF::RNode make_rnode(ROOT::RDataFrame& df)
{
    if (RUN_TEST) {
        printf("    [TEST] limiting event loop to %lld events\n", TEST_EVTS);
        return (ROOT::RDF::RNode)df.Range((ULong64_t)TEST_EVTS);
    }
    return (ROOT::RDF::RNode)df;
}

// Normalise histogram to unit area; return original integral.
Double_t normalise(TH1D* h)
{
    Double_t ig = h->Integral();
    if (ig > 0.) h->Scale(1. / ig);
    return ig;
}

// Scan reverse-cumulative histogram for first bin with content <= target.
// Sets eff_out to the content at that bin; returns the bin lower edge.
Double_t find_thr(const TH1D* cumul, Double_t target, Double_t& eff_out)
{
    for (Int_t i = 1; i <= cumul->GetNbinsX(); ++i) {
        if (cumul->GetBinContent(i) <= target) {
            eff_out = cumul->GetBinContent(i);
            return cumul->GetBinLowEdge(i);
        }
    }
    eff_out = 0.;
    return 1.;
}

// Build a pair of matched chains: NanoAOD (via xrootd) and the local deployment
// output (BDT scores).  They are filled file-by-file from the list of local
// icenet.root files so that TChain::AddFriend keeps them 1:1 in order.
// In test mode only the first file pair is added.
void fill_nano_bdt_chains(const char* modeltag, const char* dataset,
                           TChain* nano_ch, TChain* bdt_ch, bool verbose = true)
{
    TString local_dir = TString::Format("%s/%s/%s/%s",
                                         DEPLOY_BASE, modeltag, GRID_PATH, dataset);
    TChain tmp("Events");
    tmp.Add(TString::Format("%s/*/*/*icenet.root", local_dir.Data()));
    Int_t nfiles = tmp.GetListOfFiles()->GetEntries();

    if (verbose) printf("    %d deployment file(s) found\n", nfiles);
    if (nfiles == 0) { printf("    WARNING: no files found for %s\n", dataset); return; }

    TString strip_prefix = TString::Format("%s/%s/", DEPLOY_BASE, modeltag);

    Int_t limit = (RUN_TEST && nfiles > 1) ? 1 : nfiles;
    for (Int_t i = 0; i < limit; ++i) {
        TChainElement* el = (TChainElement*)tmp.GetListOfFiles()->At(i);
        TString bdt_path = el->GetTitle();

        // Derive xrootd NanoAOD URL from the local deployment path:
        //   local:  .../deploy/MODELTAG/HOST/pnfs/.../nano_N-icenet.root
        //   xrootd: root://HOST//pnfs/.../nano_N.root
        TString grid_part = bdt_path;
        grid_part.Remove(0, strip_prefix.Length());
        Int_t pnfs_pos = grid_part.Index("/pnfs/");
        TString host   = grid_part(0, pnfs_pos);
        TString path   = grid_part(pnfs_pos, grid_part.Length());
        path.ReplaceAll("-icenet.root", ".root");
        TString nano_url = TString::Format("root://%s/%s", host.Data(), path.Data());

        bdt_ch->Add(bdt_path);
        nano_ch->Add(nano_url);
    }
    if (verbose)
        printf("    nano chain: %d  |  bdt chain: %d  %s\n",
               nano_ch->GetListOfFiles()->GetEntries(),
               bdt_ch->GetListOfFiles()->GetEntries(),
               RUN_TEST ? "[TEST: 1 file]" : "");
}

// Draw one 3-column x 4-row mass-distribution canvas for a given model (m_idx)
// and mA group (mA_idx=0 -> SIG[0..3] mA=1.33, mA_idx=1 -> SIG[4..7] mA=0.40).
// Appends the canvas to the already-open PDF; closes the file if is_last=true.
// Each pad shows the leading-muonSV mass (min-chi2 candidate with dlen>0) for
// events with score>=0 (blue) and score>thr (red filled) in the given lxy bin.
// Requires grid proxy for xrootd NanoAOD access.
void make_mass_page(Int_t m_idx, Int_t sig_base, const char* group_label,
                    Float_t mass_max, const TString& pdf,
                    Bool_t is_last, Double_t thr)
{
    Int_t   mass_bins   = 60;

    printf("\n  Mass page: %s | %s\n", MODELS[m_idx].name, group_label);

    // 4 ctau cols x 3 lxy rows, landscape.
    // 2000x1500 px: each pad ~500x500 px (screen-native size; fits in PDF viewer at 1:1).
    TCanvas* cm = new TCanvas(Form("c_mass_%d_%d", m_idx, sig_base),
                               Form("Mass | %s | %s", MODELS[m_idx].name, group_label),
                               2000, 1500);
    cm->Divide(4, 3, 0.002, 0.002);

    // Outer loop: col = ctau (4 values).  Inner loop: row = lxy category (3 bins).
    // Pad numbering for Divide(4,3): pad = row*4 + col + 1.
    for (Int_t col = 0; col < 4; ++col) {          // ctau index -> canvas column
        Int_t s = sig_base + col;
        if (!run_sig(s)) {                          // skip non-selected signal points
            printf("    [col %d] %s (not run -- blank column)\n", col+1, SIG[s].label);
            continue;
        }
        printf("    [col %d] %s\n", col+1, SIG[s].label);

        TChain* nano_ch = new TChain("Events");
        TChain* bdt_ch  = new TChain("Events");
        fill_nano_bdt_chains(MODELS[m_idx].modeltag, SIG[s].dirname,
                             nano_ch, bdt_ch, /*verbose=*/false);

        if (nano_ch->GetListOfFiles()->GetEntries() == 0) {
            printf("      WARNING: no files, skipping\n");
            delete nano_ch; delete bdt_ch; continue;
        }

        nano_ch->AddFriend(bdt_ch);

        ROOT::RDataFrame df_base(*nano_ch);
        ROOT::RDF::RNode df = make_rnode(df_base);

        auto lead_idx_fn = [](const ROOT::RVec<float>& chi2,
                               const ROOT::RVec<float>& dlen) -> int {
            int best = -1; float best_c = 1e9f;
            for (int i = 0; i < (int)chi2.size(); ++i)
                if (dlen[i] > 0.f && chi2[i] < best_c) { best_c = chi2[i]; best = i; }
            return best;
        };
        auto pick_fn = [](int idx, const ROOT::RVec<float>& v) -> float {
            return idx >= 0 ? v[idx] : -1.f;
        };

        auto df2 = df
            .Define("lead_idx",  lead_idx_fn, {"muonSV_chi2", "muonSV_dlen"})
            .Define("lead_dxy",  pick_fn,     {"lead_idx", "muonSV_dxy"})
            .Define("lead_mass", pick_fn,     {"lead_idx", "muonSV_mass"});

        for (Int_t row = 0; row < 3; ++row) {       // lxy category -> canvas row
            Int_t pad = row * 4 + col + 1;
            printf("      [row %d] lxy in [%.0f,%.0f] cm  -> pad %d\n",
                   row+1, (double)LXY_LO[row], (double)LXY_HI[row], pad);

            TString f_lxy = TString::Format(
                "lead_dxy >= %.1ff && lead_dxy < %.1ff",
                LXY_LO[row], LXY_HI[row]);
            TString f_sel = TString::Format(
                "%s >= 0.f && lead_dxy >= %.1ff && lead_dxy < %.1ff",
                BRANCH, LXY_LO[row], LXY_HI[row]);
            TString f_bdt = TString::Format(
                "%s > %.8ff && lead_dxy >= %.1ff && lead_dxy < %.1ff",
                BRANCH, (float)thr, LXY_LO[row], LXY_HI[row]);

            TString hn_s = Form("hm_s_%d_%d_%d_%d", m_idx, sig_base, col, row);
            TString hn_b = Form("hm_b_%d_%d_%d_%d", m_idx, sig_base, col, row);

            printf("      running event loop...\n");
            auto r_sel = df2.Filter(f_sel.Data())
                            .Histo1D({hn_s.Data(), "", mass_bins, 0.f, mass_max}, "lead_mass");
            auto r_bdt = df2.Filter(f_bdt.Data())
                            .Histo1D({hn_b.Data(), "", mass_bins, 0.f, mass_max}, "lead_mass");
            auto c_lxy = df2.Filter(f_lxy.Data()).Count();
            auto c_sel = df2.Filter(f_sel.Data()).Count();
            auto c_bdt = df2.Filter(f_bdt.Data()).Count();

            TH1D* h_sel = (TH1D*)r_sel->Clone(hn_s.Data());
            TH1D* h_bdt = (TH1D*)r_bdt->Clone(hn_b.Data());
            h_sel->SetDirectory(nullptr); h_bdt->SetDirectory(nullptr);

            Long64_t n_lxy = *c_lxy, n_sel = *c_sel, n_bdt = *c_bdt;
            Double_t e_sel = n_lxy > 0 ? (double)n_sel/n_lxy : 0.;
            Double_t e_bdt = n_sel > 0 ? (double)n_bdt/n_sel : 0.;
            Double_t e_tot = n_lxy > 0 ? (double)n_bdt/n_lxy : 0.;
            printf("      N_lxy=%lld  N_sel=%lld  N_bdt=%lld  "
                   "es=%.4f  eb=%.4f  et=%.4f\n",
                   n_lxy, n_sel, n_bdt, e_sel, e_bdt, e_tot);

            cm->cd(pad);
            gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.16);
            gPad->SetTopMargin(0.18);  gPad->SetRightMargin(0.05);

            h_sel->SetLineColor(kBlue+1); h_sel->SetLineWidth(1);
            h_bdt->SetLineColor(kRed+1);  h_bdt->SetLineWidth(1);
            h_bdt->SetFillColorAlpha(kRed-9, 0.5); h_bdt->SetFillStyle(1001);

            h_sel->SetTitle("");
            h_sel->GetXaxis()->SetTitle("muonSV mass (min. #chi^{2}) [GeV]");
            Float_t bin_width = mass_max / mass_bins;
            h_sel->GetYaxis()->SetTitle(Form("Events / %.3g GeV", bin_width));
            h_sel->GetXaxis()->SetTitleFont(43); h_sel->GetXaxis()->SetTitleSize(16);
            h_sel->GetYaxis()->SetTitleFont(43); h_sel->GetYaxis()->SetTitleSize(16);
            h_sel->GetXaxis()->SetLabelFont(43); h_sel->GetXaxis()->SetLabelSize(14);
            h_sel->GetYaxis()->SetLabelFont(43); h_sel->GetYaxis()->SetLabelSize(14);
            h_sel->GetXaxis()->SetTitleOffset(1.1);
            h_sel->GetYaxis()->SetTitleOffset(1.5);
            h_sel->GetXaxis()->SetNdivisions(505);
            gPad->SetLogy();
            h_sel->SetMinimum(0.5);
            h_sel->Draw("hist");
            if (h_bdt->GetMaximum() > 0) h_bdt->Draw("hist same");

            TLatex tl; tl.SetNDC(); tl.SetTextFont(43); tl.SetTextSize(14);
            tl.SetTextAlign(11);
            tl.DrawLatex(0.16, 0.85, SIG[s].label);
            tl.SetTextAlign(31);
            tl.DrawLatex(0.95, 0.85, LXY_LABEL[row]);
            tl.SetTextAlign(11);

            TLegend* lg = new TLegend(0.62, 0.65, 0.87, 0.8);
            lg->SetFillStyle(0); lg->SetBorderSize(0);
            lg->SetTextFont(43); lg->SetTextSize(13);
            lg->AddEntry(h_sel, "presel. (score #geq0)", "l");
            lg->AddEntry(h_bdt,
                Form("BDT cut (thr=%.4f)", thr), "lf");
            lg->AddEntry((TObject*)nullptr, Form("#varepsilon_{sel} = %.4f", e_sel), "");
            lg->AddEntry((TObject*)nullptr, Form("#varepsilon_{BDT} = %.4f", e_bdt), "");
            lg->AddEntry((TObject*)nullptr, Form("#varepsilon_{tot}  = %.4f", e_tot), "");
            lg->Draw("same");
        }

        delete nano_ch; delete bdt_ch;
    }

    cm->Print(is_last ? (pdf + ")").Data() : pdf.Data());
    printf("    page saved (%s)\n", is_last ? "PDF closed" : "page appended");
    delete cm;
}

// --- efficiency-tables (LaTeX) ------------------------------------------------

// Human-readable cut descriptions per model (ROOT TLatex format, for mass-page legends).
static const char* KIN_LABEL[3] = {
    "(any muonSV cand.) Mu10 leg: (p_{T}^{#mu_{1}}>10,|#eta_{1}|<1.2)"
    " OR (p_{T}^{#mu_{2}}>10,|#eta_{2}|<1.2);"
    "  DoubleMu leg: max(p_{T})>4 AND min(p_{T})>3 GeV",
    "(any muonSV cand.) (p_{T}^{#mu_{1}}>10 GeV,|#eta_{1}|<1.2)"
    " OR (p_{T}^{#mu_{2}}>10 GeV,|#eta_{2}|<1.2)",
    "(any muonSV cand.) max(p_{T}^{#mu_{1}},p_{T}^{#mu_{2}})>4 GeV"
    " AND min(p_{T}^{#mu_{1}},p_{T}^{#mu_{2}})>3 GeV"
};
static const char* HLT_LABEL[3] = {
    "HLT_Mu10_Barrel_L1HP11_IP6  OR  HLT_DoubleMu4_3_LowMass",
    "HLT_Mu10_Barrel_L1HP11_IP6",
    "HLT_DoubleMu4_3_LowMass"
};

// Generate one tables PDF + one CSV per model (3 of each) in out_dir.
// Receives the N_tot/N_pre/N_bdt arrays already computed in Phase 2, and internally
// collects N_kin (NanoAOD+BDT) and lxy-split counts for signal, plus QCD counts
// from deployment output (no NanoAOD needed for QCD).
//
// Per-model output:
//   tables_<model>.tex / .pdf -- 3-page LaTeX document compiled with pdflatex
//   tables_<model>.csv        -- all tables in CSV format with quoted section headers
void make_eff_tables_pdf(
    const TString& out_dir,
    Double_t       bdt_thr[N_MODELS],
    Long64_t       N_tot_s[N_SIG][N_MODELS],
    Long64_t       N_pre_s[N_SIG][N_MODELS],   // score >= 0 (kin+trig combined)
    Long64_t       N_bdt_s[N_SIG][N_MODELS])
{
    printf("\n======== Tables: collecting additional data ========\n");

    Long64_t N_kin_s [N_SIG][N_MODELS]    = {};
    Long64_t N_sv_s  [N_SIG][N_MODELS]    = {};  // presel + SV quality (chi2<10, dR<1.2)
    Long64_t N_lxy   [N_SIG][N_MODELS][3] = {};
    Long64_t N_k_lxy [N_SIG][N_MODELS][3] = {};
    Long64_t N_p_lxy [N_SIG][N_MODELS][3] = {};
    Long64_t N_sv_lxy[N_SIG][N_MODELS][3] = {};  // presel + SV quality, per lxy bin
    Long64_t N_b_lxy [N_SIG][N_MODELS][3] = {};
    Long64_t N_tot_q [N_QCD][N_MODELS]    = {};
    Long64_t N_kin_q [N_QCD][N_MODELS]    = {};
    Long64_t N_pre_q [N_QCD][N_MODELS]    = {};
    Long64_t N_sv_q  [N_QCD][N_MODELS]    = {};
    Long64_t N_bdt_q [N_QCD][N_MODELS]    = {};

    // --- significance study (presel+SV+BDT) at 4 FPR working points ----------
    // The existing tables keep using bdt_thr[m] (FPR=1e-4); these extra thresholds
    // are used ONLY for the Asimov-significance table.
    const Double_t FPR_TGT[4] = { 1e-1, 1e-2, 1e-3, 1e-4 };
    TH1D*    sig_sv_cumul[N_SIG][N_MODELS] = {};  // reverse-cumulative BDT score, SV-selected signal
    TH1D*    h_bkg_sv    [N_MODELS]        = {};  // xs/N_tot-weighted SV-selected QCD score
    Double_t sig_thr     [N_MODELS][4]     = {};  // BDT threshold at each FPR target, per model

    auto f_kin_or = [](const ROOT::RVec<float>& m1pt, const ROOT::RVec<float>& m1eta,
                        const ROOT::RVec<float>& m2pt, const ROOT::RVec<float>& m2eta) -> bool {
        bool mu10 = Any((m1pt>10.f && abs(m1eta)<1.2f) || (m2pt>10.f && abs(m2eta)<1.2f));
        auto mx   = Map(m1pt, m2pt, [](float a, float b){ return a>b?a:b; });
        auto mn   = Map(m1pt, m2pt, [](float a, float b){ return a<b?a:b; });
        return mu10 || Any(mx > 4.f && mn > 3.f);
    };
    auto f_kin_mu10 = [](const ROOT::RVec<float>& m1pt, const ROOT::RVec<float>& m1eta,
                          const ROOT::RVec<float>& m2pt, const ROOT::RVec<float>& m2eta) -> bool {
        return Any((m1pt>10.f && abs(m1eta)<1.2f) || (m2pt>10.f && abs(m2eta)<1.2f));
    };
    auto f_kin_dmu = [](const ROOT::RVec<float>& m1pt, const ROOT::RVec<float>& /*m1eta*/,
                         const ROOT::RVec<float>& m2pt, const ROOT::RVec<float>& /*m2eta*/) -> bool {
        auto mx = Map(m1pt, m2pt, [](float a, float b){ return a>b?a:b; });
        auto mn = Map(m1pt, m2pt, [](float a, float b){ return a<b?a:b; });
        return Any(mx > 4.f && mn > 3.f);
    };
    // SV quality: at least one muonSV with chi2 < 10 and dR(mu1,mu2) < 1.2
    auto f_sv = [](const ROOT::RVec<float>& chi2,
                   const ROOT::RVec<float>& m1eta, const ROOT::RVec<float>& m2eta,
                   const ROOT::RVec<float>& m1phi, const ROOT::RVec<float>& m2phi) -> bool {
        for (int i = 0; i < (int)chi2.size(); ++i) {
            if (chi2[i] >= 10.f) continue;
            float deta = m1eta[i] - m2eta[i];
            float dphi = m1phi[i] - m2phi[i];
            while (dphi >  3.14159265f) dphi -= 6.28318530f;
            while (dphi < -3.14159265f) dphi += 6.28318530f;
            if (deta*deta + dphi*dphi < 1.44f) return true;  // dR^2 < 1.2^2
        }
        return false;
    };
    auto lead_fn = [](const ROOT::RVec<float>& chi2, const ROOT::RVec<float>& dlen) -> int {
        int best=-1; float bc=1e9f;
        for (int i=0;i<(int)chi2.size();++i) if(dlen[i]>0.f&&chi2[i]<bc){bc=chi2[i];best=i;}
        return best;
    };
    auto pick_fn = [](int idx, const ROOT::RVec<float>& v) -> float {
        return idx>=0 ? v[idx] : -1.f;
    };

    // collect signal kin + lxy counts via NanoAOD+BDT AddFriend
    for (Int_t m = 0; m < N_MODELS; ++m) {
        if (!run_model(m)) continue;
        printf("\n  -- Model %d/%d: %s (kin+lxy data via NanoAOD) --\n",
               m+1, N_MODELS, MODELS[m].name);
        for (Int_t s = 0; s < N_SIG; ++s) {
            if (!run_sig(s)) continue;
            printf("  [sig %d/%d] %s\n", s+1, N_SIG, SIG[s].label);
            TChain* nano_ch = new TChain("Events");
            TChain* bdt_ch  = new TChain("Events");
            fill_nano_bdt_chains(MODELS[m].modeltag, SIG[s].dirname,
                                  nano_ch, bdt_ch, /*verbose=*/false);
            if (nano_ch->GetListOfFiles()->GetEntries() == 0) {
                printf("    WARNING: no files, skipping\n");
                delete nano_ch; delete bdt_ch; continue;
            }
            nano_ch->AddFriend(bdt_ch);
            ROOT::RDataFrame df_base(*nano_ch);
            ROOT::RDF::RNode df_lxy = make_rnode(df_base)
                .Define("lead_idx", lead_fn, {"muonSV_chi2","muonSV_dlen"})
                .Define("lead_dxy", pick_fn, {"lead_idx","muonSV_dxy"});
            ROOT::RDF::RNode df2 = [&]() -> ROOT::RDF::RNode {
                if (m==0) return df_lxy.Define("kin_pass", f_kin_or,
                    {"muonSV_mu1pt","muonSV_mu1eta","muonSV_mu2pt","muonSV_mu2eta"});
                if (m==1) return df_lxy.Define("kin_pass", f_kin_mu10,
                    {"muonSV_mu1pt","muonSV_mu1eta","muonSV_mu2pt","muonSV_mu2eta"});
                return df_lxy.Define("kin_pass", f_kin_dmu,
                    {"muonSV_mu1pt","muonSV_mu1eta","muonSV_mu2pt","muonSV_mu2eta"});
            }();
            ROOT::RDF::RNode df3 = df2.Define("sv_pass", f_sv,
                {"muonSV_chi2","muonSV_mu1eta","muonSV_mu2eta",
                 "muonSV_mu1phi","muonSV_mu2phi"});
            TString fpre = TString::Format("%s >= 0.", BRANCH);
            TString fsv  = TString::Format("sv_pass && %s >= 0.", BRANCH);
            TString fbdt = TString::Format("%s > %.8f", BRANCH, bdt_thr[m]);
            auto c_kin = df3.Filter("kin_pass").Count();
            auto c_sv  = df3.Filter(fsv.Data()).Count();
            ROOT::RDF::RResultPtr<ULong64_t> c_lxy[3], c_klxy[3], c_plxy[3], c_svlxy[3], c_blxy[3];
            for (Int_t k = 0; k < 3; ++k) {
                TString flxy = TString::Format(
                    "lead_dxy >= %.1ff && lead_dxy < %.1ff", LXY_LO[k], LXY_HI[k]);
                c_lxy  [k] = df3.Filter(flxy.Data()).Count();
                c_klxy [k] = df3.Filter(("kin_pass && " + flxy).Data()).Count();
                c_plxy [k] = df3.Filter((fpre + " && " + flxy).Data()).Count();
                c_svlxy[k] = df3.Filter((fsv  + " && " + flxy).Data()).Count();
                c_blxy [k] = df3.Filter((fbdt + " && " + flxy).Data()).Count();
            }
            // SV-selected BDT score shape (for the significance study)
            auto h_ssv = df3.Filter(fsv.Data())
                            .Histo1D({Form("h_ssv_%d_%d",s,m),"",NBINS,0.,1.}, BRANCH);
            printf("    running event loop (kin + sv + lxy counts)...\n");
            N_kin_s[s][m] = (Long64_t)*c_kin;
            N_sv_s [s][m] = (Long64_t)*c_sv;
            for (Int_t k = 0; k < 3; ++k) {
                N_lxy  [s][m][k]  = (Long64_t)*c_lxy  [k];
                N_k_lxy[s][m][k]  = (Long64_t)*c_klxy [k];
                N_p_lxy[s][m][k]  = (Long64_t)*c_plxy [k];
                N_sv_lxy[s][m][k] = (Long64_t)*c_svlxy[k];
                N_b_lxy[s][m][k]  = (Long64_t)*c_blxy [k];
            }
            printf("    N_kin=%lld  N_sv=%lld  N_lxy=[%lld,%lld,%lld]\n",
                   N_kin_s[s][m], N_sv_s[s][m],
                   N_lxy[s][m][0], N_lxy[s][m][1], N_lxy[s][m][2]);
            // store reverse-cumulative of the SV-selected signal score
            TH1D* hssv = (TH1D*)h_ssv->Clone(Form("hssv_%d_%d",s,m));
            hssv->SetDirectory(nullptr);
            sig_sv_cumul[s][m] = (TH1D*)hssv->GetCumulative(kFALSE);
            sig_sv_cumul[s][m]->SetDirectory(nullptr);
            delete hssv;
            delete nano_ch; delete bdt_ch;
        }
    }

    // collect QCD counts via NanoAOD+BDT (same mechanism as signal)
    printf("\n  Collecting QCD counts via NanoAOD+BDT...\n");
    for (Int_t m = 0; m < N_MODELS; ++m) {
        if (!run_model(m)) continue;
        printf("  -- Model %d/%d: %s --\n", m+1, N_MODELS, MODELS[m].name);
        h_bkg_sv[m] = new TH1D(Form("h_bkg_sv_%d",m),"",NBINS,0.,1.);
        h_bkg_sv[m]->SetDirectory(nullptr);
        for (Int_t q = 0; q < N_QCD; ++q) {
            TChain* nano_ch = new TChain("Events");
            TChain* bdt_ch  = new TChain("Events");
            fill_nano_bdt_chains(MODELS[m].modeltag, QCD[q].dirname,
                                  nano_ch, bdt_ch, /*verbose=*/false);
            if (nano_ch->GetListOfFiles()->GetEntries() == 0) {
                printf("    [QCD %d] WARNING: no files\n", q);
                delete nano_ch; delete bdt_ch; continue;
            }
            nano_ch->AddFriend(bdt_ch);
            ROOT::RDataFrame df_base(*nano_ch);
            ROOT::RDF::RNode df_kin = [&]() -> ROOT::RDF::RNode {
                ROOT::RDF::RNode base = make_rnode(df_base);
                if (m==0) return base.Define("kin_pass", f_kin_or,
                    {"muonSV_mu1pt","muonSV_mu1eta","muonSV_mu2pt","muonSV_mu2eta"});
                if (m==1) return base.Define("kin_pass", f_kin_mu10,
                    {"muonSV_mu1pt","muonSV_mu1eta","muonSV_mu2pt","muonSV_mu2eta"});
                return base.Define("kin_pass", f_kin_dmu,
                    {"muonSV_mu1pt","muonSV_mu1eta","muonSV_mu2pt","muonSV_mu2eta"});
            }();
            ROOT::RDF::RNode df = df_kin.Define("sv_pass", f_sv,
                {"muonSV_chi2","muonSV_mu1eta","muonSV_mu2eta",
                 "muonSV_mu1phi","muonSV_mu2phi"});
            TString fpre = TString::Format("%s >= 0.", BRANCH);
            TString fsv  = TString::Format("sv_pass && %s >= 0.", BRANCH);
            TString fbdt = TString::Format("%s > %.8f", BRANCH, bdt_thr[m]);
            auto c_tot = df.Count();
            auto c_kin = df.Filter("kin_pass").Count();
            auto c_pre = df.Filter(fpre.Data()).Count();
            auto c_sv  = df.Filter(fsv.Data()).Count();
            auto c_bdt = df.Filter(fbdt.Data()).Count();
            auto h_qsv = df.Filter(fsv.Data())
                           .Histo1D({Form("h_qsv_%d_%d",m,q),"",NBINS,0.,1.}, BRANCH);
            N_tot_q[q][m] = (Long64_t)*c_tot;
            N_kin_q[q][m] = (Long64_t)*c_kin;
            N_pre_q[q][m] = (Long64_t)*c_pre;
            N_sv_q [q][m] = (Long64_t)*c_sv;
            N_bdt_q[q][m] = (Long64_t)*c_bdt;
            // accumulate luminosity-yield-weighted (xs/N_tot) SV-selected score
            if (N_tot_q[q][m] > 0) {
                TH1D* hq = (TH1D*)h_qsv->Clone(Form("hq_%d_%d",m,q));
                hq->SetDirectory(nullptr);
                hq->Scale(QCD[q].xs / (Double_t)N_tot_q[q][m]);
                h_bkg_sv[m]->Add(hq);
                delete hq;
            }
            delete nano_ch; delete bdt_ch;
        }
        // background score shape -> reverse cumulative -> FPR thresholds
        if (h_bkg_sv[m]->Integral() > 0.) {
            normalise(h_bkg_sv[m]);
            TH1D* bc = (TH1D*)h_bkg_sv[m]->GetCumulative(kFALSE);
            bc->SetDirectory(nullptr);
            for (Int_t k=0;k<4;++k) { Double_t dum; sig_thr[m][k]=find_thr(bc,FPR_TGT[k],dum); }
            printf("  [signif] thresholds (FPR 1e-1..1e-4): %.4f %.4f %.4f %.4f\n",
                   sig_thr[m][0],sig_thr[m][1],sig_thr[m][2],sig_thr[m][3]);
            delete bc;
        } else {
            for (Int_t k=0;k<4;++k) sig_thr[m][k]=1.;
        }
    }

    // --- output helpers -------------------------------------------------------
    auto eff = [](Long64_t num, Long64_t den) -> TString {
        return den>0 ? TString::Format("%.4f",(Double_t)num/den) : TString("--");
    };
    auto cnt = [](Long64_t n) -> TString { return TString::Format("%lld", n); };

    // Asimov discovery significance Z = sqrt(2[(S+B)ln(1+S/B)-S]); -1 if undefined.
    auto asimov_Z = [](Double_t S, Double_t B) -> Double_t {
        if (B <= 0.) return -1.;
        if (S <= 0.) return 0.;
        return sqrt(2.*((S+B)*log(1.+S/B) - S));
    };
    auto ztex = [](Double_t Z) -> TString {
        return Z < 0. ? TString("--") : TString::Format("%.3g", Z);
    };

    // A ratio in LaTeX scientific notation, 2 decimals (for tiny QCD eps_BDT):
    //   1.23e-4 -> "$1.23\\times10^{-4}$"; 0 -> "$0$"; den<=0 -> "--".
    auto sci_tex = [](Double_t r) -> TString {
        if (r < 0.) return TString("--");
        if (r == 0.) return TString("$0$");
        Int_t ex = (Int_t)floor(log10(r));
        Double_t mn = r / pow(10., ex);
        return TString::Format("$%.2f\\times10^{%d}$", mn, ex);
    };
    auto eff_sci = [&](Long64_t num, Long64_t den) -> TString {
        return den>0 ? sci_tex((Double_t)num/den) : TString("--");
    };
    auto eff_sci_csv = [](Long64_t num, Long64_t den) -> TString {
        return den>0 ? TString::Format("%.2e",(Double_t)num/den) : TString("--");
    };

    Double_t xs_total = 0.;
    for (Int_t q=0; q<N_QCD; ++q) xs_total += QCD[q].xs;

    // 3-line LaTeX column header for signal: $m_\pi=X$ GeV / $m_A=Y$ GeV / $c\tau=Z$ mm
    auto tex_sig_hdr = [](Int_t s) -> TString {
        TString lbl(SIG[s].label);
        Int_t c1=lbl.Index(", "), c2=c1>=0?lbl.Index(", ",c1+1):-1;
        if (c1<0||c2<0) return TString(SIG[s].label);
        // extract numeric values
        TString v1=lbl(lbl.Index("=")+1, c1-lbl.Index("=")-1);
        TString p2=lbl(c1+2,c2-c1-2);
        TString v2=p2(p2.Index("=")+1, p2.Length()-p2.Index("=")-1);
        TString p3=lbl(c2+2,lbl.Length()-c2-2);
        TString v3=p3(p3.Index("=")+1, p3.Length()-p3.Index("=")-1);
        v3.ReplaceAll("mm","");
        return TString::Format(
            "\\shortstack{$m_{\\pi}=%s\\,\\mathrm{GeV}$\\\\$m_{A}=%s\\,\\mathrm{GeV}$\\\\$c\\tau=%s\\,\\mathrm{mm}$}",
            v1.Data(), v2.Data(), v3.Data());
    };

    // single-line LaTeX label for a signal mass set (page titles): m_pi, m_A only.
    // Takes any signal index in the set; ctau is dropped.
    auto tex_group_title = [](Int_t s) -> TString {
        TString lbl(SIG[s].label);
        Int_t c1=lbl.Index(", "), c2=c1>=0?lbl.Index(", ",c1+1):-1;
        if (c1<0||c2<0) return TString(SIG[s].label);
        TString v1=lbl(lbl.Index("=")+1, c1-lbl.Index("=")-1);
        TString p2=lbl(c1+2,c2-c1-2);
        TString v2=p2(p2.Index("=")+1, p2.Length()-p2.Index("=")-1);
        return TString::Format(
            "$m_{\\pi}=%s\\,\\mathrm{GeV}$, $m_{A}=%s\\,\\mathrm{GeV}$",
            v1.Data(), v2.Data());
    };

    // plain-text set label "mpi=X, mA=Y" (no ctau) for CSV section headers
    auto set_label = [](Int_t s) -> TString {
        TString lbl(SIG[s].label);
        Int_t c1=lbl.Index(", ");
        Int_t c2=c1>=0?lbl.Index(", ",c1+1):-1;
        return c2>0 ? lbl(0,c2) : lbl;
    };

    // xs in LaTeX scientific notation: 2799000 -> "2.80\times10^{6}"
    auto fmt_xs = [](Double_t xs) -> TString {
        if (xs <= 0.) return TString("0");
        int ex = (int)floor(log10(xs));
        double mn = xs / pow(10., ex);
        if (ex == 0) return TString::Format("%.3g", xs);
        return TString::Format("%.2f\\times10^{%d}", mn, ex);
    };

    // 2-line LaTeX column header for QCD: p_T range / sigma=X pb
    auto tex_qcd_hdr = [&](Int_t q) -> TString {
        TString s(QCD[q].dirname); Int_t p1=s.Index("PT-"), p2=s.Index("_Fil");
        TString range_tex;
        if (p1>=0&&p2>p1) {
            TString r = s(p1+3, p2-p1-3);
            Int_t pos = r.Index("to");
            if (pos >= 0) {
                TString lo = r(0, pos), hi = r(pos+2, r.Length()-pos-2);
                range_tex = TString::Format(
                    "$p_T\\in[%s\\text{--}%s]\\,\\mathrm{GeV}$", lo.Data(), hi.Data());
            } else {
                range_tex = TString::Format("$p_T>%s\\,\\mathrm{GeV}$", r.Data());
            }
        } else { range_tex = TString(QCD[q].dirname); }
        return TString::Format("\\shortstack{%s\\\\$\\sigma=%s\\,\\mathrm{pb}$}",
                               range_tex.Data(), fmt_xs(QCD[q].xs).Data());
    };

    // Short CSV labels (no LaTeX markup)
    auto csv_qhdr = [](Int_t q) -> TString {
        TString s(QCD[q].dirname); Int_t p1=s.Index("PT-"), p2=s.Index("_Fil");
        if (p1<0||p2<=p1) return TString(QCD[q].dirname);
        TString t=s(p1+3,p2-p1-3); t.ReplaceAll("to","-"); return t;
    };

    // single merged column at end of QCD tables: counts table shows raw total, eff table xs-weighted
    TString total_cnt_hdr =
        "\\shortstack{$\\mathrm{Total}$\\\\$\\mathrm{(weighted)}$}";
    TString total_eff_hdr = TString::Format(
        "\\shortstack{$\\mathrm{Weighted}\\ \\varepsilon$\\\\$\\sigma_{\\mathrm{tot}}=%s\\,\\mathrm{pb}$}",
        fmt_xs(xs_total).Data());

    // Kinematic and HLT cut descriptions in LaTeX (mirrors KIN_LABEL / HLT_LABEL)
    const char* KIN_TEX[3] = {
        "any muonSV cand.: Mu10 leg "
        "$(p_T^{\\mu_1}{>}10,\\,|\\eta_1|{<}1.2)$ or "
        "$(p_T^{\\mu_2}{>}10,\\,|\\eta_2|{<}1.2)$; "
        "DoubleMu leg: $\\max(p_T^{\\mu}){>}4\\,\\mathrm{GeV}$ and "
        "$\\min(p_T^{\\mu}){>}3\\,\\mathrm{GeV}$",
        "any muonSV cand.: "
        "$(p_T^{\\mu_1}{>}10\\,\\mathrm{GeV},\\,|\\eta_1|{<}1.2)$ or "
        "$(p_T^{\\mu_2}{>}10\\,\\mathrm{GeV},\\,|\\eta_2|{<}1.2)$",
        "any muonSV cand.: "
        "$\\max(p_T^{\\mu_1},p_T^{\\mu_2}){>}4\\,\\mathrm{GeV}$ and "
        "$\\min(p_T^{\\mu_1},p_T^{\\mu_2}){>}3\\,\\mathrm{GeV}$"
    };
    const char* HLT_TEX[3] = {
        "\\texttt{HLT\\_Mu10\\_Barrel\\_L1HP11\\_IP6}"
        " OR \\texttt{HLT\\_DoubleMu4\\_3\\_LowMass}",
        "\\texttt{HLT\\_Mu10\\_Barrel\\_L1HP11\\_IP6}",
        "\\texttt{HLT\\_DoubleMu4\\_3\\_LowMass}"
    };

    // --- one LaTeX PDF + CSV per model ----------------------------------------
    for (Int_t m = 0; m < N_MODELS; ++m) {
        if (!run_model(m)) continue;
        TString model_safe = TString(MODELS[m].name); model_safe.ReplaceAll(" ","_");
        // In test mode prefix the file stem with "test_" (after the directory path).
        TString stem = TString(RUN_TEST ? "test_" : "") + "tables_" + model_safe;
        TString tex_path = out_dir + "/" + stem + ".tex";
        TString csv_path = out_dir + "/" + stem + ".csv";

        FILE* fcsv = fopen(csv_path.Data(), "w");
        FILE* ftex = fopen(tex_path.Data(), "w");
        if (!ftex) { printf("ERROR: cannot open %s\n", tex_path.Data()); continue; }

        // Escaped BRANCH name for LaTeX (\texttt{xgb01\_NOJETS\_NOMET})
        TString brtex(BRANCH); brtex.ReplaceAll("_","\\_");
        brtex = TString("\\texttt{") + brtex + "}";

        // xs-weighted QCD efficiencies and xs-weighted average counts.
        //   weighted eff  [step] = Sigma_q (xs_q/xs_total) * eps_q
        //   weighted count[step] = Sigma_q (xs_q/xs_total) * N_step_q
        // The weighted count is a cross-section-weighted average over pT bins:
        // if every bin had the same N, the weighted total equals that N.
        Double_t w_kin=0., w_pre=0., w_sv=0., w_bdt=0.;          // eff wrt N_total (eps_q)
        Double_t wp_hlt=0., wp_sv=0., wp_bdt=0.;                  // eff wrt previous step
        Double_t cw_tot=0., cw_kin=0., cw_pre=0., cw_sv=0., cw_bdt=0.; // for counts
        Long64_t tot_tot=0, tot_kin=0, tot_pre=0, tot_sv=0, tot_bdt=0; // raw sums (unused in tables)
        for (Int_t q=0;q<N_QCD;++q) {
            tot_tot+=N_tot_q[q][m]; tot_kin+=N_kin_q[q][m];
            tot_pre+=N_pre_q[q][m]; tot_sv +=N_sv_q [q][m]; tot_bdt+=N_bdt_q[q][m];
            cw_tot+=QCD[q].xs*(Double_t)N_tot_q[q][m];
            cw_kin+=QCD[q].xs*(Double_t)N_kin_q[q][m];
            cw_pre+=QCD[q].xs*(Double_t)N_pre_q[q][m];
            cw_sv +=QCD[q].xs*(Double_t)N_sv_q [q][m];
            cw_bdt+=QCD[q].xs*(Double_t)N_bdt_q[q][m];
            if (N_tot_q[q][m]>0) {
                w_kin +=QCD[q].xs*(Double_t)N_kin_q[q][m]/N_tot_q[q][m];
                w_pre +=QCD[q].xs*(Double_t)N_pre_q[q][m]/N_tot_q[q][m];
                w_sv  +=QCD[q].xs*(Double_t)N_sv_q [q][m]/N_tot_q[q][m];
                w_bdt +=QCD[q].xs*(Double_t)N_bdt_q[q][m]/N_tot_q[q][m];
            }
            // weighted efficiency relative to the previous cut step
            if (N_kin_q[q][m]>0) wp_hlt+=QCD[q].xs*(Double_t)N_pre_q[q][m]/N_kin_q[q][m];
            if (N_pre_q[q][m]>0) wp_sv +=QCD[q].xs*(Double_t)N_sv_q [q][m]/N_pre_q[q][m];
            if (N_sv_q [q][m]>0) wp_bdt+=QCD[q].xs*(Double_t)N_bdt_q[q][m]/N_sv_q [q][m];
        }
        TString wstr_kin=TString::Format("%.4f",w_kin/xs_total);
        TString wstr_pre=TString::Format("%.4f",w_pre/xs_total);
        TString wstr_sv =TString::Format("%.4f",w_sv /xs_total);
        // eps_BDT is tiny -> scientific notation (2 decimals) for LaTeX and CSV
        TString wstr_bdt     = sci_tex(w_bdt/xs_total);
        TString wstr_bdt_csv = TString::Format("%.2e",w_bdt/xs_total);
        // weighted eff wrt previous step (kin's previous step is N_total, so == wstr_kin)
        TString wprev_hlt    =TString::Format("%.4f",wp_hlt/xs_total);
        TString wprev_sv     =TString::Format("%.4f",wp_sv /xs_total);
        TString wprev_bdt    = sci_tex(wp_bdt/xs_total);
        TString wprev_bdt_csv= TString::Format("%.2e",wp_bdt/xs_total);

        // xs-weighted average counts (shown to 1 decimal)
        Double_t y_tot = cw_tot/xs_total, y_kin = cw_kin/xs_total,
                 y_pre = cw_pre/xs_total, y_sv  = cw_sv /xs_total, y_bdt = cw_bdt/xs_total;
        auto ytex = [ ](Double_t y) -> TString { return TString::Format("%.1f", y); };
        auto ycsv = [ ](Double_t y) -> TString { return TString::Format("%.1f", y); };

        // LaTeX preamble
        fprintf(ftex,
            "\\documentclass[10pt,a4paper]{article}\n"
            "\\usepackage[a4paper,top=1.5cm,bottom=1.5cm,left=1.5cm,right=1.5cm]{geometry}\n"
            "\\usepackage{booktabs}\n"
            "\\usepackage{pdflscape}\n"
            "\\usepackage{graphicx}\n"
            "\\usepackage{amsmath,amssymb}\n"
            "\\usepackage[T1]{fontenc}\n"
            "\\usepackage{helvet}\n"
            "\\renewcommand{\\familydefault}{\\sfdefault}\n"
            "\\setlength{\\parindent}{0pt}\n"
            "\\setlength{\\tabcolsep}{8pt}\n"
            "\\renewcommand{\\arraystretch}{1.4}\n"
            "\\begin{document}\n\n");

        // QCD reference column headers (shown next to the signal columns)
        TString qcd_cnt_hdr =
            "\\shortstack{$\\mathrm{QCD}$\\\\$\\mathrm{Total\\ (weighted)}$}";
        TString qcd_eff_hdr =
            "\\shortstack{$\\mathrm{QCD}$\\\\$\\mathrm{Weighted}\\ \\varepsilon$}";

        // lxy landscape header: 3 lxy group-cells, then 4 ctau headers under each
        auto write_lxy_header = [&](FILE* f, Int_t gb, const char* first_col_label) {
            fprintf(f," & \\multicolumn{4}{c|}"
                "{$l_{xy}\\in[%.0f,\\,%.0f]\\,\\mathrm{cm}$}"
                " & \\multicolumn{4}{c|}"
                "{$l_{xy}\\in[%.0f,\\,%.0f]\\,\\mathrm{cm}$}"
                " & \\multicolumn{4}{c}"
                "{$l_{xy}\\in[%.0f,\\,%.0f]\\,\\mathrm{cm}$}\\\\\n",
                (Double_t)LXY_LO[0],(Double_t)LXY_HI[0],
                (Double_t)LXY_LO[1],(Double_t)LXY_HI[1],
                (Double_t)LXY_LO[2],(Double_t)LXY_HI[2]);
            fprintf(f,"%s", first_col_label);
            for (Int_t k=0;k<3;++k)
                for (Int_t c=0;c<N_CTAU;++c) fprintf(f," & %s",tex_sig_hdr(gb+c).Data());
            fprintf(f," \\\\\n\\midrule\n");
        };

        // ---- One cutflow page + one lxy page per signal mass set ------------
        for (Int_t g = 0; g < N_SETS; ++g) {
            if (!run_set(g)) continue;
            Int_t gb = g * N_CTAU;   // base index of this set's 4 ctau points

            // ===== Cutflow page (4 ctau columns + QCD weighted reference) =====
            fprintf(ftex,
                "{\\large Signal efficiency cutflow"
                " $|$ \\texttt{%s} $|$ %s"
                " $|$ 2024, $13.6\\,\\mathrm{TeV}$}\n\n"
                "\\vspace{0.8em}\n\n"
                "Step 1 -- Kinematic selection: %s\n\n"
                "Step 2 -- Trigger selection: %s\n\n"
                "Step 3 -- SV quality: at least one muonSV with"
                " $\\chi^2<10$ and $\\Delta R(\\mu_1,\\mu_2)<1.2$\n\n"
                "Step 4 -- BDT selection: %s $>%.4f$"
                " (target FPR $=10^{-4}$)\n\n"
                "\\vspace{0.8em}\n\n",
                MODELS[m].name, tex_group_title(gb).Data(), KIN_TEX[m], HLT_TEX[m],
                brtex.Data(), bdt_thr[m]);

            // absolute counts (4 ctau columns + QCD weighted total)
            fprintf(ftex,
                "\\vspace{18pt}\n\n"
                "Absolute event counts\n\n"
                "\\begin{tabular}{cccccc}\n\\toprule\n$\\mathrm{Cut\\ step}$");
            for (Int_t c=0;c<N_CTAU;++c) fprintf(ftex," & %s",tex_sig_hdr(gb+c).Data());
            fprintf(ftex," & %s \\\\\n\\midrule\n$N_{\\mathrm{total}}$", qcd_cnt_hdr.Data());
            for (Int_t c=0;c<N_CTAU;++c) fprintf(ftex," & %lld",N_tot_s[gb+c][m]);
            fprintf(ftex," & %s\\\\\n$\\mathrm{Kin.\\ sel.}$",ytex(y_tot).Data());
            for (Int_t c=0;c<N_CTAU;++c) fprintf(ftex," & %lld",N_kin_s[gb+c][m]);
            fprintf(ftex," & %s\\\\\n$\\mathrm{HLT\\ presel.}$",ytex(y_kin).Data());
            for (Int_t c=0;c<N_CTAU;++c) fprintf(ftex," & %lld",N_pre_s[gb+c][m]);
            fprintf(ftex," & %s\\\\\n$\\mathrm{SV\\ sel.}$ ($\\chi^2{<}10$, $\\Delta R{<}1.2$)",ytex(y_pre).Data());
            for (Int_t c=0;c<N_CTAU;++c) fprintf(ftex," & %lld",N_sv_s[gb+c][m]);
            fprintf(ftex," & %s\\\\\n$\\mathrm{BDT\\ sel.}$",ytex(y_sv).Data());
            for (Int_t c=0;c<N_CTAU;++c) fprintf(ftex," & %lld",N_bdt_s[gb+c][m]);
            fprintf(ftex," & %s\\\\\n\\bottomrule\n\\end{tabular}\n\n",ytex(y_bdt).Data());

            // relative efficiencies (4 ctau columns + QCD weighted eff)
            fprintf(ftex,
                "\\vspace{18pt}\n\n"
                "Relative efficiencies (w.r.t.\\ $N_{\\mathrm{total}}$)\n\n"
                "\\begin{tabular}{cccccc}\n\\toprule\n$\\mathrm{Efficiency}$");
            for (Int_t c=0;c<N_CTAU;++c) fprintf(ftex," & %s",tex_sig_hdr(gb+c).Data());
            fprintf(ftex," & %s \\\\\n\\midrule\n$\\varepsilon_{\\mathrm{kin}}$", qcd_eff_hdr.Data());
            for (Int_t c=0;c<N_CTAU;++c) fprintf(ftex," & %s",eff(N_kin_s[gb+c][m],N_tot_s[gb+c][m]).Data());
            fprintf(ftex," & %s\\\\\n$\\varepsilon_{\\mathrm{HLT}}$", wstr_kin.Data());
            for (Int_t c=0;c<N_CTAU;++c) fprintf(ftex," & %s",eff(N_pre_s[gb+c][m],N_tot_s[gb+c][m]).Data());
            fprintf(ftex," & %s\\\\\n$\\varepsilon_{\\mathrm{SV}}$", wstr_pre.Data());
            for (Int_t c=0;c<N_CTAU;++c) fprintf(ftex," & %s",eff(N_sv_s[gb+c][m],N_tot_s[gb+c][m]).Data());
            fprintf(ftex," & %s\\\\\n$\\varepsilon_{\\mathrm{BDT}}$", wstr_sv.Data());
            for (Int_t c=0;c<N_CTAU;++c) fprintf(ftex," & %s",eff(N_bdt_s[gb+c][m],N_tot_s[gb+c][m]).Data());
            fprintf(ftex," & %s\\\\\n\\bottomrule\n\\end{tabular}\n\n", wstr_bdt.Data());

            // relative efficiencies w.r.t. the previous cut step (4 ctau + QCD weighted)
            fprintf(ftex,
                "\\vspace{18pt}\n\n"
                "Relative efficiencies (w.r.t.\\ previous step)\n\n"
                "\\begin{tabular}{cccccc}\n\\toprule\n$\\mathrm{Efficiency}$");
            for (Int_t c=0;c<N_CTAU;++c) fprintf(ftex," & %s",tex_sig_hdr(gb+c).Data());
            fprintf(ftex," & %s \\\\\n\\midrule\n$\\varepsilon_{\\mathrm{kin}}$", qcd_eff_hdr.Data());
            for (Int_t c=0;c<N_CTAU;++c) fprintf(ftex," & %s",eff(N_kin_s[gb+c][m],N_tot_s[gb+c][m]).Data());
            fprintf(ftex," & %s\\\\\n$\\varepsilon_{\\mathrm{HLT}}$", wstr_kin.Data());
            for (Int_t c=0;c<N_CTAU;++c) fprintf(ftex," & %s",eff(N_pre_s[gb+c][m],N_kin_s[gb+c][m]).Data());
            fprintf(ftex," & %s\\\\\n$\\varepsilon_{\\mathrm{SV}}$", wprev_hlt.Data());
            for (Int_t c=0;c<N_CTAU;++c) fprintf(ftex," & %s",eff(N_sv_s[gb+c][m],N_pre_s[gb+c][m]).Data());
            fprintf(ftex," & %s\\\\\n$\\varepsilon_{\\mathrm{BDT}}$", wprev_sv.Data());
            for (Int_t c=0;c<N_CTAU;++c) fprintf(ftex," & %s",eff(N_bdt_s[gb+c][m],N_sv_s[gb+c][m]).Data());
            fprintf(ftex," & %s\\\\\n\\bottomrule\n\\end{tabular}\n\n", wprev_bdt.Data());

            // 4th table: Asimov significance Z. Col1 = FPR (scientific), col2 = BDT
            // threshold, then one Z column per ctau point (rows = FPR working points).
            fprintf(ftex,
                "\\vspace{18pt}\n\n"
                "Asymptotic significance $Z$ (Asimov; presel.\\,$+$\\,SV\\,$+$\\,BDT,"
                " $\\sigma_{\\mathrm{ggH}}\\!\\times\\!\\mathrm{BR}=%.3g\\,\\mathrm{pb}$,"
                " $L=%.2f\\,\\mathrm{fb}^{-1}$)\n\n"
                "\\begin{tabular}{cccccc}\n\\toprule\n"
                "$\\mathrm{FPR}$ & $\\mathrm{BDT\\ threshold}$",
                SIG_XS*SIG_BR, LUMI/1e3);
            for (Int_t c=0;c<N_CTAU;++c) fprintf(ftex," & %s",tex_sig_hdr(gb+c).Data());
            fprintf(ftex," \\\\\n\\midrule\n");
            for (Int_t k=0;k<4;++k) {
                Double_t B = FPR_TGT[k]*LUMI*w_sv;
                fprintf(ftex,"$10^{-%d}$ & $%.4f$", k+1, sig_thr[m][k]);
                for (Int_t c=0;c<N_CTAU;++c) {
                    Int_t s = gb+c;
                    Double_t epsS = (N_tot_s[s][m]>0 && sig_sv_cumul[s][m])
                        ? sig_sv_cumul[s][m]->GetBinContent(
                              sig_sv_cumul[s][m]->GetXaxis()->FindBin(sig_thr[m][k]))
                          / (Double_t)N_tot_s[s][m]
                        : 0.;
                    Double_t S = LUMI*SIG_XS*SIG_BR*epsS;
                    fprintf(ftex," & %s", ztex(asimov_Z(S,B)).Data());
                }
                fprintf(ftex,"\\\\\n");
            }
            fprintf(ftex,"\\bottomrule\n\\end{tabular}\n\n");
            fprintf(ftex,
                "{\\small Per signal point, $Z=\\sqrt{2\\left[(S{+}B)\\ln(1{+}S/B)-S\\right]}$"
                " with signal yield $S=L\\,\\sigma_{\\mathrm{ggH}}\\,\\mathrm{BR}\\,\\varepsilon_S$"
                " ($\\varepsilon_S=$ fraction of produced signal passing presel.\\,$+$\\,SV\\,$+$\\,BDT)"
                " and background yield $B=\\mathrm{FPR}\\times L\\sum_q\\sigma_q\\,\\varepsilon_q^{\\mathrm{SV}}$"
                " (SV-selected QCD). For each FPR the BDT threshold is the one at which the"
                " cross-section-weighted SV-selected QCD has that false-positive rate.}\n\n");

            if (fcsv) {
                fprintf(fcsv,"\n\"# TABLE: Signal absolute counts | %s | %s\"\n",MODELS[m].name,set_label(gb).Data());
                fprintf(fcsv,"Cut step");
                for (Int_t c=0;c<N_CTAU;++c) fprintf(fcsv,",\"%s\"",SIG[gb+c].label);
                fprintf(fcsv,",QCD Total (weighted)\nN_total");
                for (Int_t c=0;c<N_CTAU;++c) fprintf(fcsv,",%lld",N_tot_s[gb+c][m]);
                fprintf(fcsv,",%s\nKin sel.",ycsv(y_tot).Data());
                for (Int_t c=0;c<N_CTAU;++c) fprintf(fcsv,",%lld",N_kin_s[gb+c][m]);
                fprintf(fcsv,",%s\nHLT presel.",ycsv(y_kin).Data());
                for (Int_t c=0;c<N_CTAU;++c) fprintf(fcsv,",%lld",N_pre_s[gb+c][m]);
                fprintf(fcsv,",%s\nSV sel. (chi2<10 dR<1.2)",ycsv(y_pre).Data());
                for (Int_t c=0;c<N_CTAU;++c) fprintf(fcsv,",%lld",N_sv_s[gb+c][m]);
                fprintf(fcsv,",%s\nBDT sel.",ycsv(y_sv).Data());
                for (Int_t c=0;c<N_CTAU;++c) fprintf(fcsv,",%lld",N_bdt_s[gb+c][m]);
                fprintf(fcsv,",%s\n",ycsv(y_bdt).Data());

                fprintf(fcsv,"\n\"# TABLE: Signal relative effs wrt N_total | %s | %s\"\n",MODELS[m].name,set_label(gb).Data());
                fprintf(fcsv,"Efficiency");
                for (Int_t c=0;c<N_CTAU;++c) fprintf(fcsv,",\"%s\"",SIG[gb+c].label);
                fprintf(fcsv,",QCD Weighted.eff\neps_kin");
                for (Int_t c=0;c<N_CTAU;++c) fprintf(fcsv,",%s",eff(N_kin_s[gb+c][m],N_tot_s[gb+c][m]).Data());
                fprintf(fcsv,",%s\neps_HLT",wstr_kin.Data());
                for (Int_t c=0;c<N_CTAU;++c) fprintf(fcsv,",%s",eff(N_pre_s[gb+c][m],N_tot_s[gb+c][m]).Data());
                fprintf(fcsv,",%s\neps_SV",wstr_pre.Data());
                for (Int_t c=0;c<N_CTAU;++c) fprintf(fcsv,",%s",eff(N_sv_s[gb+c][m],N_tot_s[gb+c][m]).Data());
                fprintf(fcsv,",%s\neps_BDT",wstr_sv.Data());
                for (Int_t c=0;c<N_CTAU;++c) fprintf(fcsv,",%s",eff(N_bdt_s[gb+c][m],N_tot_s[gb+c][m]).Data());
                fprintf(fcsv,",%s\n",wstr_bdt_csv.Data());

                fprintf(fcsv,"\n\"# TABLE: Signal relative effs wrt previous step | %s | %s\"\n",MODELS[m].name,set_label(gb).Data());
                fprintf(fcsv,"Efficiency");
                for (Int_t c=0;c<N_CTAU;++c) fprintf(fcsv,",\"%s\"",SIG[gb+c].label);
                fprintf(fcsv,",QCD Weighted.eff\neps_kin");
                for (Int_t c=0;c<N_CTAU;++c) fprintf(fcsv,",%s",eff(N_kin_s[gb+c][m],N_tot_s[gb+c][m]).Data());
                fprintf(fcsv,",%s\neps_HLT",wstr_kin.Data());
                for (Int_t c=0;c<N_CTAU;++c) fprintf(fcsv,",%s",eff(N_pre_s[gb+c][m],N_kin_s[gb+c][m]).Data());
                fprintf(fcsv,",%s\neps_SV",wprev_hlt.Data());
                for (Int_t c=0;c<N_CTAU;++c) fprintf(fcsv,",%s",eff(N_sv_s[gb+c][m],N_pre_s[gb+c][m]).Data());
                fprintf(fcsv,",%s\neps_BDT",wprev_sv.Data());
                for (Int_t c=0;c<N_CTAU;++c) fprintf(fcsv,",%s",eff(N_bdt_s[gb+c][m],N_sv_s[gb+c][m]).Data());
                fprintf(fcsv,",%s\n",wprev_bdt_csv.Data());

                fprintf(fcsv,"\n\"# TABLE: Asymptotic significance (Asimov) | %s | %s | sigmaxBR=%.4g pb, L=%.2f/fb\"\n",
                    MODELS[m].name,set_label(gb).Data(),SIG_XS*SIG_BR,LUMI/1e3);
                fprintf(fcsv,"FPR,BDT threshold");
                for (Int_t c=0;c<N_CTAU;++c) fprintf(fcsv,",\"%s\"",SIG[gb+c].label);
                fprintf(fcsv,"\n");
                for (Int_t k=0;k<4;++k) {
                    Double_t B = FPR_TGT[k]*LUMI*w_sv;
                    fprintf(fcsv,"%.0e,%.4f",FPR_TGT[k],sig_thr[m][k]);
                    for (Int_t c=0;c<N_CTAU;++c) {
                        Int_t s = gb+c;
                        Double_t epsS = (N_tot_s[s][m]>0 && sig_sv_cumul[s][m])
                            ? sig_sv_cumul[s][m]->GetBinContent(
                                  sig_sv_cumul[s][m]->GetXaxis()->FindBin(sig_thr[m][k]))
                              / (Double_t)N_tot_s[s][m]
                            : 0.;
                        Double_t Z = asimov_Z(LUMI*SIG_XS*SIG_BR*epsS, B);
                        fprintf(fcsv,",%s", Z<0.?TString("--").Data():TString::Format("%.4g",Z).Data());
                    }
                    fprintf(fcsv,"\n");
                }
            }

            fprintf(ftex,"\\clearpage\n\n");

            // ===== lxy page (landscape, 4 ctau x 3 lxy bins) ==================
            fprintf(ftex,
                "\\begin{landscape}\n\n"
                "{\\large Signal efficiencies by $l_{xy}$ bin $|$"
                " \\texttt{%s} $|$ %s $|$ 2024}\n\n"
                "\\vspace{0.6em}\n\n"
                "Efficiencies relative to $N$ in each $l_{xy}$ bin.\n\n",
                MODELS[m].name, tex_group_title(gb).Data());

            fprintf(ftex,"\\vspace{18pt}\n\nAbsolute event counts\n\n"
                "\\resizebox{\\linewidth}{!}{%%\n"
                "\\begin{tabular}{c|cccc|cccc|cccc}\n\\toprule\n");
            write_lxy_header(ftex, gb, "$\\mathrm{Step}$");
            fprintf(ftex,"$N_{l_{xy}}$");
            for (Int_t k=0;k<3;++k) for (Int_t c=0;c<N_CTAU;++c) fprintf(ftex," & %lld",N_lxy[gb+c][m][k]);
            fprintf(ftex,"\\\\\n$\\mathrm{Kin.\\ sel.}$");
            for (Int_t k=0;k<3;++k) for (Int_t c=0;c<N_CTAU;++c) fprintf(ftex," & %lld",N_k_lxy[gb+c][m][k]);
            fprintf(ftex,"\\\\\n$\\mathrm{HLT\\ presel.}$");
            for (Int_t k=0;k<3;++k) for (Int_t c=0;c<N_CTAU;++c) fprintf(ftex," & %lld",N_p_lxy[gb+c][m][k]);
            fprintf(ftex,"\\\\\n$\\mathrm{SV\\ sel.}$ ($\\chi^2{<}10$, $\\Delta R{<}1.2$)");
            for (Int_t k=0;k<3;++k) for (Int_t c=0;c<N_CTAU;++c) fprintf(ftex," & %lld",N_sv_lxy[gb+c][m][k]);
            fprintf(ftex,"\\\\\n$\\mathrm{BDT\\ sel.}$");
            for (Int_t k=0;k<3;++k) for (Int_t c=0;c<N_CTAU;++c) fprintf(ftex," & %lld",N_b_lxy[gb+c][m][k]);
            fprintf(ftex,"\\\\\n\\bottomrule\n\\end{tabular}}\n\n");

            fprintf(ftex,"\\vspace{18pt}\n\nRelative efficiencies (w.r.t.\\ $N$ in $l_{xy}$ bin)\n\n"
                "\\resizebox{\\linewidth}{!}{%%\n"
                "\\begin{tabular}{c|cccc|cccc|cccc}\n\\toprule\n");
            write_lxy_header(ftex, gb, "$\\mathrm{Efficiency}$");
            fprintf(ftex,"$\\varepsilon_{\\mathrm{kin}}$");
            for (Int_t k=0;k<3;++k) for (Int_t c=0;c<N_CTAU;++c) fprintf(ftex," & %s",eff(N_k_lxy[gb+c][m][k],N_lxy[gb+c][m][k]).Data());
            fprintf(ftex,"\\\\\n$\\varepsilon_{\\mathrm{HLT}}$");
            for (Int_t k=0;k<3;++k) for (Int_t c=0;c<N_CTAU;++c) fprintf(ftex," & %s",eff(N_p_lxy[gb+c][m][k],N_lxy[gb+c][m][k]).Data());
            fprintf(ftex,"\\\\\n$\\varepsilon_{\\mathrm{SV}}$");
            for (Int_t k=0;k<3;++k) for (Int_t c=0;c<N_CTAU;++c) fprintf(ftex," & %s",eff(N_sv_lxy[gb+c][m][k],N_lxy[gb+c][m][k]).Data());
            fprintf(ftex,"\\\\\n$\\varepsilon_{\\mathrm{BDT}}$");
            for (Int_t k=0;k<3;++k) for (Int_t c=0;c<N_CTAU;++c) fprintf(ftex," & %s",eff(N_b_lxy[gb+c][m][k],N_lxy[gb+c][m][k]).Data());
            fprintf(ftex,"\\\\\n\\bottomrule\n\\end{tabular}}\n\n");

            fprintf(ftex,"\\vspace{18pt}\n\nRelative efficiencies (w.r.t.\\ previous step)\n\n"
                "\\resizebox{\\linewidth}{!}{%%\n"
                "\\begin{tabular}{c|cccc|cccc|cccc}\n\\toprule\n");
            write_lxy_header(ftex, gb, "$\\mathrm{Efficiency}$");
            fprintf(ftex,"$\\varepsilon_{\\mathrm{kin}}$");
            for (Int_t k=0;k<3;++k) for (Int_t c=0;c<N_CTAU;++c) fprintf(ftex," & %s",eff(N_k_lxy[gb+c][m][k],N_lxy[gb+c][m][k]).Data());
            fprintf(ftex,"\\\\\n$\\varepsilon_{\\mathrm{HLT}}$");
            for (Int_t k=0;k<3;++k) for (Int_t c=0;c<N_CTAU;++c) fprintf(ftex," & %s",eff(N_p_lxy[gb+c][m][k],N_k_lxy[gb+c][m][k]).Data());
            fprintf(ftex,"\\\\\n$\\varepsilon_{\\mathrm{SV}}$");
            for (Int_t k=0;k<3;++k) for (Int_t c=0;c<N_CTAU;++c) fprintf(ftex," & %s",eff(N_sv_lxy[gb+c][m][k],N_p_lxy[gb+c][m][k]).Data());
            fprintf(ftex,"\\\\\n$\\varepsilon_{\\mathrm{BDT}}$");
            for (Int_t k=0;k<3;++k) for (Int_t c=0;c<N_CTAU;++c) fprintf(ftex," & %s",eff(N_b_lxy[gb+c][m][k],N_sv_lxy[gb+c][m][k]).Data());
            fprintf(ftex,"\\\\\n\\bottomrule\n\\end{tabular}}\n\n");

            fprintf(ftex,"\\end{landscape}\n\\clearpage\n\n");

            // CSV: 3 lxy-bin tables for this set (4 ctau columns each)
            if (fcsv) {
                for (Int_t k=0;k<3;++k) {
                    fprintf(fcsv,"\n\"# TABLE: lxy=[%.0f,%.0f]cm efficiencies | %s | %s\"\n",
                        (Double_t)LXY_LO[k],(Double_t)LXY_HI[k],MODELS[m].name,set_label(gb).Data());
                    fprintf(fcsv,"Step");
                    for (Int_t c=0;c<N_CTAU;++c) fprintf(fcsv,",\"%s\"",SIG[gb+c].label);
                    fprintf(fcsv,"\nN_lxy");
                    for (Int_t c=0;c<N_CTAU;++c) fprintf(fcsv,",%lld",N_lxy[gb+c][m][k]);
                    fprintf(fcsv,"\neps_kin");
                    for (Int_t c=0;c<N_CTAU;++c) fprintf(fcsv,",%s",eff(N_k_lxy[gb+c][m][k],N_lxy[gb+c][m][k]).Data());
                    fprintf(fcsv,"\neps_HLT");
                    for (Int_t c=0;c<N_CTAU;++c) fprintf(fcsv,",%s",eff(N_p_lxy[gb+c][m][k],N_lxy[gb+c][m][k]).Data());
                    fprintf(fcsv,"\neps_SV");
                    for (Int_t c=0;c<N_CTAU;++c) fprintf(fcsv,",%s",eff(N_sv_lxy[gb+c][m][k],N_lxy[gb+c][m][k]).Data());
                    fprintf(fcsv,"\neps_BDT");
                    for (Int_t c=0;c<N_CTAU;++c) fprintf(fcsv,",%s",eff(N_b_lxy[gb+c][m][k],N_lxy[gb+c][m][k]).Data());
                    // efficiencies relative to the previous step (same lxy bin)
                    fprintf(fcsv,"\neps_kin_prev");
                    for (Int_t c=0;c<N_CTAU;++c) fprintf(fcsv,",%s",eff(N_k_lxy[gb+c][m][k],N_lxy[gb+c][m][k]).Data());
                    fprintf(fcsv,"\neps_HLT_prev");
                    for (Int_t c=0;c<N_CTAU;++c) fprintf(fcsv,",%s",eff(N_p_lxy[gb+c][m][k],N_k_lxy[gb+c][m][k]).Data());
                    fprintf(fcsv,"\neps_SV_prev");
                    for (Int_t c=0;c<N_CTAU;++c) fprintf(fcsv,",%s",eff(N_sv_lxy[gb+c][m][k],N_p_lxy[gb+c][m][k]).Data());
                    fprintf(fcsv,"\neps_BDT_prev");
                    for (Int_t c=0;c<N_CTAU;++c) fprintf(fcsv,",%s",eff(N_b_lxy[gb+c][m][k],N_sv_lxy[gb+c][m][k]).Data());
                    fprintf(fcsv,"\n");
                }
            }
        }

        // ---- QCD background page (landscape, shared, last page) -------------
        TString qcd_colspec = "c";          // all cells centered
        for (Int_t q=0;q<N_QCD+1;++q) qcd_colspec += "c";

        fprintf(ftex,
            "\\begin{landscape}\n\n"
            "{\\large QCD background cutflow $|$ \\texttt{%s}"
            " $|$ 2024, $13.6\\,\\mathrm{TeV}$}\n\n"
            "\\vspace{0.8em}\n\n"
            "\\textit{Column headers: $p_T$ range [GeV] and $\\sigma$ [pb].}\n\n"
            "Preselection: %s $\\geq 0$"
            "\\quad BDT: %s $>%.4f$ (FPR$=10^{-4}$)\n\n"
            "\\vspace{0.8em}\n\n",
            MODELS[m].name, brtex.Data(), brtex.Data(), bdt_thr[m]);

        // QCD absolute counts
        fprintf(ftex,
            "\\vspace{18pt}\n\n"
            "Absolute event counts;"
            " Total\\ (weighted)$=\\sum_q(\\sigma_q N_q)/\\sum_q\\sigma_q$\n\n"
            "\\resizebox{\\linewidth}{!}{%%\n"
            "\\begin{tabular}{%s}\n\\toprule\nCut step",
            qcd_colspec.Data());
        for (Int_t q=0;q<N_QCD;++q) fprintf(ftex," & %s",tex_qcd_hdr(q).Data());
        fprintf(ftex," & %s\\\\\n\\midrule\n", total_cnt_hdr.Data());
        fprintf(ftex,"$N_{\\mathrm{total}}$");
        for (Int_t q=0;q<N_QCD;++q) fprintf(ftex," & %lld",N_tot_q[q][m]);
        fprintf(ftex," & %s\\\\\n",ytex(y_tot).Data());
        fprintf(ftex,"$\\mathrm{Kin.\\ sel.}$");
        for (Int_t q=0;q<N_QCD;++q) fprintf(ftex," & %lld",N_kin_q[q][m]);
        fprintf(ftex," & %s\\\\\n",ytex(y_kin).Data());
        fprintf(ftex,"$\\mathrm{HLT\\ presel.}$");
        for (Int_t q=0;q<N_QCD;++q) fprintf(ftex," & %lld",N_pre_q[q][m]);
        fprintf(ftex," & %s\\\\\n",ytex(y_pre).Data());
        fprintf(ftex,"$\\mathrm{SV\\ sel.}$ ($\\chi^2{<}10$, $\\Delta R{<}1.2$)");
        for (Int_t q=0;q<N_QCD;++q) fprintf(ftex," & %lld",N_sv_q[q][m]);
        fprintf(ftex," & %s\\\\\n",ytex(y_sv).Data());
        fprintf(ftex,"$\\mathrm{BDT\\ sel.}$");
        for (Int_t q=0;q<N_QCD;++q) fprintf(ftex," & %lld",N_bdt_q[q][m]);
        fprintf(ftex," & %s\\\\\n\\bottomrule\n\\end{tabular}}\n\n",ytex(y_bdt).Data());

        if (fcsv) {
            fprintf(fcsv,"\n\"# TABLE: QCD absolute counts | %s\"\n",MODELS[m].name);
            fprintf(fcsv,"Cut step");
            for (Int_t q=0;q<N_QCD;++q) fprintf(fcsv,",%s",csv_qhdr(q).Data());
            fprintf(fcsv,",Total (weighted)\nN_total");
            for (Int_t q=0;q<N_QCD;++q) fprintf(fcsv,",%lld",N_tot_q[q][m]);
            fprintf(fcsv,",%s\nKin sel.",ycsv(y_tot).Data());
            for (Int_t q=0;q<N_QCD;++q) fprintf(fcsv,",%lld",N_kin_q[q][m]);
            fprintf(fcsv,",%s\nHLT presel.",ycsv(y_kin).Data());
            for (Int_t q=0;q<N_QCD;++q) fprintf(fcsv,",%lld",N_pre_q[q][m]);
            fprintf(fcsv,",%s\nSV sel.",ycsv(y_pre).Data());
            for (Int_t q=0;q<N_QCD;++q) fprintf(fcsv,",%lld",N_sv_q[q][m]);
            fprintf(fcsv,",%s\nBDT sel.",ycsv(y_sv).Data());
            for (Int_t q=0;q<N_QCD;++q) fprintf(fcsv,",%lld",N_bdt_q[q][m]);
            fprintf(fcsv,",%s\n",ycsv(y_bdt).Data());
        }

        // QCD relative efficiencies
        fprintf(ftex,
            "\\vspace{18pt}\n\n"
            "Relative efficiencies (w.r.t.\\ $N_{\\mathrm{total}}$);"
            " Weighted\\ $\\varepsilon=\\sum(\\sigma_q\\varepsilon_q)/\\sum\\sigma_q$\n\n"
            "\\resizebox{\\linewidth}{!}{%%\n"
            "\\begin{tabular}{%s}\n\\toprule\n$\\mathrm{Efficiency}$",
            qcd_colspec.Data());
        for (Int_t q=0;q<N_QCD;++q) fprintf(ftex," & %s",tex_qcd_hdr(q).Data());
        fprintf(ftex," & %s\\\\\n\\midrule\n", total_eff_hdr.Data());
        fprintf(ftex,"$\\varepsilon_{\\mathrm{kin}}$");
        for (Int_t q=0;q<N_QCD;++q) fprintf(ftex," & %s",eff(N_kin_q[q][m],N_tot_q[q][m]).Data());
        fprintf(ftex," & %s\\\\\n",wstr_kin.Data());
        fprintf(ftex,"$\\varepsilon_{\\mathrm{HLT}}$");
        for (Int_t q=0;q<N_QCD;++q) fprintf(ftex," & %s",eff(N_pre_q[q][m],N_tot_q[q][m]).Data());
        fprintf(ftex," & %s\\\\\n",wstr_pre.Data());
        fprintf(ftex,"$\\varepsilon_{\\mathrm{SV}}$");
        for (Int_t q=0;q<N_QCD;++q) fprintf(ftex," & %s",eff(N_sv_q[q][m],N_tot_q[q][m]).Data());
        fprintf(ftex," & %s\\\\\n",wstr_sv.Data());
        fprintf(ftex,"$\\varepsilon_{\\mathrm{BDT}}$");
        for (Int_t q=0;q<N_QCD;++q) fprintf(ftex," & %s",eff_sci(N_bdt_q[q][m],N_tot_q[q][m]).Data());
        fprintf(ftex," & %s\\\\\n\\bottomrule\n\\end{tabular}}\n\n", wstr_bdt.Data());

        // QCD relative efficiencies w.r.t. previous step
        fprintf(ftex,
            "\\vspace{18pt}\n\n"
            "Relative efficiencies (w.r.t.\\ previous step);"
            " Weighted\\ $\\varepsilon=\\sum(\\sigma_q\\varepsilon_q)/\\sum\\sigma_q$\n\n"
            "\\resizebox{\\linewidth}{!}{%%\n"
            "\\begin{tabular}{%s}\n\\toprule\n$\\mathrm{Efficiency}$",
            qcd_colspec.Data());
        for (Int_t q=0;q<N_QCD;++q) fprintf(ftex," & %s",tex_qcd_hdr(q).Data());
        fprintf(ftex," & %s\\\\\n\\midrule\n", total_eff_hdr.Data());
        fprintf(ftex,"$\\varepsilon_{\\mathrm{kin}}$");
        for (Int_t q=0;q<N_QCD;++q) fprintf(ftex," & %s",eff(N_kin_q[q][m],N_tot_q[q][m]).Data());
        fprintf(ftex," & %s\\\\\n",wstr_kin.Data());
        fprintf(ftex,"$\\varepsilon_{\\mathrm{HLT}}$");
        for (Int_t q=0;q<N_QCD;++q) fprintf(ftex," & %s",eff(N_pre_q[q][m],N_kin_q[q][m]).Data());
        fprintf(ftex," & %s\\\\\n",wprev_hlt.Data());
        fprintf(ftex,"$\\varepsilon_{\\mathrm{SV}}$");
        for (Int_t q=0;q<N_QCD;++q) fprintf(ftex," & %s",eff(N_sv_q[q][m],N_pre_q[q][m]).Data());
        fprintf(ftex," & %s\\\\\n",wprev_sv.Data());
        fprintf(ftex,"$\\varepsilon_{\\mathrm{BDT}}$");
        for (Int_t q=0;q<N_QCD;++q) fprintf(ftex," & %s",eff_sci(N_bdt_q[q][m],N_sv_q[q][m]).Data());
        fprintf(ftex," & %s\\\\\n\\bottomrule\n\\end{tabular}}\n\n", wprev_bdt.Data());

        if (fcsv) {
            fprintf(fcsv,"\n\"# TABLE: QCD relative efficiencies wrt N_total | %s\"\n",MODELS[m].name);
            fprintf(fcsv,"Efficiency");
            for (Int_t q=0;q<N_QCD;++q) fprintf(fcsv,",%s",csv_qhdr(q).Data());
            fprintf(fcsv,",Weighted.eff\neps_kin");
            for (Int_t q=0;q<N_QCD;++q) fprintf(fcsv,",%s",eff(N_kin_q[q][m],N_tot_q[q][m]).Data());
            fprintf(fcsv,",%s\neps_HLT",wstr_kin.Data());
            for (Int_t q=0;q<N_QCD;++q) fprintf(fcsv,",%s",eff(N_pre_q[q][m],N_tot_q[q][m]).Data());
            fprintf(fcsv,",%s\neps_SV",wstr_pre.Data());
            for (Int_t q=0;q<N_QCD;++q) fprintf(fcsv,",%s",eff(N_sv_q[q][m],N_tot_q[q][m]).Data());
            fprintf(fcsv,",%s\neps_BDT",wstr_sv.Data());
            for (Int_t q=0;q<N_QCD;++q) fprintf(fcsv,",%s",eff_sci_csv(N_bdt_q[q][m],N_tot_q[q][m]).Data());
            fprintf(fcsv,",%s\n",wstr_bdt_csv.Data());

            fprintf(fcsv,"\n\"# TABLE: QCD relative efficiencies wrt previous step | %s\"\n",MODELS[m].name);
            fprintf(fcsv,"Efficiency");
            for (Int_t q=0;q<N_QCD;++q) fprintf(fcsv,",%s",csv_qhdr(q).Data());
            fprintf(fcsv,",Weighted.eff\neps_kin");
            for (Int_t q=0;q<N_QCD;++q) fprintf(fcsv,",%s",eff(N_kin_q[q][m],N_tot_q[q][m]).Data());
            fprintf(fcsv,",%s\neps_HLT",wstr_kin.Data());
            for (Int_t q=0;q<N_QCD;++q) fprintf(fcsv,",%s",eff(N_pre_q[q][m],N_kin_q[q][m]).Data());
            fprintf(fcsv,",%s\neps_SV",wprev_hlt.Data());
            for (Int_t q=0;q<N_QCD;++q) fprintf(fcsv,",%s",eff(N_sv_q[q][m],N_pre_q[q][m]).Data());
            fprintf(fcsv,",%s\neps_BDT",wprev_sv.Data());
            for (Int_t q=0;q<N_QCD;++q) fprintf(fcsv,",%s",eff_sci_csv(N_bdt_q[q][m],N_sv_q[q][m]).Data());
            fprintf(fcsv,",%s\n",wprev_bdt_csv.Data());
        }

        fprintf(ftex,"\\end{landscape}\n\\clearpage\n\n");

        // (lxy-split tables are emitted per signal point above, before this QCD page)

        fprintf(ftex,"\\end{document}\n");
        fclose(ftex);
        if (fcsv) fclose(fcsv);

        // Compile LaTeX -> PDF
        TString compile_cmd = TString::Format(
            "cd '%s' && pdflatex -interaction=nonstopmode %s.tex"
            " > %s.log 2>&1",
            out_dir.Data(), stem.Data(), stem.Data());
        printf("  Compiling %s\n", tex_path.Data());
        Int_t ret = gSystem->Exec(compile_cmd.Data());
        if (ret != 0)
            printf("  WARNING: pdflatex failed (ret=%d) -- see %s/%s.log\n",
                   ret, out_dir.Data(), stem.Data());
        else
            printf("  -> %s/%s.pdf\n  -> %s\n",
                   out_dir.Data(), stem.Data(), csv_path.Data());
    }

    // clean up significance-study histograms
    for (Int_t m=0;m<N_MODELS;++m) {
        delete h_bkg_sv[m];
        for (Int_t s=0;s<N_SIG;++s) delete sig_sv_cumul[s][m];
    }
}

// --- main ---------------------------------------------------------------------

void run_analysis()
{
    gStyle->SetOptStat(0);
    gStyle->SetTextFont(42);

    // Phase 1: BDT threshold per model
    // By default, hardcoded thresholds (FPR=1e-4 from a prior full QCD study)
    // are used and Phase 1 is skipped.  Pass "compute_thr" to recompute from data.

    Double_t bdt_thr[N_MODELS];
    Double_t bkg_eff[N_MODELS];
    TH1D*    bkg_c  [N_MODELS];   // nullptr when Phase 1 is skipped

    for (Int_t m = 0; m < N_MODELS; ++m) {
        bkg_c[m]  = nullptr;
        bdt_thr[m] = HARDCODED_THR[m];
        bkg_eff[m] = TARGET;       // definition; overwritten if COMPUTE_THR
    }

    if (!COMPUTE_THR) {
        printf("\n======== Phase 1: using hardcoded BDT thresholds (FPR=1e-4) ========\n");
        for (Int_t m = 0; m < N_MODELS; ++m)
            printf("  %-18s  score > %.4f\n", MODELS[m].name, bdt_thr[m]);
        printf("  (pass \"compute_thr\" argument to recompute from QCD data)\n");

    } else {
        printf("\n======== Phase 1: computing BDT thresholds from QCD background ========\n");
        if (RUN_TEST)
            printf("  [TEST MODE] event loop limited to %lld events per sample\n", TEST_EVTS);

        for (Int_t m = 0; m < N_MODELS; ++m) {
            if (!run_model(m)) continue;
            printf("\n  -- Model %d/%d: %s --\n", m+1, N_MODELS, MODELS[m].name);

            TH1D* h_bkg = new TH1D(Form("h_bkg_%d", m), "bkg", NBINS, 0., 1.);
            h_bkg->SetDirectory(nullptr);

            for (Int_t q = 0; q < N_QCD; ++q) {
                printf("\n  [QCD %2d/%d] %s\n", q+1, N_QCD, QCD[q].dirname);
                TChain* ch = make_chain(MODELS[m].modeltag, QCD[q].dirname);
                if (ch->GetListOfFiles()->GetEntries() == 0) {
                    printf("    WARNING: skipping -- no files found\n");
                    delete ch; continue;
                }
                // Each event is weighted by the bin cross section sigma_q only
                // (no /N and no luminosity): the reverse-cumulative fraction then
                // reads Sigma_q sigma_q N_q(>T) / Sigma_q sigma_q N_q(presel), the
                // cross-section-weighted QCD event count -- the same weighting used
                // by the "Total (weighted)" QCD column.  GetEntries() is only needed
                // to detect empty chains (skipped in test mode, where it is slow).
                if (!RUN_TEST) {
                    Double_t n = (Double_t)ch->GetEntries();
                    if (n == 0.) {
                        printf("    WARNING: skipping -- no entries found\n");
                        delete ch; continue;
                    }
                    printf("    entries: %.0f  xs: %.3g pb  weight: sigma_q = %.3g\n",
                           n, QCD[q].xs, QCD[q].xs);
                } else {
                    printf("    xs: %.3g pb  weight: sigma_q = %.3g  (test mode)\n",
                           QCD[q].xs, QCD[q].xs);
                }

                ROOT::RDataFrame df_all(*ch);
                ROOT::RDF::RNode df = make_rnode(df_all);
                TString hn = Form("h_qcd_%d_%d", m, q);
                printf("    running event loop...\n");
                auto h = df.Histo1D({hn.Data(), hn.Data(), NBINS, 0., 1.}, BRANCH);
                TH1D* hq = (TH1D*)h->Clone(hn.Data());
                hq->SetDirectory(nullptr);
                hq->Scale(QCD[q].xs);  // cross-section weight (matches Total column scheme)
                printf("    done. weighted integral: %.3e\n", hq->Integral());
                h_bkg->Add(hq);
                delete hq; delete ch;
            }

            printf("\n  Building background cumulative histogram...\n");
            normalise(h_bkg);
            bkg_c[m] = (TH1D*)h_bkg->GetCumulative(kFALSE);
            bkg_c[m]->SetDirectory(nullptr);
            delete h_bkg;

            bdt_thr[m] = find_thr(bkg_c[m], TARGET, bkg_eff[m]);
            printf("  >>> Threshold at %.0e bkg eff: score > %.4f  (actual bkg eff = %.2e)\n",
                   TARGET, bdt_thr[m], bkg_eff[m]);
        }
    }

    // Phase 2: signal cutflow counts AND efficiency histograms

    printf("\n======== Phase 2: signal cutflow ========\n");
    if (RUN_TEST)
        printf("  [TEST MODE] N_tot is a placeholder (%lld); efficiencies are approximate\n",
               TEST_EVTS);

    // Zero-initialised so models skipped via ONLY_MODEL stay at 0 (not garbage).
    Long64_t N_tot[N_SIG][N_MODELS] = {};
    Long64_t N_sel[N_SIG][N_MODELS] = {};  // score >= 0  (passed trigger+kin selection)
    Long64_t N_bdt[N_SIG][N_MODELS] = {};  // score > threshold
    // Normalised reverse-cumulative histograms for plotting
    TH1D* sig_c[N_SIG][N_MODELS] = {};

    for (Int_t m = 0; m < N_MODELS; ++m) {
        if (!run_model(m)) continue;
        printf("\n  -- Model %d/%d: %s  (thr=%.4f) --\n",
               m+1, N_MODELS, MODELS[m].name, bdt_thr[m]);

        for (Int_t s = 0; s < N_SIG; ++s) {
            if (!run_sig(s)) continue;
            printf("\n  [sig %d/%d] %s [%s]\n", s+1, N_SIG, SIG[s].label, SIG[s].lxy);

            TChain* ch = make_chain(MODELS[m].modeltag, SIG[s].dirname);

            // In test mode skip the slow GetEntries() scan and use a placeholder
            if (RUN_TEST) {
                N_tot[s][m] = TEST_EVTS;
                printf("    [TEST] N_tot placeholder = %lld\n", N_tot[s][m]);
            } else {
                printf("    counting total entries (opens all files)...\n");
                N_tot[s][m] = ch->GetEntries();
                printf("    N_tot = %lld\n", N_tot[s][m]);
            }

            if (ch->GetListOfFiles()->GetEntries() == 0) {
                printf("    WARNING: empty chain -- skipping\n");
                N_sel[s][m] = N_bdt[s][m] = 0;
                sig_c[s][m] = nullptr;
                delete ch; continue;
            }

            ROOT::RDataFrame df_all(*ch);
            ROOT::RDF::RNode df = make_rnode(df_all);

            // Lazy event counts + BDT score histogram (single event-loop pass)
            TString fsel = TString::Format("%s >= 0.", BRANCH);
            TString fbdt = TString::Format("%s > %.8f", BRANCH, bdt_thr[m]);
            auto cnt_sel = df.Filter(fsel.Data()).Count();
            auto cnt_bdt = df.Filter(fbdt.Data()).Count();

            TString hn = Form("h_sig_%d_%d", s, m);
            auto h_raw = df.Histo1D({hn.Data(), hn.Data(), NBINS, 0., 1.}, BRANCH);

            printf("    running event loop (counts + histogram)...\n");
            // Dereferencing h_raw triggers the event loop, computing all three
            // lazy actions (cnt_sel, cnt_bdt, h_raw) in a single pass.
            TH1D* hs = (TH1D*)h_raw->Clone(hn.Data());
            hs->SetDirectory(nullptr);

            N_sel[s][m] = (Long64_t)(*cnt_sel);
            N_bdt[s][m] = (Long64_t)(*cnt_bdt);

            Double_t e_sel = N_tot[s][m] > 0 ? (Double_t)N_sel[s][m]/N_tot[s][m] : 0.;
            Double_t e_bdt = N_sel[s][m] > 0 ? (Double_t)N_bdt[s][m]/N_sel[s][m] : 0.;
            Double_t e_tot = N_tot[s][m] > 0 ? (Double_t)N_bdt[s][m]/N_tot[s][m] : 0.;
            printf("    N_sel=%lld  N_bdt=%lld  "
                   "es=%.4f  eb=%.4f  et=%.4f\n",
                   N_sel[s][m], N_bdt[s][m], e_sel, e_bdt, e_tot);

            normalise(hs);
            sig_c[s][m] = (TH1D*)hs->GetCumulative(kFALSE);
            sig_c[s][m]->SetDirectory(nullptr);
            delete hs; delete ch;
        }
    }

    // Phase 3: print formatted cutflow table

    const char* SEP = "============================================================"
                      "======================================================\n";
    const char* sep = "------------------------------------------------------------"
                      "------------------------------------------------------\n";

    printf("\n%s", SEP);
    printf("  noMET signal efficiency cutflow  |  5 signal mass sets  |  2024\n");
    printf("  L_int = %.2f fb^-1  |  BDT branch: %s  |  target bkg eff = %.0e\n",
           LUMI / 1e3, BRANCH, TARGET);
    printf("%s", SEP);
    printf("  lxy displacement categories (mean lxy ~ ctau x 4-6 for mpi=4 GeV):\n");
    printf("    LOW   ctau=0.1mm  mean lxy ~0.05cm   >99.9%% dimuons below 1 cm\n");
    printf("    MED   ctau=1mm    mean lxy ~0.5cm    transitional, mostly below 1 cm\n");
    printf("    HIGH  ctau=10mm   mean lxy ~5cm      majority above 1 cm\n");
    printf("    XHIGH ctau=100mm  mean lxy ~50cm     extreme displacement\n");
    printf("%s", sep);
    printf("  BDT thresholds:\n");
    for (Int_t m = 0; m < N_MODELS; ++m) {
        if (run_model(m))
            printf("    %-18s  score > %.4f  (bkg eff = %.2e)\n",
                   MODELS[m].name, bdt_thr[m], bkg_eff[m]);
        else
            printf("    %-18s  (model not run)\n", MODELS[m].name);
    }
    printf("%s", SEP);

    // Header
    printf("  %-22s %5s | %-28s| %-28s| %-28s\n",
           "Signal point", "lxy",
           " Mu10 OR DoubleMu", " Mu10", " DoubleMu");
    printf("  %-22s %5s |%8s%8s%8s  |%8s%8s%8s  |%8s%8s%8s\n",
           "", "", "Ntot/k","Nsel/k","Nbdt/k",
                   "Ntot/k","Nsel/k","Nbdt/k",
                   "Ntot/k","Nsel/k","Nbdt/k");
    printf("%s", sep);

    for (Int_t s = 0; s < N_SIG; ++s) {
        if (!run_sig(s)) {
            printf("  %-22s %5s |  (signal point not run)\n", SIG[s].label, SIG[s].lxy);
            continue;
        }
        // Event count row
        printf("  %-22s %5s |", SIG[s].label, SIG[s].lxy);
        for (Int_t m = 0; m < N_MODELS; ++m) {
            if (run_model(m))
                printf(" %7.1f %7.1f %7.1f  ",
                       N_tot[s][m] / 1e3,
                       N_sel[s][m] / 1e3,
                       N_bdt[s][m] / 1e3);
            else
                printf(" %7s %7s %7s  ", "--", "--", "--");
            if (m < N_MODELS - 1) printf("|");
        }
        printf("\n");

        // Efficiency row
        printf("  %-22s %5s |", "(efficiencies)", "");
        for (Int_t m = 0; m < N_MODELS; ++m) {
            if (run_model(m)) {
                Double_t e_sel = N_tot[s][m] > 0 ? (Double_t)N_sel[s][m]/N_tot[s][m] : 0.;
                Double_t e_bdt = N_sel[s][m] > 0 ? (Double_t)N_bdt[s][m]/N_sel[s][m] : 0.;
                Double_t e_tot = N_tot[s][m] > 0 ? (Double_t)N_bdt[s][m]/N_tot[s][m] : 0.;
                printf(" es=%.3f eb=%.3f et=%.4f  ", e_sel, e_bdt, e_tot);
            } else {
                printf(" %22s  ", "(model not run)");
            }
            if (m < N_MODELS - 1) printf("|");
        }
        printf("\n");
    }

    printf("%s", SEP);
    printf("  Ntot/k = total NanoAOD events [x1000] (score=-1 included in denominator)\n");
    printf("  Nsel/k = events with score >= 0 (passed trigger + kinematic selection)\n");
    printf("  Nbdt/k = events with score > threshold (passed BDT cut)\n");
    printf("  es = Nsel/Ntot (trigger+kin eff)  "
           "eb = Nbdt/Nsel (BDT eff | selection)  "
           "et = Nbdt/Ntot (total eff)\n");
    printf("%s\n", SEP);

    // Phase 4: one PDF per HLT model

    TString script_dir = gSystem->DirName(__FILE__);
    TString out_dir    = script_dir + "/output_getSignalEff_noMET";
    gSystem->mkdir(out_dir.Data(), kTRUE);

    printf("\n======== Phase 4: generating per-model PDFs ========\n");
    if (!SKIP_MASS)
        printf("  Mass pages ON (NanoAOD via xrootd -- run setproxy.sh first)\n");
    else
        printf("  Mass pages SKIPPED\n");

    for (Int_t m = 0; m < N_MODELS; ++m) {
        if (!run_model(m)) continue;
        // In plain test mode (no model explicitly selected) plot only the first model.
        if (RUN_TEST && ONLY_MODEL < 0 && m > 0) break;
        TString model_safe = MODELS[m].name;
        model_safe.ReplaceAll(" ", "_");
        // In test mode prefix the file name with "test_" (after the directory path).
        TString pdf_m = out_dir + "/" + (RUN_TEST ? "test_" : "")
                      + "getSignalEff_noMET_" + model_safe + ".pdf";
        printf("\n  -- Model %d/%d: %s -> %s --\n",
               m+1, N_MODELS, MODELS[m].name, pdf_m.Data());

        // Page 1: efficiency overview for this model
        printf("    Drawing page 1: efficiency overview...\n");
        {
            TCanvas* c1 = new TCanvas(Form("c_eff_ov_%d", m),
                                       Form("Eff overview - %s", MODELS[m].name),
                                       1800, 1400);
            c1->SetLogy(); c1->SetLogx();
            gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12);

            TH1D* fr = bkg_c[m] ? (TH1D*)bkg_c[m]->Clone(Form("fr_%d", m))
                                 : new TH1D(Form("fr_%d", m), "", NBINS, 0., 1.);
            fr->SetDirectory(nullptr);
            fr->GetXaxis()->SetRangeUser(1e-4, 1.);
            fr->GetXaxis()->SetTitle("BDT threshold");
            fr->GetYaxis()->SetTitle("Fraction of trigger-selected events above threshold");
            fr->GetXaxis()->SetTitleOffset(1.10);
            fr->GetYaxis()->SetTitleOffset(1.40);
            fr->SetMinimum(3e-5); fr->SetMaximum(1.3); fr->SetTitle("");
            fr->DrawClone("AXIS");

            TLatex hdr; hdr.SetNDC(); hdr.SetTextSize(0.030); hdr.SetTextAlign(31);
            hdr.DrawLatex(0.95, 0.93,
                Form("%s  |  2024  |  #sqrt{s} = 13.6 TeV  |  L_{int} = %.2f fb^{-1}",
                     MODELS[m].name, LUMI / 1e3));

            TLegend* leg = new TLegend(0.3, 0.2, 0.87, 0.32);
            leg->SetFillStyle(0); leg->SetBorderSize(0);
            leg->SetTextSize(0.032);
            leg->SetEntrySeparation(0.4);

            if (bkg_c[m]) {
                bkg_c[m]->SetLineColor(kBlack); bkg_c[m]->SetLineWidth(1);
                bkg_c[m]->DrawClone("hist same");
                leg->AddEntry(bkg_c[m],
                    Form("Background (thr=%.4f)", bdt_thr[m]), "l");
            } else {
                leg->AddEntry((TObject*)nullptr,
                    Form("thr=%.4f (FPR=1e-4)", bdt_thr[m]), "");
            }

            Color_t refcols[2] = { kBlue+1, kRed+1 };
            Style_t refstyles[2] = { 2, 3 };
            for (Int_t ref = 0; ref < 2; ++ref) {
                Int_t ridx = (ref == 0) ? 1 : 2;
                if (!sig_c[ridx][m]) continue;
                Int_t tbin = sig_c[ridx][m]->GetXaxis()->FindBin(bdt_thr[m]);
                Double_t eff = sig_c[ridx][m]->GetBinContent(tbin);
                sig_c[ridx][m]->SetLineColor(refcols[ref]);
                sig_c[ridx][m]->SetLineWidth(1);
                sig_c[ridx][m]->SetLineStyle(refstyles[ref]);
                sig_c[ridx][m]->DrawClone("hist same");
                leg->AddEntry(sig_c[ridx][m],
                    Form("%s  ( #varepsilon_{BDT}=%.3f)", SIG[ridx].label, eff), "l");
            }
            leg->Draw("same");

            c1->Print((pdf_m + "(").Data());
            delete fr; delete leg; delete c1;
        }

        // Pages 2+: mass distributions, one page per signal set (4 ctau cols each)
        if (!SKIP_MASS) {
            struct MassGrp { Int_t base; const char* label; Float_t mmax; };
            const MassGrp grps[N_SETS] = {
                {  0, "mpi=4,  mA=1.33 GeV", 3.0f },
                {  4, "mpi=4,  mA=0.40 GeV", 1.5f },
                {  8, "mpi=10, mA=1.00 GeV", 2.5f },
                { 12, "mpi=1,  mA=0.33 GeV", 1.0f },
                { 16, "mpi=10, mA=3.33 GeV", 7.0f },
            };
            // Last drawn page closes the PDF -> find the highest selected set.
            Int_t last_g = -1;
            for (Int_t g = 0; g < N_SETS; ++g) if (run_set(g)) last_g = g;
            for (Int_t g = 0; g < N_SETS; ++g) {
                if (!run_set(g)) continue;
                printf("    Drawing mass page (set %d): %s...\n", g, grps[g].label);
                make_mass_page(m, grps[g].base, grps[g].label,
                               grps[g].mmax, pdf_m, /*is_last=*/(g == last_g), bdt_thr[m]);
            }
        } else {
            // No mass pages: close the PDF that was opened on page 1
            TCanvas* ctmp = new TCanvas("ctmp","",10,10);
            ctmp->Print((pdf_m + ")").Data());
            delete ctmp;
        }

        printf("  -> %s\n", pdf_m.Data());
    }

    printf("\nPDFs saved (one per HLT model) in %s\n", out_dir.Data());

    // Phase 5: efficiency tables PDF
    if (MAKE_TABLES)
        make_eff_tables_pdf(out_dir, bdt_thr, N_tot, N_sel, N_bdt);
}


// ROOT macro entry point - name must match the filename.
//
// Optional argument (space-separated keywords, combinable):
//   "test"        - 1 file per sample, event loop capped at 5000 events
//   "compute_thr" - recompute thresholds from QCD data (slow; hardcoded by default).
//                   QCD events are weighted by bin cross section sigma_q (matches the
//                   "Total (weighted)" column), so the threshold is set on the
//                   xs-weighted QCD event count at FPR = 1e-4.
//   "skip_mass"   - omit mass distribution pages (default is to include them)
//   "tables"      - produce LaTeX tables PDFs + CSVs (requires setproxy.sh for NanoAOD)
//   "model0/1/2"  - process only one HLT model (0=Mu10ORDoubleMu, 1=Mu10, 2=DoubleMu)
//   "set0..set4"  - process only one signal mass set (4 ctau cols each):
//                   0=mpi4/mA1.33 1=mpi4/mA0.40 2=mpi10/mA1.00 3=mpi1/mA0.33 4=mpi10/mA3.33
//
// Examples:
//   root -l -b -q getSignalEff_noMET.C                                  # full run
//   root -l -b -q 'getSignalEff_noMET.C("test")'                        # fast debug
//   root -l -b -q 'getSignalEff_noMET.C("skip_mass")'                   # efficiency plots only
//   root -l -b -q 'getSignalEff_noMET.C("skip_mass test")'              # efficiency only, fast
//   root -l -b -q 'getSignalEff_noMET.C("compute_thr")'                 # recompute thresholds
//   root -l -b -q 'getSignalEff_noMET.C("tables skip_mass")'            # tables PDFs+CSVs (no mass pages)
//   root -l -b -q 'getSignalEff_noMET.C("tables skip_mass test")'       # tables PDFs+CSVs, fast debug
//   root -l -b -q 'getSignalEff_noMET.C("compute_thr model1 skip_mass")'# recompute thr for Mu10 only
// All output goes to _tools/output_getSignalEff_noMET/
int getSignalEff_noMET(TString mode = "")
{
    mode.ToLower();
    if (mode.Contains("test")) {
        RUN_TEST = kTRUE;
        printf("\n*** TEST MODE: 1 file per sample, event loop capped at %lld events ***\n",
               TEST_EVTS);
    }
    if (mode.Contains("compute_thr")) {
        COMPUTE_THR = kTRUE;
        printf("*** COMPUTE_THR: thresholds recomputed from QCD data ***\n");
    }
    if (mode.Contains("skip_mass")) {
        SKIP_MASS = kTRUE;
        printf("*** SKIP_MASS: mass distribution pages will not be generated ***\n");
    }
    if (mode.Contains("tables")) {
        MAKE_TABLES = kTRUE;
        printf("*** MAKE_TABLES: efficiency tables PDF will be generated ***\n");
    }
    // Single-model selector: process only one of the 3 HLT models.
    //   model0 -> Mu10ORDoubleMu, model1 -> Mu10, model2 -> DoubleMu
    if      (mode.Contains("model0")) ONLY_MODEL = 0;
    else if (mode.Contains("model1")) ONLY_MODEL = 1;
    else if (mode.Contains("model2")) ONLY_MODEL = 2;
    if (ONLY_MODEL >= 0)
        printf("*** ONLY_MODEL: processing model %d only (%s) ***\n",
               ONLY_MODEL, MODELS[ONLY_MODEL].name);
    // Single signal-set selector: process only one of the 5 mass sets (4 ctau each).
    //   set0 mpi=4/mA=1.33, set1 mpi=4/mA=0.40, set2 mpi=10/mA=1.00,
    //   set3 mpi=1/mA=0.33, set4 mpi=10/mA=3.33
    if      (mode.Contains("set0")) ONLY_SET = 0;
    else if (mode.Contains("set1")) ONLY_SET = 1;
    else if (mode.Contains("set2")) ONLY_SET = 2;
    else if (mode.Contains("set3")) ONLY_SET = 3;
    else if (mode.Contains("set4")) ONLY_SET = 4;
    if (ONLY_SET >= 0)
        printf("*** ONLY_SET: processing signal set %d only (%s) ***\n",
               ONLY_SET, SIG[ONLY_SET*4].label);
    printf("\n");
    run_analysis();
    return 0;
}
