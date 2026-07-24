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
//   Each cutflow page also carries a 95% CL upper-limit-on-B(H->psipsibar) table
//   and an Asimov-significance table at FPR = 1e-1..1e-4 (ROOT-drawn companions below).
//
// Limit PDFs ("tables" keyword required), one per model -- ROOT-drawn:
//   limits_<model>.pdf
//     page 1     - fraction of SV-selected QCD above threshold vs BDT threshold,
//                  with dashed lines marking the FPR = 1e-1..1e-4 working points
//     pages 2..6 - one per signal set: 95% CL upper limit on B(H->psipsibar) vs ctau,
//                  one coloured/markered curve per FPR working point
//   significance_<model>.pdf
//     one page per signal set: asymptotic Asimov significance Z vs ctau,
//     one coloured/markered curve per FPR working point (legends at top)
//   limits_by_lxy_<model>.pdf / significance_by_lxy_<model>.pdf
//     page 1 (limits only) - inclusive FPR calibration; then one page per signal
//     set with 3 pads (one per lxy category), each = limit/Z vs ctau for that bin
//     (each lxy bin is its own search region: own windowed-SV thresholds+weight)
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
#include "TGraph.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include "TChainElement.h"
#include <iostream>
#include <cmath>
#include <vector>

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

static const char*    BRANCH   = "xgb01_NOJETS";
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
      // "modeltag__scenarioA_all_no_DA_old_BDT_with_dRmSV_2024" },
      "modeltag__scenarioA_all_no_DA_old_BDT_2024" },
    { "Mu10",
      // "modeltag__scenarioA_all_no_DA_old_BDT_with_dRmSV_2024_Mu10" },
      "modeltag__scenarioA_all_no_DA_old_BDT_2024_Mu10" },
    { "DoubleMu",
      // "modeltag__scenarioA_all_no_DA_old_BDT_with_dRmSV_2024_DoubleMu" },
      "modeltag__scenarioA_all_no_DA_old_BDT_2024_DoubleMu" },
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

// B(A'->mumu) per signal set, extracted from the Pythia gen fragments
//   .../EXO-MCsampleRequests/genFragments/Hadronizer/13p6TeV/DQCD_Run3/
//   GluGluHToDarkShowers-ScenarioA_Par-*_cfi.py
// (normalised from the 9900015:addChannel A' decay weights). B(pi3->A'A')=1 for
// every set (4900111:addChannel = 1 1.0 91 9900015 9900015). These BRs are already
// folded into the signal efficiency by the MC; the array is for labelling/reference.
//   set0 mpi=4/mA=1.33  set1 mpi=4/mA=0.40  set2 mpi=10/mA=1.00
//   set3 mpi=1/mA=0.33  set4 mpi=10/mA=3.33
static const Double_t BR_A_MUMU[N_SETS] = { 0.317, 0.440, 0.307, 0.464, 0.193 };

// ---------------------------------------------------------------------------
// Mass-window background (prepared 2026-06-24)
//
// Rather than counting the whole QCD sample as background, restrict it to a
// dimuon-mass window around each signal's A' peak (m_A). The window is applied
// to the leading-muonSV invariant mass (lead_mass) on BOTH signal and QCD, so
// the SV-selected background yield becomes per-signal-set.
//
// Peaks (m_A) per set, matching the SIG[] set ordering:
//   set0 m_A=1.33  set1 m_A=0.40  set2 m_A=1.00  set3 m_A=0.33  set4 m_A=3.33
// Window: 0.9*m_A to 1.1*m_A (i.e. +-10% around the resonance peak), same
// fractional width for every set.
static const Double_t MA_PEAK   [N_SETS] = { 1.33, 0.40, 1.00, 0.33, 3.33 };
static const Double_t MWIN_FRAC_LO = 0.9;   // window lower edge = MWIN_FRAC_LO * m_A
static const Double_t MWIN_FRAC_HI = 1.1;   // window upper edge = MWIN_FRAC_HI * m_A
// Master switch: when false the limit/significance use the full QCD (inclusive);
// when true they use the per-set 0.9-1.1*m_A windowed-SV QCD background.
static const Bool_t   USE_MASS_WINDOW = true;
// When true (and USE_MASS_WINDOW), the SAME window is also applied to the signal
// (eps_S measured in the window). Set false to window ONLY the QCD background.
static const Bool_t   WINDOW_SIGNAL   = true;

// Per-set leading-muonSV mass-window cut string (for RDataFrame Filter):
//   0.9*m_A <= lead_mass < 1.1*m_A
// 'set' is the signal-set index (= signal index / N_CTAU). Requires a
// 'lead_mass' column to have been Define()'d on the node.
static inline TString mass_window_cut(Int_t set) {
    return TString::Format("lead_mass >= %.4ff && lead_mass < %.4ff",
                           MWIN_FRAC_LO * MA_PEAK[set],
                           MWIN_FRAC_HI * MA_PEAK[set]);
}

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

#include "getSignalEff_cache.h"
// Caches filled during the read phase and written to a ROOT file at the end.
static CutflowCache g_cut  = {};
static MassCache    g_mass = {};
static TablesCache  g_tab  = {};

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

    // EXTRACT: no canvas -- histograms/counts are stored into g_mass below.
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

            // ---- EXTRACT: store mass histograms + counts into the cache ----
            (void)pad; (void)e_sel; (void)e_bdt; (void)e_tot;
            g_mass.h_sel[m_idx][s][row] = h_sel;
            g_mass.h_bdt[m_idx][s][row] = h_bdt;
            g_mass.n_lxy[m_idx][s][row] = n_lxy;
            g_mass.n_sel[m_idx][s][row] = n_sel;
            g_mass.n_bdt[m_idx][s][row] = n_bdt;
            g_mass.has [m_idx][s][row] = 1;
        }

        delete nano_ch; delete bdt_ch;
    }
    (void)group_label; (void)mass_max; (void)pdf; (void)is_last;
    printf("    mass data collected for %s\n", MODELS[m_idx].name);
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
    Long64_t N_sv_s  [N_SIG][N_MODELS]    = {};  // presel + SV quality (chi2<10)
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
    TH1D*    sig_sv_cumul_win[N_SIG][N_MODELS] = {};  // same, SV+mass-window selected (WINDOW_SIGNAL)
    TH1D*    h_bkg_sv    [N_MODELS]        = {};  // xs/N_tot-weighted SV-selected QCD score
    TH1D*    bkg_sv_c    [N_MODELS]        = {};  // its normalised reverse-cumulative (= FPR vs threshold)
    Double_t sig_thr     [N_MODELS][4]     = {};  // BDT threshold at each FPR target, per model
    // Mass-window background (USE_MASS_WINDOW): per-(model,set) windowed-SV QCD.
    // When the window is on, the limit/significance use these instead of the
    // per-model arrays above (which stay as the inclusive cutflow reference).
    TH1D*    bkg_sv_c_win[N_MODELS][N_SETS] = {};  // FPR-vs-threshold on windowed-SV QCD
    Double_t sig_thr_win [N_MODELS][N_SETS][4] = {};  // BDT threshold at each FPR, per (model,set)
    Double_t w_sv_win    [N_MODELS][N_SETS] = {};  // windowed SV bkg weight sum_q xs_q*eps_q^(SV+win)
    // Per-lxy-category versions (for the "by_lxy" limit/significance PDFs): each
    // lxy bin is its own search region -> own windowed-SV thresholds + weight, and
    // signal efficiency measured in that lxy bin.
    TH1D*    sig_sv_cumul_lxy[N_SIG][N_MODELS][3] = {};   // windowed-SV signal score, per lxy
    Double_t sig_thr_win_lxy [N_MODELS][N_SETS][3][4] = {};  // FPR thresholds per (model,set,lxy)
    Double_t w_sv_win_lxy    [N_MODELS][N_SETS][3]    = {};  // windowed-SV bkg weight per (model,set,lxy)
    // RAW (unweighted) QCD count reverse-cumulatives -> #events above threshold,
    // for the event-count tables / QCD-lxy plots / CSVs (absolute MC counts).
    TH1D*    bkg_raw_c_win    [N_MODELS][N_SETS]    = {};  // windowed-SV QCD raw count cumul
    TH1D*    bkg_raw_c_win_lxy[N_MODELS][N_SETS][3] = {};  // windowed-SV QCD raw count cumul, per lxy
    // xs-weighted (by xs_q, NOT xs_q/N_tot) QCD count reverse-cumulatives. Dividing the
    // value above threshold by xs_total gives the sigma-weighted AVERAGE MC count passing
    // = same definition as the cutflow "QCD Total (weighted)" column (no luminosity).
    TH1D*    bkg_xsw_c_win    [N_MODELS][N_SETS]    = {};  // windowed-SV QCD xs-weighted count cumul
    TH1D*    bkg_xsw_c_win_lxy[N_MODELS][N_SETS][3] = {};  // windowed-SV QCD xs-weighted count cumul, per lxy

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
    // SV quality: at least one muonSV with chi2 < 10
    // (dR(mu1,mu2) < 1.2 cut removed 2026-06-24)
    auto f_sv = [](const ROOT::RVec<float>& chi2) -> bool {
        for (int i = 0; i < (int)chi2.size(); ++i) {
            if (chi2[i] < 10.f) return true;
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
                .Define("lead_dxy", pick_fn, {"lead_idx","muonSV_dxy"})
                .Define("lead_mass",pick_fn, {"lead_idx","muonSV_mass"});
            ROOT::RDF::RNode df2 = [&]() -> ROOT::RDF::RNode {
                if (m==0) return df_lxy.Define("kin_pass", f_kin_or,
                    {"muonSV_mu1pt","muonSV_mu1eta","muonSV_mu2pt","muonSV_mu2eta"});
                if (m==1) return df_lxy.Define("kin_pass", f_kin_mu10,
                    {"muonSV_mu1pt","muonSV_mu1eta","muonSV_mu2pt","muonSV_mu2eta"});
                return df_lxy.Define("kin_pass", f_kin_dmu,
                    {"muonSV_mu1pt","muonSV_mu1eta","muonSV_mu2pt","muonSV_mu2eta"});
            }();
            ROOT::RDF::RNode df3 = df2.Define("sv_pass", f_sv, {"muonSV_chi2"});
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
            // SV+mass-window-selected signal score (windowed eps_S; set = s/N_CTAU)
            TString fsw = fsv + " && " + mass_window_cut(s / N_CTAU);
            auto h_ssv_w = df3.Filter(fsw.Data())
                            .Histo1D({Form("h_ssvw_%d_%d",s,m),"",NBINS,0.,1.}, BRANCH);
            // per-lxy windowed-SV signal score (for the by-lxy plots); same window
            // semantics as eps_S (window only if USE_MASS_WINDOW && WINDOW_SIGNAL)
            ROOT::RDF::RResultPtr<TH1D> h_ssv_lxy[3];
            for (Int_t c=0;c<3;++c) {
                TString fsl = fsv;
                if (USE_MASS_WINDOW && WINDOW_SIGNAL) fsl += " && " + mass_window_cut(s / N_CTAU);
                fsl += TString::Format(" && lead_dxy >= %.1ff && lead_dxy < %.1ff", LXY_LO[c], LXY_HI[c]);
                h_ssv_lxy[c] = df3.Filter(fsl.Data())
                    .Histo1D({Form("h_ssvlxy_%d_%d_%d",s,m,c),"",NBINS,0.,1.}, BRANCH);
            }
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
            // windowed-SV signal reverse-cumulative (for WINDOW_SIGNAL eps_S)
            TH1D* hssvw = (TH1D*)h_ssv_w->Clone(Form("hssvw_%d_%d",s,m));
            hssvw->SetDirectory(nullptr);
            sig_sv_cumul_win[s][m] = (TH1D*)hssvw->GetCumulative(kFALSE);
            sig_sv_cumul_win[s][m]->SetDirectory(nullptr);
            delete hssvw;
            // per-lxy windowed-SV signal reverse-cumulatives (by-lxy eps_S)
            for (Int_t c=0;c<3;++c) {
                TH1D* hl = (TH1D*)h_ssv_lxy[c]->Clone(Form("hssvlxy_%d_%d_%d",s,m,c));
                hl->SetDirectory(nullptr);
                sig_sv_cumul_lxy[s][m][c] = (TH1D*)hl->GetCumulative(kFALSE);
                sig_sv_cumul_lxy[s][m][c]->SetDirectory(nullptr);
                delete hl;
            }
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
        // per-set windowed-SV background score accumulators (mass-window option)
        TH1D* hbw[N_SETS] = {};
        for (Int_t s=0; s<N_SETS; ++s) {
            hbw[s] = new TH1D(Form("h_bkg_svwin_%d_%d",m,s),"",NBINS,0.,1.);
            hbw[s]->SetDirectory(nullptr);
        }
        // per-(set,lxy) windowed-SV background accumulators (by-lxy plots)
        TH1D* hbwl[N_SETS][3] = {};
        for (Int_t s=0; s<N_SETS; ++s) for (Int_t c=0;c<3;++c) {
            hbwl[s][c] = new TH1D(Form("h_bkg_svwinlxy_%d_%d_%d",m,s,c),"",NBINS,0.,1.);
            hbwl[s][c]->SetDirectory(nullptr);
        }
        // RAW (unweighted) windowed-SV QCD count accumulators (for absolute counts)
        TH1D* hbw_raw[N_SETS] = {}; TH1D* hbwl_raw[N_SETS][3] = {};
        for (Int_t s=0; s<N_SETS; ++s) {
            hbw_raw[s] = new TH1D(Form("h_bkg_rawwin_%d_%d",m,s),"",NBINS,0.,1.);
            hbw_raw[s]->SetDirectory(nullptr);
            for (Int_t c=0;c<3;++c) {
                hbwl_raw[s][c] = new TH1D(Form("h_bkg_rawwinlxy_%d_%d_%d",m,s,c),"",NBINS,0.,1.);
                hbwl_raw[s][c]->SetDirectory(nullptr);
            }
        }
        // xs-weighted (by xs_q) windowed-SV QCD count accumulators (sigma-weighted avg count)
        TH1D* hbw_xs[N_SETS] = {}; TH1D* hbwl_xs[N_SETS][3] = {};
        for (Int_t s=0; s<N_SETS; ++s) {
            hbw_xs[s] = new TH1D(Form("h_bkg_xswin_%d_%d",m,s),"",NBINS,0.,1.);
            hbw_xs[s]->SetDirectory(nullptr);
            for (Int_t c=0;c<3;++c) {
                hbwl_xs[s][c] = new TH1D(Form("h_bkg_xswinlxy_%d_%d_%d",m,s,c),"",NBINS,0.,1.);
                hbwl_xs[s][c]->SetDirectory(nullptr);
            }
        }
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
            ROOT::RDF::RNode df = df_kin.Define("sv_pass", f_sv, {"muonSV_chi2"})
                .Define("lead_idx",  lead_fn, {"muonSV_chi2","muonSV_dlen"})
                .Define("lead_mass", pick_fn, {"lead_idx","muonSV_mass"})
                .Define("lead_dxy",  pick_fn, {"lead_idx","muonSV_dxy"});
            TString fpre = TString::Format("%s >= 0.", BRANCH);
            TString fsv  = TString::Format("sv_pass && %s >= 0.", BRANCH);
            TString fbdt = TString::Format("%s > %.8f", BRANCH, bdt_thr[m]);
            auto c_tot = df.Count();
            auto c_kin = df.Filter("kin_pass").Count();
            auto c_pre = df.Filter(fpre.Data()).Count();
            auto c_sv  = df.Filter(fsv.Data()).Count();
            auto c_bdt = df.Filter(fbdt.Data()).Count();
            // Pre-filter SV+score>=0 ONCE (single jitted filter); the per-(set,lxy)
            // windowed histograms below then use *compiled* lambda filters (much
            // faster than ~20 jitted string filters per QCD bin).
            ROOT::RDF::RNode df_sv = df.Filter(fsv.Data());
            auto h_qsv = df_sv.Histo1D({Form("h_qsv_%d_%d",m,q),"",NBINS,0.,1.}, BRANCH);
            // per-set windowed-SV score histograms (compiled mass-window filter)
            ROOT::RDF::RResultPtr<TH1D> h_qsv_w[N_SETS];
            for (Int_t s=0; s<N_SETS; ++s) {
                double mlo=MWIN_FRAC_LO*MA_PEAK[s], mhi=MWIN_FRAC_HI*MA_PEAK[s];
                h_qsv_w[s] = df_sv.Filter([mlo,mhi](float lm){ return lm>=mlo && lm<mhi; },
                                          {"lead_mass"})
                    .Histo1D({Form("h_qsvw_%d_%d_%d",m,q,s),"",NBINS,0.,1.}, BRANCH);
            }
            // per-(set,lxy) windowed-SV score histograms (compiled mass+lxy filter)
            ROOT::RDF::RResultPtr<TH1D> h_qsv_wl[N_SETS][3];
            for (Int_t s=0; s<N_SETS; ++s) for (Int_t c=0;c<3;++c) {
                double mlo=MWIN_FRAC_LO*MA_PEAK[s], mhi=MWIN_FRAC_HI*MA_PEAK[s];
                float  llo=LXY_LO[c], lhi=LXY_HI[c];
                h_qsv_wl[s][c] = df_sv.Filter(
                    [mlo,mhi,llo,lhi](float lm, float ld){ return lm>=mlo && lm<mhi && ld>=llo && ld<lhi; },
                    {"lead_mass","lead_dxy"})
                    .Histo1D({Form("h_qsvwl_%d_%d_%d_%d",m,q,s,c),"",NBINS,0.,1.}, BRANCH);
            }
            N_tot_q[q][m] = (Long64_t)*c_tot;
            N_kin_q[q][m] = (Long64_t)*c_kin;
            N_pre_q[q][m] = (Long64_t)*c_pre;
            N_sv_q [q][m] = (Long64_t)*c_sv;
            N_bdt_q[q][m] = (Long64_t)*c_bdt;
            // accumulate luminosity-yield-weighted (xs/N_tot) SV-selected score
            if (N_tot_q[q][m] > 0) {
                Double_t wnorm = QCD[q].xs / (Double_t)N_tot_q[q][m];
                TH1D* hq = (TH1D*)h_qsv->Clone(Form("hq_%d_%d",m,q));
                hq->SetDirectory(nullptr); hq->Scale(wnorm);
                h_bkg_sv[m]->Add(hq); delete hq;
                // per-set windowed SV background: score shape + weight sum_q xs*N_win/N_tot
                for (Int_t s=0; s<N_SETS; ++s) {
                    TH1D* hqw = (TH1D*)h_qsv_w[s]->Clone(Form("hqw_%d_%d_%d",m,q,s));
                    hqw->SetDirectory(nullptr);
                    w_sv_win[m][s] += QCD[q].xs * hqw->GetEntries() / (Double_t)N_tot_q[q][m];
                    hqw->Scale(wnorm);
                    hbw[s]->Add(hqw); delete hqw;
                    hbw_raw[s]->Add(h_qsv_w[s].GetPtr());   // raw (unweighted) count
                    hbw_xs[s]->Add(h_qsv_w[s].GetPtr(), QCD[q].xs); // xs-weighted (by xs_q) count
                    // per-(set,lxy)
                    for (Int_t c=0;c<3;++c) {
                        TH1D* hqwl = (TH1D*)h_qsv_wl[s][c]->Clone(Form("hqwl_%d_%d_%d_%d",m,q,s,c));
                        hqwl->SetDirectory(nullptr);
                        w_sv_win_lxy[m][s][c] += QCD[q].xs * hqwl->GetEntries() / (Double_t)N_tot_q[q][m];
                        hqwl->Scale(wnorm);
                        hbwl[s][c]->Add(hqwl); delete hqwl;
                        hbwl_raw[s][c]->Add(h_qsv_wl[s][c].GetPtr());  // raw count
                        hbwl_xs[s][c]->Add(h_qsv_wl[s][c].GetPtr(), QCD[q].xs); // xs-weighted count
                    }
                }
            }
            delete nano_ch; delete bdt_ch;
        }
        // inclusive background score shape -> reverse cumulative -> FPR thresholds
        if (h_bkg_sv[m]->Integral() > 0.) {
            normalise(h_bkg_sv[m]);
            bkg_sv_c[m] = (TH1D*)h_bkg_sv[m]->GetCumulative(kFALSE);  // kept for the limit PDF
            bkg_sv_c[m]->SetDirectory(nullptr);
            for (Int_t k=0;k<4;++k) { Double_t dum; sig_thr[m][k]=find_thr(bkg_sv_c[m],FPR_TGT[k],dum); }
            printf("  [signif] thresholds (FPR 1e-1..1e-4): %.4f %.4f %.4f %.4f\n",
                   sig_thr[m][0],sig_thr[m][1],sig_thr[m][2],sig_thr[m][3]);
        } else {
            for (Int_t k=0;k<4;++k) sig_thr[m][k]=1.;
        }
        // per-set windowed background shape -> reverse cumulative -> FPR thresholds
        for (Int_t s=0; s<N_SETS; ++s) {
            if (hbw[s]->Integral() > 0.) {
                normalise(hbw[s]);
                bkg_sv_c_win[m][s] = (TH1D*)hbw[s]->GetCumulative(kFALSE);
                bkg_sv_c_win[m][s]->SetDirectory(nullptr);
                for (Int_t k=0;k<4;++k) { Double_t dum; sig_thr_win[m][s][k]=find_thr(bkg_sv_c_win[m][s],FPR_TGT[k],dum); }
            } else {
                for (Int_t k=0;k<4;++k) sig_thr_win[m][s][k]=1.;
            }
            delete hbw[s];
            // RAW (unnormalized) windowed-SV QCD count reverse-cumulative
            bkg_raw_c_win[m][s] = (TH1D*)hbw_raw[s]->GetCumulative(kFALSE);
            bkg_raw_c_win[m][s]->SetDirectory(nullptr); delete hbw_raw[s];
            // xs-weighted (by xs_q) windowed-SV QCD count reverse-cumulative
            bkg_xsw_c_win[m][s] = (TH1D*)hbw_xs[s]->GetCumulative(kFALSE);
            bkg_xsw_c_win[m][s]->SetDirectory(nullptr); delete hbw_xs[s];
            // per-(set,lxy) windowed background -> reverse cumulative -> FPR thresholds
            for (Int_t c=0;c<3;++c) {
                if (hbwl[s][c]->Integral() > 0.) {
                    normalise(hbwl[s][c]);
                    TH1D* cwl = (TH1D*)hbwl[s][c]->GetCumulative(kFALSE); cwl->SetDirectory(nullptr);
                    for (Int_t k=0;k<4;++k) { Double_t dum; sig_thr_win_lxy[m][s][c][k]=find_thr(cwl,FPR_TGT[k],dum); }
                    delete cwl;
                } else {
                    for (Int_t k=0;k<4;++k) sig_thr_win_lxy[m][s][c][k]=1.;
                }
                delete hbwl[s][c];
                // RAW windowed+lxy QCD count reverse-cumulative
                bkg_raw_c_win_lxy[m][s][c] = (TH1D*)hbwl_raw[s][c]->GetCumulative(kFALSE);
                bkg_raw_c_win_lxy[m][s][c]->SetDirectory(nullptr); delete hbwl_raw[s][c];
                // xs-weighted windowed+lxy QCD count reverse-cumulative
                bkg_xsw_c_win_lxy[m][s][c] = (TH1D*)hbwl_xs[s][c]->GetCumulative(kFALSE);
                bkg_xsw_c_win_lxy[m][s][c]->SetDirectory(nullptr); delete hbwl_xs[s][c];
            }
        }
    }

    // ---- EXTRACT: copy computed tables data into the global cache struct ----
    memcpy(g_tab.N_kin_s,   N_kin_s,   sizeof(N_kin_s));
    memcpy(g_tab.N_sv_s,    N_sv_s,    sizeof(N_sv_s));
    memcpy(g_tab.N_lxy,     N_lxy,     sizeof(N_lxy));
    memcpy(g_tab.N_k_lxy,   N_k_lxy,   sizeof(N_k_lxy));
    memcpy(g_tab.N_p_lxy,   N_p_lxy,   sizeof(N_p_lxy));
    memcpy(g_tab.N_sv_lxy,  N_sv_lxy,  sizeof(N_sv_lxy));
    memcpy(g_tab.N_b_lxy,   N_b_lxy,   sizeof(N_b_lxy));
    memcpy(g_tab.N_tot_q,   N_tot_q,   sizeof(N_tot_q));
    memcpy(g_tab.N_kin_q,   N_kin_q,   sizeof(N_kin_q));
    memcpy(g_tab.N_pre_q,   N_pre_q,   sizeof(N_pre_q));
    memcpy(g_tab.N_sv_q,    N_sv_q,    sizeof(N_sv_q));
    memcpy(g_tab.N_bdt_q,   N_bdt_q,   sizeof(N_bdt_q));
    memcpy(g_tab.sig_thr,          sig_thr,          sizeof(sig_thr));
    memcpy(g_tab.sig_thr_win,      sig_thr_win,      sizeof(sig_thr_win));
    memcpy(g_tab.sig_thr_win_lxy,  sig_thr_win_lxy,  sizeof(sig_thr_win_lxy));
    memcpy(g_tab.w_sv_win,         w_sv_win,         sizeof(w_sv_win));
    memcpy(g_tab.w_sv_win_lxy,     w_sv_win_lxy,     sizeof(w_sv_win_lxy));
    memcpy(g_tab.sig_sv_cumul,     sig_sv_cumul,     sizeof(sig_sv_cumul));
    memcpy(g_tab.sig_sv_cumul_win, sig_sv_cumul_win, sizeof(sig_sv_cumul_win));
    memcpy(g_tab.sig_sv_cumul_lxy, sig_sv_cumul_lxy, sizeof(sig_sv_cumul_lxy));
    memcpy(g_tab.bkg_sv_c,         bkg_sv_c,         sizeof(bkg_sv_c));
    memcpy(g_tab.bkg_sv_c_win,     bkg_sv_c_win,     sizeof(bkg_sv_c_win));
    memcpy(g_tab.bkg_raw_c_win,    bkg_raw_c_win,    sizeof(bkg_raw_c_win));
    memcpy(g_tab.bkg_raw_c_win_lxy,bkg_raw_c_win_lxy,sizeof(bkg_raw_c_win_lxy));
    memcpy(g_tab.bkg_xsw_c_win,    bkg_xsw_c_win,    sizeof(bkg_xsw_c_win));
    memcpy(g_tab.bkg_xsw_c_win_lxy,bkg_xsw_c_win_lxy,sizeof(bkg_xsw_c_win_lxy));
    printf("  [extract] tables/significance data collected\n");
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

        // EXTRACT: cutflow-overview is not drawn (bkg_c/sig_c are cached instead).

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
        } // (SKIP_MASS: nothing to collect)

        (void)pdf_m;
    }

    // Phase 5: collect tables/significance data (fills g_tab)
    if (MAKE_TABLES)
        make_eff_tables_pdf(out_dir, bdt_thr, N_tot, N_sel, N_bdt);

    // ---- EXTRACT: fill the cutflow cache and write everything to disk ----
    memcpy(g_cut.bdt_thr, bdt_thr, sizeof(bdt_thr));
    memcpy(g_cut.bkg_eff, bkg_eff, sizeof(bkg_eff));
    memcpy(g_cut.N_tot,   N_tot,   sizeof(N_tot));
    memcpy(g_cut.N_sel,   N_sel,   sizeof(N_sel));
    memcpy(g_cut.N_bdt,   N_bdt,   sizeof(N_bdt));
    memcpy(g_cut.bkg_c,   bkg_c,   sizeof(bkg_c));
    memcpy(g_cut.sig_c,   sig_c,   sizeof(sig_c));
    for (Int_t m=0;m<N_MODELS;++m) g_cut.has_bkg_c[m] = bkg_c[m] ? 1 : 0;
    for (Int_t s=0;s<N_SIG;++s) for (Int_t m=0;m<N_MODELS;++m)
        g_cut.has_sig_c[s][m] = sig_c[s][m] ? 1 : 0;

    TString cache = gse_cache_path(out_dir, RUN_TEST, ONLY_MODEL, ONLY_SET);
    TFile* fout = TFile::Open(cache, "RECREATE");
    gse_save_cutflow(fout, g_cut);
    gse_save_mass   (fout, g_mass);
    gse_save_tables (fout, g_tab);
    fout->Close(); delete fout;
    printf("\n[extract] wrote cache -> %s\n", cache.Data());
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
int getSignalEff_extract(TString mode = "")
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
