// ===========================================================================
//  getSignalEff_cache.h  --  serialization contract shared by
//    getSignalEff_extract.C  (reads data  -> writes cache)
//    getSignalEff_render.C   (reads cache -> makes plots/tables)
//
//  DO NOT include before the config constants (N_SIG, N_MODELS, N_QCD, N_SETS)
//  are defined -- this header uses them for the array dimensions.  It is meant
//  to be #included inside the two driver macros, which are copies of
//  getSignalEff_noMET.C and therefore already define those constants.
//
//  The cache is a single ROOT file holding:
//    - every histogram (cutflow cumulatives, signal/QCD score cumulatives,
//      mass-distribution histograms) written under a deterministic name;
//    - one TTree "scalars" with a single entry whose fixed-size array branches
//      carry all the count / threshold / weight arrays.
//
//  getSignalEff_noMET.C is NOT touched by any of this.
// ===========================================================================
#ifndef GETSIGNALEFF_CACHE_H
#define GETSIGNALEFF_CACHE_H

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TString.h"
#include <cstdio>

// Deterministic cache-file path (extract writes it, render reads it). Both
// derive the same name from the run mode so partial (test / single-model /
// single-set) extracts don't collide and render picks the matching cache.
static inline TString gse_cache_path(const TString& out_dir, Bool_t run_test,
                                     Int_t only_model, Int_t only_set)
{
    TString n = out_dir + "/" + (run_test ? "test_" : "") + "getSignalEff_cache";
    if (only_model >= 0) n += Form("_model%d", only_model);
    if (only_set   >= 0) n += Form("_set%d",   only_set);
    return n + ".root";
}

// --- low-level histogram helpers -------------------------------------------

// Write a histogram under 'name' (skipped if null). A companion 0/1 flag branch
// in the scalar tree records presence so the reader knows whether to expect it.
static inline void gse_writeH(TFile* f, const TH1D* h, const char* name)
{
    if (!h) return;
    f->cd();
    TH1D* c = (TH1D*)h->Clone(name);
    c->SetDirectory(f);
    c->Write(name, TObject::kOverwrite);
}

// Read a histogram by name; returns a detached clone (SetDirectory(0)) or null.
static inline TH1D* gse_readH(TFile* f, const char* name)
{
    TH1D* h = (TH1D*)f->Get(name);
    if (!h) return nullptr;
    TH1D* c = (TH1D*)h->Clone(name);
    c->SetDirectory(nullptr);
    return c;
}

// ===========================================================================
//  CUTFLOW cache  (produced in run_analysis Phase 1/2)
// ===========================================================================
struct CutflowCache {
    Double_t bdt_thr[N_MODELS];
    Double_t bkg_eff[N_MODELS];
    Long64_t N_tot[N_SIG][N_MODELS];
    Long64_t N_sel[N_SIG][N_MODELS];
    Long64_t N_bdt[N_SIG][N_MODELS];
    Int_t    has_bkg_c[N_MODELS];          // presence flags
    Int_t    has_sig_c[N_SIG][N_MODELS];
    TH1D*    bkg_c[N_MODELS];
    TH1D*    sig_c[N_SIG][N_MODELS];
};

static inline void gse_save_cutflow(TFile* f, CutflowCache& C)
{
    f->cd();
    TTree* t = new TTree("cutflow_scalars", "");
    t->Branch("bdt_thr",   C.bdt_thr,   Form("bdt_thr[%d]/D",   N_MODELS));
    t->Branch("bkg_eff",   C.bkg_eff,   Form("bkg_eff[%d]/D",   N_MODELS));
    t->Branch("N_tot",     C.N_tot,     Form("N_tot[%d][%d]/L", N_SIG, N_MODELS));
    t->Branch("N_sel",     C.N_sel,     Form("N_sel[%d][%d]/L", N_SIG, N_MODELS));
    t->Branch("N_bdt",     C.N_bdt,     Form("N_bdt[%d][%d]/L", N_SIG, N_MODELS));
    t->Branch("has_bkg_c", C.has_bkg_c, Form("has_bkg_c[%d]/I", N_MODELS));
    t->Branch("has_sig_c", C.has_sig_c, Form("has_sig_c[%d][%d]/I", N_SIG, N_MODELS));
    t->Fill();
    t->Write("cutflow_scalars", TObject::kOverwrite);
    for (Int_t m=0;m<N_MODELS;++m) gse_writeH(f, C.bkg_c[m], Form("bkg_c_%d",m));
    for (Int_t s=0;s<N_SIG;++s) for (Int_t m=0;m<N_MODELS;++m)
        gse_writeH(f, C.sig_c[s][m], Form("sig_c_%d_%d",s,m));
}

static inline void gse_load_cutflow(TFile* f, CutflowCache& C)
{
    TTree* t = (TTree*)f->Get("cutflow_scalars");
    t->SetBranchAddress("bdt_thr",   C.bdt_thr);
    t->SetBranchAddress("bkg_eff",   C.bkg_eff);
    t->SetBranchAddress("N_tot",     C.N_tot);
    t->SetBranchAddress("N_sel",     C.N_sel);
    t->SetBranchAddress("N_bdt",     C.N_bdt);
    t->SetBranchAddress("has_bkg_c", C.has_bkg_c);
    t->SetBranchAddress("has_sig_c", C.has_sig_c);
    t->GetEntry(0);
    for (Int_t m=0;m<N_MODELS;++m)
        C.bkg_c[m] = C.has_bkg_c[m] ? gse_readH(f, Form("bkg_c_%d",m)) : nullptr;
    for (Int_t s=0;s<N_SIG;++s) for (Int_t m=0;m<N_MODELS;++m)
        C.sig_c[s][m] = C.has_sig_c[s][m] ? gse_readH(f, Form("sig_c_%d_%d",s,m)) : nullptr;
}

// ===========================================================================
//  MASS-PAGE cache  (per model m, signal s, lxy row 0..2)
//  histograms hm_s_<m>_<s>_<row> / hm_b_<m>_<s>_<row> + the three counts.
// ===========================================================================
struct MassCache {
    // counts indexed [m][s][row]; presence flag has_[m][s][row]
    Long64_t n_lxy[N_MODELS][N_SIG][3];
    Long64_t n_sel[N_MODELS][N_SIG][3];
    Long64_t n_bdt[N_MODELS][N_SIG][3];
    Int_t    has [N_MODELS][N_SIG][3];
    TH1D*    h_sel[N_MODELS][N_SIG][3];
    TH1D*    h_bdt[N_MODELS][N_SIG][3];
};

static inline void gse_save_mass(TFile* f, MassCache& M)
{
    f->cd();
    TTree* t = new TTree("mass_scalars", "");
    t->Branch("n_lxy", M.n_lxy, Form("n_lxy[%d][%d][3]/L", N_MODELS, N_SIG));
    t->Branch("n_sel", M.n_sel, Form("n_sel[%d][%d][3]/L", N_MODELS, N_SIG));
    t->Branch("n_bdt", M.n_bdt, Form("n_bdt[%d][%d][3]/L", N_MODELS, N_SIG));
    t->Branch("has",   M.has,   Form("has[%d][%d][3]/I",   N_MODELS, N_SIG));
    t->Fill();
    t->Write("mass_scalars", TObject::kOverwrite);
    for (Int_t m=0;m<N_MODELS;++m) for (Int_t s=0;s<N_SIG;++s) for (Int_t r=0;r<3;++r) {
        gse_writeH(f, M.h_sel[m][s][r], Form("hm_s_%d_%d_%d",m,s,r));
        gse_writeH(f, M.h_bdt[m][s][r], Form("hm_b_%d_%d_%d",m,s,r));
    }
}

static inline void gse_load_mass(TFile* f, MassCache& M)
{
    TTree* t = (TTree*)f->Get("mass_scalars");
    t->SetBranchAddress("n_lxy", M.n_lxy);
    t->SetBranchAddress("n_sel", M.n_sel);
    t->SetBranchAddress("n_bdt", M.n_bdt);
    t->SetBranchAddress("has",   M.has);
    t->GetEntry(0);
    for (Int_t m=0;m<N_MODELS;++m) for (Int_t s=0;s<N_SIG;++s) for (Int_t r=0;r<3;++r) {
        M.h_sel[m][s][r] = M.has[m][s][r] ? gse_readH(f, Form("hm_s_%d_%d_%d",m,s,r)) : nullptr;
        M.h_bdt[m][s][r] = M.has[m][s][r] ? gse_readH(f, Form("hm_b_%d_%d_%d",m,s,r)) : nullptr;
    }
}

// ===========================================================================
//  TABLES / SIGNIFICANCE cache  (everything make_eff_tables_pdf computes from
//  data before it starts drawing, i.e. lines 599-952 of getSignalEff_noMET.C)
// ===========================================================================
struct TablesCache {
    // signal counts
    Long64_t N_kin_s [N_SIG][N_MODELS];
    Long64_t N_sv_s  [N_SIG][N_MODELS];
    Long64_t N_lxy   [N_SIG][N_MODELS][3];
    Long64_t N_k_lxy [N_SIG][N_MODELS][3];
    Long64_t N_p_lxy [N_SIG][N_MODELS][3];
    Long64_t N_sv_lxy[N_SIG][N_MODELS][3];
    Long64_t N_b_lxy [N_SIG][N_MODELS][3];
    // QCD counts
    Long64_t N_tot_q [N_QCD][N_MODELS];
    Long64_t N_kin_q [N_QCD][N_MODELS];
    Long64_t N_pre_q [N_QCD][N_MODELS];
    Long64_t N_sv_q  [N_QCD][N_MODELS];
    Long64_t N_bdt_q [N_QCD][N_MODELS];
    // thresholds + weights
    Double_t sig_thr        [N_MODELS][4];
    Double_t sig_thr_win    [N_MODELS][N_SETS][4];
    Double_t sig_thr_win_lxy [N_MODELS][N_SETS][3][4];
    Double_t w_sv_win       [N_MODELS][N_SETS];
    Double_t w_sv_win_lxy   [N_MODELS][N_SETS][3];
    // histograms (presence flags packed alongside)
    TH1D* sig_sv_cumul    [N_SIG][N_MODELS];
    TH1D* sig_sv_cumul_win[N_SIG][N_MODELS];
    TH1D* sig_sv_cumul_lxy[N_SIG][N_MODELS][3];
    TH1D* bkg_sv_c        [N_MODELS];
    TH1D* bkg_sv_c_win    [N_MODELS][N_SETS];
    TH1D* bkg_raw_c_win   [N_MODELS][N_SETS];
    TH1D* bkg_raw_c_win_lxy[N_MODELS][N_SETS][3];
    TH1D* bkg_xsw_c_win   [N_MODELS][N_SETS];
    TH1D* bkg_xsw_c_win_lxy[N_MODELS][N_SETS][3];
};

static inline void gse_save_tables(TFile* f, TablesCache& T)
{
    f->cd();
    TTree* t = new TTree("tables_scalars", "");
    t->Branch("N_kin_s",  T.N_kin_s,  Form("N_kin_s[%d][%d]/L",  N_SIG, N_MODELS));
    t->Branch("N_sv_s",   T.N_sv_s,   Form("N_sv_s[%d][%d]/L",   N_SIG, N_MODELS));
    t->Branch("N_lxy",    T.N_lxy,    Form("N_lxy[%d][%d][3]/L",    N_SIG, N_MODELS));
    t->Branch("N_k_lxy",  T.N_k_lxy,  Form("N_k_lxy[%d][%d][3]/L",  N_SIG, N_MODELS));
    t->Branch("N_p_lxy",  T.N_p_lxy,  Form("N_p_lxy[%d][%d][3]/L",  N_SIG, N_MODELS));
    t->Branch("N_sv_lxy", T.N_sv_lxy, Form("N_sv_lxy[%d][%d][3]/L", N_SIG, N_MODELS));
    t->Branch("N_b_lxy",  T.N_b_lxy,  Form("N_b_lxy[%d][%d][3]/L",  N_SIG, N_MODELS));
    t->Branch("N_tot_q",  T.N_tot_q,  Form("N_tot_q[%d][%d]/L", N_QCD, N_MODELS));
    t->Branch("N_kin_q",  T.N_kin_q,  Form("N_kin_q[%d][%d]/L", N_QCD, N_MODELS));
    t->Branch("N_pre_q",  T.N_pre_q,  Form("N_pre_q[%d][%d]/L", N_QCD, N_MODELS));
    t->Branch("N_sv_q",   T.N_sv_q,   Form("N_sv_q[%d][%d]/L",  N_QCD, N_MODELS));
    t->Branch("N_bdt_q",  T.N_bdt_q,  Form("N_bdt_q[%d][%d]/L", N_QCD, N_MODELS));
    t->Branch("sig_thr",        T.sig_thr,         Form("sig_thr[%d][4]/D", N_MODELS));
    t->Branch("sig_thr_win",    T.sig_thr_win,     Form("sig_thr_win[%d][%d][4]/D", N_MODELS, N_SETS));
    t->Branch("sig_thr_win_lxy",T.sig_thr_win_lxy, Form("sig_thr_win_lxy[%d][%d][3][4]/D", N_MODELS, N_SETS));
    t->Branch("w_sv_win",       T.w_sv_win,        Form("w_sv_win[%d][%d]/D", N_MODELS, N_SETS));
    t->Branch("w_sv_win_lxy",   T.w_sv_win_lxy,    Form("w_sv_win_lxy[%d][%d][3]/D", N_MODELS, N_SETS));
    t->Fill();
    t->Write("tables_scalars", TObject::kOverwrite);
    // histograms
    for (Int_t s=0;s<N_SIG;++s) for (Int_t m=0;m<N_MODELS;++m) {
        gse_writeH(f, T.sig_sv_cumul[s][m],     Form("sig_sv_cumul_%d_%d",s,m));
        gse_writeH(f, T.sig_sv_cumul_win[s][m], Form("sig_sv_cumul_win_%d_%d",s,m));
        for (Int_t c=0;c<3;++c)
            gse_writeH(f, T.sig_sv_cumul_lxy[s][m][c], Form("sig_sv_cumul_lxy_%d_%d_%d",s,m,c));
    }
    for (Int_t m=0;m<N_MODELS;++m) {
        gse_writeH(f, T.bkg_sv_c[m], Form("bkg_sv_c_%d",m));
        for (Int_t s=0;s<N_SETS;++s) {
            gse_writeH(f, T.bkg_sv_c_win[m][s],  Form("bkg_sv_c_win_%d_%d",m,s));
            gse_writeH(f, T.bkg_raw_c_win[m][s], Form("bkg_raw_c_win_%d_%d",m,s));
            gse_writeH(f, T.bkg_xsw_c_win[m][s], Form("bkg_xsw_c_win_%d_%d",m,s));
            for (Int_t c=0;c<3;++c) {
                gse_writeH(f, T.bkg_raw_c_win_lxy[m][s][c], Form("bkg_raw_c_win_lxy_%d_%d_%d",m,s,c));
                gse_writeH(f, T.bkg_xsw_c_win_lxy[m][s][c], Form("bkg_xsw_c_win_lxy_%d_%d_%d",m,s,c));
            }
        }
    }
}

static inline void gse_load_tables(TFile* f, TablesCache& T)
{
    TTree* t = (TTree*)f->Get("tables_scalars");
    t->SetBranchAddress("N_kin_s",  T.N_kin_s);
    t->SetBranchAddress("N_sv_s",   T.N_sv_s);
    t->SetBranchAddress("N_lxy",    T.N_lxy);
    t->SetBranchAddress("N_k_lxy",  T.N_k_lxy);
    t->SetBranchAddress("N_p_lxy",  T.N_p_lxy);
    t->SetBranchAddress("N_sv_lxy", T.N_sv_lxy);
    t->SetBranchAddress("N_b_lxy",  T.N_b_lxy);
    t->SetBranchAddress("N_tot_q",  T.N_tot_q);
    t->SetBranchAddress("N_kin_q",  T.N_kin_q);
    t->SetBranchAddress("N_pre_q",  T.N_pre_q);
    t->SetBranchAddress("N_sv_q",   T.N_sv_q);
    t->SetBranchAddress("N_bdt_q",  T.N_bdt_q);
    t->SetBranchAddress("sig_thr",         T.sig_thr);
    t->SetBranchAddress("sig_thr_win",     T.sig_thr_win);
    t->SetBranchAddress("sig_thr_win_lxy", T.sig_thr_win_lxy);
    t->SetBranchAddress("w_sv_win",        T.w_sv_win);
    t->SetBranchAddress("w_sv_win_lxy",    T.w_sv_win_lxy);
    t->GetEntry(0);
    for (Int_t s=0;s<N_SIG;++s) for (Int_t m=0;m<N_MODELS;++m) {
        T.sig_sv_cumul[s][m]     = gse_readH(f, Form("sig_sv_cumul_%d_%d",s,m));
        T.sig_sv_cumul_win[s][m] = gse_readH(f, Form("sig_sv_cumul_win_%d_%d",s,m));
        for (Int_t c=0;c<3;++c)
            T.sig_sv_cumul_lxy[s][m][c] = gse_readH(f, Form("sig_sv_cumul_lxy_%d_%d_%d",s,m,c));
    }
    for (Int_t m=0;m<N_MODELS;++m) {
        T.bkg_sv_c[m] = gse_readH(f, Form("bkg_sv_c_%d",m));
        for (Int_t s=0;s<N_SETS;++s) {
            T.bkg_sv_c_win[m][s]  = gse_readH(f, Form("bkg_sv_c_win_%d_%d",m,s));
            T.bkg_raw_c_win[m][s] = gse_readH(f, Form("bkg_raw_c_win_%d_%d",m,s));
            T.bkg_xsw_c_win[m][s] = gse_readH(f, Form("bkg_xsw_c_win_%d_%d",m,s));
            for (Int_t c=0;c<3;++c) {
                T.bkg_raw_c_win_lxy[m][s][c] = gse_readH(f, Form("bkg_raw_c_win_lxy_%d_%d_%d",m,s,c));
                T.bkg_xsw_c_win_lxy[m][s][c] = gse_readH(f, Form("bkg_xsw_c_win_lxy_%d_%d_%d",m,s,c));
            }
        }
    }
}

#endif // GETSIGNALEFF_CACHE_H
