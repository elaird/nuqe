#ifndef TT_generator2_h
#define TT_generator2_h

#include <iostream>
#include "Rtypes.h"
#include <math.h>
#include "TT_params.h"
#include "TT_nucleus.h"
#include "TT_event.h"
#include "TH1F.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TT_drawer.h"

using namespace std;

class TT_generator2 {
 public :

  TT_params  *f_params;
  TT_nucleus *f_nucleus;
  TH1F       *f_flux_histo;

  TFile      *f_file;
  TTree      *f_tree;

  Float_t f_tree_xs;
  Float_t f_tree_k[4];
  Float_t f_tree_kprime[4];
  Float_t f_tree_q[4];
  Float_t f_tree_p[4];
  Float_t f_tree_pprime[4];
  Float_t f_tree_mag_p;
  Float_t f_tree_cth_pq;
  Float_t f_tree_phi_p;
  Float_t f_tree_E;
  Int_t   f_tree_process;

  static const Int_t f_N_Processes=TT_params::f_N_Processes;
  Long64_t   f_accepted_points[f_N_Processes];
  Long64_t   f_total_points[f_N_Processes];
  Long64_t   f_total_accepted;
  Int_t      f_process[f_N_Processes];
  Double_t   f_rate_est[f_N_Processes];
  Double_t   f_rate_err[f_N_Processes];
  TT_drawer *f_drawer[f_N_Processes];

  Long64_t   f_total_to_accept;
  Bool_t     f_keep_going;
  Long64_t   f_event_display;

  TT_generator2(TH1F *flux_histo,TT_params *params,TT_nucleus *nucleus);
  virtual ~TT_generator2();
  void Setup_processes();
  Bool_t Generate_events();
  Bool_t Reject_process(Int_t);
  Bool_t Keep_going();
  void Report_progress();
  void Fill_tree(TT_event*);
  void Update_rates(Int_t);
  void Finish_up();
};

#endif

#ifdef TT_generator2_cxx
TT_generator2::TT_generator2(TH1F *flux_histo,TT_params *params,TT_nucleus *nucleus)
{
  f_params=params;
  f_nucleus=nucleus;
  f_flux_histo=flux_histo;
  f_total_to_accept=params->f_N_Events;

  //set up file and tree
  f_file = new TFile(f_params->f_filename,"RECREATE");
  f_tree = new TTree("tree","");
  
  f_tree->Branch("xs_branch",&f_tree_xs,"xs/F");
  f_tree->Branch("k_branch",f_tree_k,"k[4]/F");
  f_tree->Branch("kprime_branch",f_tree_kprime,"kprime[4]/F");
  f_tree->Branch("q_branch",f_tree_q,"q[4]/F");
  f_tree->Branch("p_branch",f_tree_p,"p[4]/F");
  f_tree->Branch("pprime_branch",f_tree_pprime,"pprime[4]/F");
  f_tree->Branch("mag_p_branch",&f_tree_mag_p,"mag_p/F");
  f_tree->Branch("cth_branch",&f_tree_cth_pq,"cth_pq/F");
  f_tree->Branch("phi_branch",&f_tree_phi_p,"phi/F");
  f_tree->Branch("E_branch",&f_tree_E,"E/F");
  f_tree->Branch("process_branch",&f_tree_process,"process/I");
}

TT_generator2::~TT_generator2() {
  f_file->Close();

  delete f_file;
  //delete f_tree;

  for (Int_t iProcess=0;iProcess<f_N_Processes;iProcess++) {
    if (f_drawer[iProcess]) delete f_drawer[iProcess];
  }
}

#endif // #ifdef TT_generator2_cxx

