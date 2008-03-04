#ifndef TT_generator_h
#define TT_generator_h

#include <iostream>
#include "Rtypes.h"
#include <math.h>
#include "TT_params.h"
#include "TT_nucleus.h"
#include "TT_event.h"
#include "TT_drawer.h"
#include "TH1F.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"

using namespace std;

class TT_generator {
 public :

  TT_params  *f_params;
  TT_nucleus *f_nucleus;

  TFile     *f_files[2];
  TTree     *f_trees[2];

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
  Int_t   f_tree_N_MA;
  Float_t *f_tree_MA;
  Float_t *f_tree_MA_weights;
  Int_t   f_tree_process;
  Float_t f_tree_enuqe;
  Float_t f_tree_Q2qe;
  
  static const Int_t f_N_Processes=TT_params::f_N_Processes;
  Double_t   *f_accepted_points[f_N_Processes];
  Double_t   *f_total_points[f_N_Processes];
  Double_t    f_total_accepted;
  Int_t       f_process[f_N_Processes];
  Double_t   *f_rate_est[f_N_Processes];
  Double_t   *f_rate_err[f_N_Processes];
  TT_drawer **f_drawer[f_N_Processes];

  Double_t   f_total_to_accept;
  Bool_t     f_keep_going;
  Double_t   f_event_display;

  Int_t   f_N_rate_bins;
  //Bool_t *f_flux_bins_on;

  TT_generator(TT_params *params,TT_nucleus *nucleus);
  virtual ~TT_generator();
  void Setup_processes();
  Bool_t Generate_events();
  Bool_t Reject_process(Int_t,Int_t);
  Bool_t Keep_going();
  void Report_progress();
  void Fill_tree(TT_event*);
  void Update_rates(Int_t,Int_t);
  void Make_graphs();
  void Write_shuffled_tree();
};

#endif

#ifdef TT_generator_cxx
TT_generator::TT_generator(TT_params *params,TT_nucleus *nucleus)
{
  f_params=params;
  f_nucleus=nucleus;
  f_total_to_accept=params->f_N_Events;

  f_N_rate_bins=f_params->f_N_rate_bins;
  //f_flux_bins_on=f_params->f_flux_bins_on;

  f_tree_N_MA      =f_params->f_N_MA;
  f_tree_MA        =new Float_t[f_tree_N_MA];
  f_tree_MA_weights=new Float_t[f_tree_N_MA];
  for (Int_t iMA=0;iMA<f_tree_N_MA;iMA++) {
    f_tree_MA[iMA]=f_params->f_MA[iMA];
  }

  //set up files and trees
  TString file_names[2];
  TString tree_names[2];
  file_names[0]=f_params->f_filename+".tmp";
  file_names[1]=f_params->f_filename;
  tree_names[0]="do_not_use_this_tree";
  tree_names[1]="tree";
  for (Int_t i=0;i<2;i++) {
    f_files[i]=new TFile(file_names[i],"RECREATE");
    f_files[i]->cd();
    f_trees[i]=new TTree(tree_names[i],"");
    
    f_trees[i]->Branch("xs_branch",&f_tree_xs,"xs/F");
    f_trees[i]->Branch("k_branch",f_tree_k,"k[4]/F");
    f_trees[i]->Branch("kprime_branch",f_tree_kprime,"kprime[4]/F");
    f_trees[i]->Branch("q_branch",f_tree_q,"q[4]/F");
    f_trees[i]->Branch("p_branch",f_tree_p,"p[4]/F");
    f_trees[i]->Branch("pprime_branch",f_tree_pprime,"pprime[4]/F");
    f_trees[i]->Branch("mag_p_branch",&f_tree_mag_p,"mag_p/F");
    f_trees[i]->Branch("cth_branch",&f_tree_cth_pq,"cth_pq/F");
    f_trees[i]->Branch("phi_branch",&f_tree_phi_p,"phi/F");
    f_trees[i]->Branch("E_branch",&f_tree_E,"E/F");
    f_trees[i]->Branch("N_MA_branch",&f_tree_N_MA,"N_MA/I");
    f_trees[i]->Branch("MA_branch",f_tree_MA,"MA[N_MA]/F");
    f_trees[i]->Branch("MA_weights_branch",f_tree_MA_weights,"MA_weights[N_MA]/F");
    f_trees[i]->Branch("process_branch",&f_tree_process,"process/I");
    f_trees[i]->Branch("enuqe_branch",&f_tree_enuqe,"enuqe/F");
    f_trees[i]->Branch("Q2qe_branch",&f_tree_Q2qe,"Q2qe/F");
  }

  f_total_accepted=0;
  for (Int_t iProcess=0;iProcess<f_N_Processes;iProcess++) {
    f_accepted_points[iProcess]=new Double_t  [f_N_rate_bins];
    f_total_points   [iProcess]=new Double_t  [f_N_rate_bins];
    f_rate_est       [iProcess]=new Double_t  [f_N_rate_bins];
    f_rate_err       [iProcess]=new Double_t  [f_N_rate_bins];
    f_drawer         [iProcess]=new TT_drawer*[f_N_rate_bins];
    f_process        [iProcess]=iProcess;

    for (Int_t iBin=0;iBin<f_N_rate_bins;iBin++) {
      f_accepted_points[iProcess][iBin]=0;
      f_total_points   [iProcess][iBin]=0;
      f_rate_est       [iProcess][iBin]=0.0;
      f_rate_err       [iProcess][iBin]=0.0;
      f_drawer         [iProcess][iBin]=0;
    }
  }
}

TT_generator::~TT_generator() {
  if (f_tree_MA)         delete [] f_tree_MA;
  if (f_tree_MA_weights) delete [] f_tree_MA_weights;

  for (Int_t i=0;i<2;i++) {
    //delete f_trees[i];
    delete f_files[i];
  }
  for (Int_t iProcess=0;iProcess<f_N_Processes;iProcess++) {
    if (f_accepted_points[iProcess]) delete [] f_accepted_points[iProcess];
    if (f_total_points   [iProcess]) delete [] f_total_points   [iProcess];
    if (f_rate_est       [iProcess]) delete [] f_rate_est       [iProcess];
    if (f_rate_err       [iProcess]) delete [] f_rate_err       [iProcess];

    for (Int_t iBin=0;iBin<f_N_rate_bins;iBin++) {
      if (f_drawer[iProcess][iBin]) delete f_drawer[iProcess][iBin];
    }

    if (f_drawer         [iProcess]) delete [] f_drawer         [iProcess];
  }
  
  //remove temporary file
  gROOT->ProcessLine(".! rm "+f_params->f_filename+".tmp");
}

#endif // #ifdef TT_generator_cxx

