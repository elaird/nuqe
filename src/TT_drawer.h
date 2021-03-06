#ifndef TT_drawer_h
#define TT_drawer_h

#include <iostream>
#include "Rtypes.h"
#include <math.h>
#include "TT_params.h"
#include "TT_nucleus.h"
#include "TT_event.h"
#include "TRandom2.h"
#include "TH1F.h"
#include "TMath.h"

class TT_drawer {
 public :

  TT_params  *f_params;
  TT_nucleus *f_nucleus;
  TT_event   *f_event;

  Double_t f_Enu_min;
  Double_t f_Enu_max;
  Double_t f_w_min;
  Double_t f_w_max;
  Double_t f_qbold_min;
  Double_t f_qbold_max;
  Double_t f_mag_p_min;
  Double_t f_mag_p_max;
  Double_t f_cos_theta_pq_min;
  Double_t f_cos_theta_pq_max;
  Double_t f_phi_p_min;
  Double_t f_phi_p_max;

  Double_t f_height;
  Double_t f_integral;

  Int_t f_process;
  Int_t f_bin;

  Bool_t f_initializing;
  Bool_t f_DIE_NOW;

  TH1D *f_init_histo;

  TT_drawer(TT_params *params,TT_nucleus *nucleus,Int_t process,Int_t bin);
  virtual ~TT_drawer();
  Double_t Flux_histo_height();
  void Init_randomly();
  void Compute_integral();
  Bool_t Event_goodness();
  Bool_t Draw_point();
  Bool_t false1();
  Bool_t false2();
  Bool_t false3();
  Bool_t false4();
  Bool_t false5();

};

#endif

#ifdef TT_drawer_cxx
TT_drawer::TT_drawer(TT_params *params,TT_nucleus *nucleus,Int_t process,Int_t bin)
{
  f_params=params;
  f_nucleus=nucleus;
  f_event = new TT_event(params,nucleus,process);

  f_process=process;
  f_bin=bin;

  f_initializing=kFALSE;
  f_DIE_NOW=kFALSE;

  f_init_histo=0;

  //f_w_min=0.0;
  //f_w_max=Enu-params->f_m_lep;
  //
  //Double_t p_lep=sqrt(pow(f_Enu-f_w_min,2)-pow(f_params->f_m_lep,2));
  //f_qbold_min=f_Enu-p_lep;
  //f_qbold_max=f_Enu+p_lep;

  //f_Enu_min=f_flux_histo->GetXaxis()->GetXmin();
  //f_Enu_max=f_flux_histo->GetXaxis()->GetXmax();

  //f_Enu_min=f_params->f_flux_histo->GetBinLowEdge(f_bin);
  //f_Enu_max=f_params->f_flux_histo->GetBinWidth(f_bin)+f_Enu_min;

  //TAxis *axis=f_params->f_flux_histo->GetXaxis();
  //Double_t axis_length=axis->GetXmax()-axis->GetXmin();;
  //f_Enu_min=axis->GetXmin() + axis_length*(f_bin+0.0)/f_params->f_N_rate_bins;
  //f_Enu_max=axis->GetXmin() + axis_length*(f_bin+1.0)/f_params->f_N_rate_bins;

  f_Enu_min=f_params->f_rate_regions[f_bin];
  f_Enu_max=f_params->f_rate_regions[f_bin+1];

  f_w_min=0.0;
  f_w_max=f_Enu_max-params->f_m_lepton;

  Double_t p_lep=sqrt(pow(f_Enu_max-f_w_min,2)-pow(f_params->f_m_lepton,2));
  f_qbold_min=f_Enu_max-p_lep;
  f_qbold_max=f_Enu_max+p_lep;

  f_mag_p_min=0.0;
  if      (f_process==1) f_mag_p_max=f_nucleus->f_SM_p_fermi;
  else if (f_process==2) f_mag_p_max=f_nucleus->f_AS_MF_gen_p_max;
  else if (f_process==3) f_mag_p_max=f_nucleus->f_AS_corr_gen_p_max;
  else if (f_process==4) f_mag_p_max=f_nucleus->f_SM_p_fermi;

  f_cos_theta_pq_min=-1.0;
  f_cos_theta_pq_max=1.0;

  f_phi_p_min=-TMath::Pi();
  f_phi_p_max=TMath::Pi();

}

TT_drawer::~TT_drawer() {
  if (f_init_histo) delete f_init_histo;
  delete f_event;
}

#endif // #ifdef TT_drawer_cxx

