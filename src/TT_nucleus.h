#ifndef TT_nucleus_h
#define TT_nucleus_h

#include <iostream>
#include "Rtypes.h"
#include <math.h>

class TT_nucleus {
 public :

  static const Double_t f_amu_conv=0.9314940;//Gev/amu
  static const Double_t f_fm_conv=0.197327;//GeV.fm
  static const Double_t f_M_p=0.93827;
  static const Double_t f_M_n=0.93956;

  Double_t f_M_target_nucleus;
  Double_t f_M_residual_nucleus;
  Double_t f_AS_SF_p_fermi;
  Double_t f_AS_MF_gen_p_max;
  Double_t f_AS_corr_gen_p_max;
  Double_t f_SM_p_fermi;
  Double_t f_SM_e_bind;
  Int_t    f_N;
  Int_t    f_Z;

  Double_t f_n_mf_A;
  Double_t f_n_mf_B;
  Double_t f_n_mf_C;
  Double_t f_n_mf_D;
  Double_t f_n_mf_E;
  Double_t f_n_mf_F;

  Double_t f_n_corr_A;
  Double_t f_n_corr_B;
  Double_t f_n_corr_C;
  Double_t f_n_corr_D;

  Double_t f_E1;
  Double_t f_E1_p;
  Double_t f_E1_n;

  Double_t f_E2;
  Double_t f_avg_p2;

  TT_nucleus(Int_t Z, Int_t N);
  virtual ~TT_nucleus();

  Double_t n_sm(Double_t p);
  Double_t n_mf(Double_t p);
  Double_t n_corr(Double_t p);
  Double_t P_corr(Double_t p_mag,Double_t E);
  void n_plot();

 private :
  void compute_p2_avg();
  void init_neutron();
  void init_oxygen16();
  void init_carbon12();

  void compute_e1_oxygen16();
  void compute_e1_carbon12();

  void compute_e2_oxygen16();
  void compute_e2_carbon12();
};

#endif

#ifdef TT_nucleus_cxx
TT_nucleus::TT_nucleus(Int_t Z=0, Int_t N=0)
{
  if (Z==0 && N==1) init_neutron();
  if (Z==8 && N==8) init_oxygen16();
  if (Z==6 && N==6) init_carbon12();
}

TT_nucleus::~TT_nucleus() {
}

#endif // #ifdef TT_nucleus_cxx

