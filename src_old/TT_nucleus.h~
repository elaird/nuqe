#ifndef TT_nucleus_h
#define TT_nucleus_h

#include <iostream>
#include "Rtypes.h"
#include <math.h>

using namespace std;

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
  Double_t f_E2;
  Double_t f_avg_p2;

  TT_nucleus(Int_t Z, Int_t N);
  virtual ~TT_nucleus();

  Double_t n_sm(Double_t p);
  Double_t n_mf(Double_t p);
  Double_t n_corr(Double_t p);
  Double_t P_corr(Double_t p_mag,Double_t E);
  void n_plot();
};

#endif

#ifdef TT_nucleus_cxx
TT_nucleus::TT_nucleus(Int_t Z=0, Int_t N=0)
{

  if (Z==0 && N==1) {
    f_M_target_nucleus=f_M_n;
    f_M_residual_nucleus=0.0;

    f_AS_SF_p_fermi=0.209;
    //actually compute these
    f_AS_MF_gen_p_max=0.4;
    f_AS_corr_gen_p_max=1.0;

    f_SM_p_fermi=0.005;
    f_SM_e_bind=0.000;
    f_N=1;
    f_Z=0;

    f_n_mf_A=pow(f_fm_conv,-3)*2.74;
    f_n_mf_B=pow(f_fm_conv,-2)*3.33;
    f_n_mf_C=pow(f_fm_conv,-2)*6.66;
    f_n_mf_D=pow(f_fm_conv,-4)*0.00;
    f_n_mf_E=pow(f_fm_conv,-6)*0.00;
    f_n_mf_F=pow(f_fm_conv,-8)*0.00;

    f_n_corr_A=pow(f_fm_conv,-3)*0.326;
    f_n_corr_B=pow(f_fm_conv,-2)*1.40;
    f_n_corr_C=pow(f_fm_conv,-3)*0.0263;
    f_n_corr_D=pow(f_fm_conv,-2)*0.22;

    f_E1=0.01918;
    f_E2=0.02633;//actually compute this
    f_avg_p2=pow(0.1622,2);//actually compute this
  }

  if (Z==8 && N==8) {
    f_M_target_nucleus = 15.9949146*f_amu_conv;
    f_M_residual_nucleus = 15.0030654*f_amu_conv;

    f_AS_SF_p_fermi=0.209;
    //actually compute these
    //f_AS_MF_gen_p_max=0.4;
    //f_AS_corr_gen_p_max=1.0;
    f_AS_MF_gen_p_max=0.6;
    f_AS_corr_gen_p_max=1.5;

    f_SM_p_fermi=0.225;
    f_SM_e_bind=0.027;
    //f_SM_e_bind=0.0;
    //f_SM_p_fermi=0.001;
    //f_SM_e_bind=0.000;
    f_N=8;
    f_Z=8;

    f_n_mf_A=pow(f_fm_conv,-3)*2.74;
    f_n_mf_B=pow(f_fm_conv,-2)*3.33;
    f_n_mf_C=pow(f_fm_conv,-2)*6.66;
    f_n_mf_D=pow(f_fm_conv,-4)*0.00;
    f_n_mf_E=pow(f_fm_conv,-6)*0.00;
    f_n_mf_F=pow(f_fm_conv,-8)*0.00;

    f_n_corr_A=pow(f_fm_conv,-3)*0.326;
    f_n_corr_B=pow(f_fm_conv,-2)*1.40;
    f_n_corr_C=pow(f_fm_conv,-3)*0.0263;
    f_n_corr_D=pow(f_fm_conv,-2)*0.22;

    //f_E1=0.01918;
    f_E1=0.02232;
    f_E2=0.02633;//actually compute this
    f_avg_p2=pow(0.1622,2);//actually compute this
  }

}

TT_nucleus::~TT_nucleus() {
}

#endif // #ifdef TT_nucleus_cxx

