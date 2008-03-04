#ifndef TT_class_h
#define TT_class_h

#include <iostream>
#include <complex>
#include "Rtypes.h"
#include "TT_param_class.h"

using namespace std;

class TT_class {
 public :

  TT_param_class *f_params;

  Double_t f_k_lower[4];
  Double_t f_k_upper[4];

  Double_t f_q_lower[4];
  Double_t f_q_upper[4];

  Double_t f_kprime_lower[4];
  Double_t f_kprime_upper[4];

  Double_t f_Re_L_lower[4][4];  
  Double_t f_Im_L_lower[4][4];  
  
  Double_t f_d3sigma_d3kprime;
  Double_t f_d3sigma_d3kprime_err;

  Double_t f_d2sigma_de_dcostheta;
  Double_t f_d2sigma_de_dcostheta_err;

  TT_class(TT_param_class *params,Double_t Enu,Double_t kprime[4]);
  virtual ~TT_class();

  virtual void Init_vectors(Double_t Enu,Double_t kprime[3]);
  virtual void Init_L();

  virtual void     Compute_XS();
  virtual void     Integrate();
  virtual Double_t Evaluate_integrand(Double_t E,Double_t p_lower[4]);
  virtual Double_t SF_SM(Double_t E,Double_t p_lower[4]);
};

#endif

#ifdef TT_class_cxx
TT_class::TT_class(TT_param_class *params,Double_t Enu,Double_t kprime[4])
{
  f_params = params;
  Init_vectors(Enu,kprime);
  Init_L();
}

TT_class::~TT_class() {
}


////////////////////////////////////////////////////////////////////////
void TT_class::Init_vectors(Double_t Enu,Double_t kprime[4])
{
  f_k_lower[0]=Enu;
  f_k_lower[1]=0.0;
  f_k_lower[2]=0.0;
  f_k_lower[3]=Enu;

  for (Int_t i=0;i<4;i++) {
    f_kprime_lower[i]=kprime[i];
    f_q_lower[i]=f_k_lower[i]-f_kprime_lower[i];
  }
  
  f_params->Upper(f_params->f_g,f_k_lower,f_k_upper);
  f_params->Upper(f_params->f_g,f_kprime_lower,f_kprime_upper);
  f_params->Upper(f_params->f_g,f_q_lower,f_q_upper);
}

////////////////////////////////////////////////////////////////////////
void TT_class::Init_L()
{

  Double_t kk=0.0;
  for (Int_t i=0;i<4;i++) kk+= f_k_lower[i]*f_kprime_upper[i];

  for (Int_t mu=0;mu<4;mu++) {
    for (Int_t nu=0;nu<4;nu++) {
      f_Re_L_lower[mu][nu] = 2.0*(f_k_lower[mu]*f_kprime_lower[nu]+f_kprime_lower[mu]*f_k_lower[nu] - kk*f_params->f_g[mu][nu]);
      for (Int_t rho=0;rho<4;rho++) {
        for (Int_t sigma=0;sigma<4;sigma++) {
          f_Im_L_lower[mu][nu] = -2.0*f_params->f_epsilon_lower[mu][nu][rho][sigma]*f_k_upper[rho]*f_kprime_upper[sigma];
        }
      }
      //cout << f_Im_L_lower[mu][nu] << endl;
    }
  }
}

#endif // #ifdef TT_class_cxx

