#ifndef TV_class_h
#define TV_class_h

#include <iostream>
#include "Rtypes.h"
#include "TT_param_class.h"
#include "TH1D.h"
#include "TT_nucleus.h"

using namespace std;

class TV_class {
 public :

  TT_param_class *f_params;
  TT_nucleus *f_nucleus;

  Bool_t f_bad_event;

  TH1D *f_histo;
  Double_t f_k_lower[4];
  Double_t f_k_upper[4];

  Double_t f_q_lower[4];
  Double_t f_q_upper[4];

  Double_t f_kprime_lower[4];
  Double_t f_kprime_upper[4];

  Double_t f_Re_L_lower[4][4];  
  Double_t f_Im_L_lower[4][4];  

  Double_t f_cos_theta_p;
  Double_t f_phi_p;

  Double_t f_d2sigma_dw_dq;
  Double_t f_d2sigma_dw_dq_err;
  Double_t f_d4sigma_dw_dq_bold_domega_p;

  TV_class(TT_param_class *params,TT_nucleus *nucleus,Double_t Enu,Double_t q_lower[4]);
  TV_class(TT_param_class *params,TT_nucleus *nucleus,Double_t Enu,Double_t w,Double_t q_bold,Double_t cth,Double_t phi);
  virtual ~TV_class();

  void Init_vectors(Double_t Enu,Double_t q[4]);
  void Init_L();

  void     Compute_XS();
  void     Evaluate_d4sigma_dw_dq_bold_domega_p();
  void     Single_nucleon_xs();
  void     SM_generate_events();
  void     SM_integrate();
  void     SM_differential();
  void     Phony_differential();
  Double_t SM_evaluate_integrand(Double_t p_lower[4],Double_t E_pprime);
  void     Get_contraction(Double_t p_lower[4],Double_t& Re_ans,Double_t& Im_ans);
  Double_t SM_determine_p(Double_t p_lower[4]);
  Double_t SM_argument_of_delta(Double_t p_lower[4],Double_t factor);
  Double_t SM_argument_of_delta_prime(Double_t p_try[4]);
  Double_t SM_argument_of_delta_prime_old(Double_t p_try[4]);
  void     SM_determine_coeffs(Double_t qmag,Double_t cos_theta_pq,Double_t& D,Double_t& a,Double_t& b,Double_t& c);
  void     Determine_stuff(Double_t p_lower[4],Double_t& q_mag,Double_t& p_mag_old,Double_t& cos_theta_pq);
  Double_t Pprime_squared(Double_t p_lower[4],Double_t factor);
  Double_t SM_compute_enuqe();

  void     AS_SF_integrate();
  Double_t AS_SF_determine_p(Double_t p_lower[4]);
};

#endif

#ifdef TV_class_cxx
TV_class::TV_class(TT_param_class *params,TT_nucleus *nucleus,Double_t Enu,Double_t q_lower[4])
{
  f_histo=new TH1D("histo","",100,-1.0,1.0);
  f_params = params;
  f_nucleus= nucleus;
  Init_vectors(Enu,q_lower);
  Init_L();
}

TV_class::TV_class(TT_param_class *params,TT_nucleus *nucleus,Double_t Enu,Double_t w,Double_t q_bold,Double_t cth,Double_t phi)
{
  f_histo=new TH1D("histo","",100,-1.0,1.0);
  f_params = params;
  f_nucleus= nucleus;
  f_bad_event=0;

  Double_t p_lep_sq=pow(Enu-w,2)-pow(f_params->f_m_lep,2);
  Double_t cos_theta_q  =(pow(Enu,2)-p_lep_sq+pow(q_bold,2))/(2.0*Enu*q_bold);
  Double_t cos_theta_lep=(pow(Enu,2)+p_lep_sq-pow(q_bold,2))/(2.0*Enu*sqrt(p_lep_sq));
  if (TMath::Abs(cos_theta_q)>1.0 || TMath::Abs(cos_theta_lep)>1.0) f_bad_event=1;
  if (f_params->f_do_phony_gen) {
    f_bad_event=0;
    cos_theta_q=0.5;
  }

  if (!f_bad_event) {
    Double_t q_lower[4];
    q_lower[0]=w;
    q_lower[1]=0.0;
    q_lower[2]=-q_bold*sqrt(1.0-pow(cos_theta_q,2));
    q_lower[3]=q_bold*cos_theta_q;
	  
    Init_vectors(Enu,q_lower);
    f_cos_theta_p=cth;
    f_phi_p=phi;
  }
}

TV_class::~TV_class()
{
  delete f_histo;
}

////////////////////////////////////////////////////////////////////////
void TV_class::Init_vectors(Double_t Enu,Double_t q[4])
{
  f_k_lower[0]=Enu;
  f_k_lower[1]=0.0;
  f_k_lower[2]=0.0;
  f_k_lower[3]=Enu;

  for (Int_t i=0;i<4;i++) {
    f_q_lower[i]=q[i];
    f_kprime_lower[i]=f_k_lower[i]-f_q_lower[i];
  }

  f_params->Upper(f_params->f_g,f_k_lower,f_k_upper);
  f_params->Upper(f_params->f_g,f_kprime_lower,f_kprime_upper);
  f_params->Upper(f_params->f_g,f_q_lower,f_q_upper);
}

////////////////////////////////////////////////////////////////////////
void TV_class::Init_L()
{

  Double_t kk=0.0;
  for (Int_t i=0;i<4;i++) kk+= f_k_lower[i]*f_kprime_upper[i];

  for (Int_t mu=0;mu<4;mu++) {
    for (Int_t nu=0;nu<4;nu++) {
      f_Re_L_lower[mu][nu] = 2.0*(f_k_lower[mu]*f_kprime_lower[nu]+f_kprime_lower[mu]*f_k_lower[nu] - kk*f_params->f_g[mu][nu]);
      f_Im_L_lower[mu][nu] = 0.0;

      if (mu==nu) continue;
      for (Int_t rho=0;rho<4;rho++) {
	if (rho==mu || rho==nu) continue;
	for (Int_t sigma=0;sigma<4;sigma++) {
	  if (sigma==mu || sigma==nu || sigma==rho) continue;
	  f_Im_L_lower[mu][nu] += -2.0*f_params->f_epsilon_lower[mu][nu][rho][sigma]*f_k_upper[rho]*f_kprime_upper[sigma];
        }
      }
      //cout << "f_Im_L_lower[" << mu << "][" << nu << "]: " << f_Im_L_lower[mu][nu] << endl;
    }
  }
}

#endif // #ifdef TV_class_cxx

