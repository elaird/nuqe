#ifndef TV_class_h
#define TV_class_h

#include <iostream>
#include "Rtypes.h"
#include "TT_params.h"
#include "TH1D.h"
#include "TT_nucleus.h"

using namespace std;

class TV_class {
 public :

  TT_params *f_params;
  TT_nucleus *f_nucleus;

  Bool_t f_bad_event;

  Double_t f_k_lower[4];
  Double_t f_k_upper[4];

  Double_t f_q_lower[4];
  Double_t f_q_upper[4];

  Double_t f_kprime_lower[4];
  Double_t f_kprime_upper[4];

  Double_t f_Re_L_lower[4][4];  
  Double_t f_Im_L_lower[4][4];  

  Double_t f_cos_theta_q;
  Double_t f_cos_theta_p;
  Double_t f_phi_p;
  Double_t f_mag_p;

  Double_t f_d2sigma_dw_dq;
  Double_t f_d2sigma_dw_dq_err;
  Double_t f_d4sigma_dw_dq_bold_domega_p;

  TV_class(TT_params *params,TT_nucleus *nucleus,Double_t Enu,Double_t q_lower[4]);
  //TV_class(TT_params *params,TT_nucleus *nucleus,Double_t Enu,Double_t w,Double_t cth_q_bold,Double_t cth,Double_t phi);
  TV_class();
  virtual ~TV_class();

  void Init(TT_params *params,TT_nucleus *nucleus,Double_t Enu,Double_t w,Double_t q_bold,Double_t cth,Double_t phi);
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
TV_class::TV_class(TT_params *params,TT_nucleus *nucleus,Double_t Enu,Double_t q_lower[4])
{
  f_params = params;
  f_nucleus= nucleus;
  Init_vectors(Enu,q_lower);
  Init_L();
}

//TV_class::TV_class(TT_params *params,TT_nucleus *nucleus,Double_t Enu,Double_t w,Double_t cth_qbold,Double_t cth_p,Double_t phi_p)
//{
//  f_histo=new TH1D("histo","",100,-1.0,1.0);
//  f_params = params;
//  f_nucleus= nucleus;
//  f_bad_event=0;
//
//  Double_t p_lep_sq=pow(Enu-w,2)-pow(f_params->f_m_lep,2);
//
//  Double_t a=1.0;
//  Double_t b=-2.0*Enu*cth_qbold;
//  Double_t c=Enu*Enu-p_lep_sq;
//
//  Double_t disc=b*b-4.0*a*c;
//  Double_t root1,root2;
//  if (disc<0.0) { 
//    f_bad_event=1;
//    cout << "d<0!" << endl;
//  }
//  else {
//    root1=(-b+sqrt(disc))/(2.0*a);
//    root2=(-b-sqrt(disc))/(2.0*a);
//    Double_t ans1=a*root1*root1+b*root1+c;
//    Double_t ans2=a*root2*root2+b*root2+c;
//    
//    Double_t thing1=(Enu-root1*cth_qbold)/sqrt(p_lep_sq);
//    Double_t thing2=(Enu-root2*cth_qbold)/sqrt(p_lep_sq);
//    if (root1>0.0 && root2>0.0) {
//      printf("root1=%8.6g, ans1=%8.6g, thing1=%8.6g, Q21=%8.6g\n",root1,ans1,thing1,root1*root1-w*w);
//      printf("root2=%8.6g, ans2=%8.6g, thing2=%8.6g, Q22=%8.6g\n",root2,ans2,thing2.deq,root2*root2-w*w);
//      cout << endl;
//    }
//  }
//
//  Double_t q_bold=1.0;
//
//  //if (f_params->f_do_phony_gen) {
//  //  f_bad_event=0;
//  //  cos_theta_q=0.5;
//  //}
//
//  if (!f_bad_event) {
//    Double_t q_lower[4];
//    q_lower[0]=w;
//    q_lower[1]=0.0;
//    q_lower[2]=-q_bold*sqrt(1.0-pow(cth_qbold,2));
//    q_lower[3]=q_bold*cth_qbold;
//	  
//    Init_vectors(Enu,q_lower);
//    f_cos_theta_p=cth_p;
//    f_phi_p=phi_p;
//  }
//}

TV_class::TV_class()
{
}

TV_class::~TV_class()
{
}


#endif // #ifdef TV_class_cxx

