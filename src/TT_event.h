#ifndef TT_event_h
#define TT_event_h

#include <iostream>
#include "Rtypes.h"
#include "TT_params.h"
#include "TT_nucleus.h"
#include "TString.h"

class TT_event {
 public :

  TT_params *f_params;
  TT_nucleus *f_nucleus;

  Double_t f_k_lower[4];
  Double_t f_k_upper[4];

  Double_t f_q_lower[4];
  Double_t f_q_upper[4];

  Double_t f_kprime_lower[4];
  Double_t f_kprime_upper[4];

  Double_t f_p_lower[4];
  Double_t f_p_upper[4];

  Double_t f_pprime_lower[4];
  Double_t f_pprime_upper[4];

  Double_t f_dsigma_dall;
  Double_t *f_MA_weights;

  //redundant but useful
  Double_t f_E;
  Double_t f_mag_q;
  Double_t f_cos_theta_q;
  Double_t f_mag_p;
  Double_t f_cos_theta_pq;
  Double_t f_phi_p;
  Double_t f_enuqe;
  Double_t f_Q2qe;

  Double_t f_kappa;
  Double_t f_lambda;
  Double_t f_tau;
  Double_t f_psi;

  Int_t f_process;

  TT_event(TT_params *params,TT_nucleus *nucleus,Int_t process);
  virtual ~TT_event();

  Bool_t Init(Double_t Enu,Double_t w);
  Bool_t Init(Double_t Enu,Double_t w,Double_t q_bold,Double_t mag_p,Double_t phi_p);
  Bool_t Init(Double_t Enu,Double_t w,Double_t q_bold,Double_t mag_p,Double_t cos_theta_pq,Double_t phi_p);
  void Init_k(Double_t Enu);
  void Init_q(Double_t q[4]);

  void     Print_Stuff();
  void     Evaluate_dsigma_dall();

  Bool_t   Setup_kinematics1();
  Bool_t   Setup_kinematics2();
  Double_t Evaluate_integrand();
  void     Compute_p_and_pprime();
  Double_t f_block_momentum;

  Double_t Get_contraction(Double_t target_axial_mass);
  void     AS_SM_determine_cos_theta_pq_and_E();
  void     SM_SM_determine_cos_theta_pq_and_E();
  void     AS_MF_determine_cos_theta_pq_and_E();

  Double_t SM_argument_of_delta();
  Double_t AS_MF_argument_of_delta();

  Double_t SM_argument_of_delta_prime();
  Double_t AS_MF_argument_of_delta_prime();

  void     SM_compute_enuqe_and_Q2qe();

};

#endif

#ifdef TT_event_cxx

TT_event::TT_event(TT_params *params,TT_nucleus *nucleus,Int_t process)
{
  f_process=process;
  f_params=params;
  f_nucleus=nucleus;
  f_MA_weights=new Double_t[f_params->f_N_MA];
}

TT_event::~TT_event()
{
  delete [] f_MA_weights;
}


#endif // #ifdef TT_event_cxx

