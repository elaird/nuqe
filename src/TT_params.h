#ifndef TT_params_h
#define TT_params_h

#include <iostream>
#include "Rtypes.h"
#include "TString.h"
#include "TH1F.h"
#include "TRandom2.h"
#include "TRandom3.h"

using namespace std;

class TT_params {
  public :

  TRandom* f_rand;
  
  //from PDG
  static const Double_t f_M_proton=0.93827;
  static const Double_t f_M_neutron=0.93956;
  static const Double_t f_GF=1.1664e-5;
  static const Double_t f_coscab=0.974;
  static const Double_t f_unit_conv=0.389379e-27;  
  
  //source?
  static const Double_t f_Mpi=0.13957;
  static const Double_t f_target_gA=-1.267;
  //static const Double_t f_target_axial_mass=1.03;
  static const Double_t f_target_vector_mass=0.84;
  static const Double_t f_xi=3.706;

  Int_t f_N_MA;
  Double_t *f_MA;
  Int_t f_MA_base_index;

  Double_t f_g[4][4];
  Double_t f_epsilon_lower[4][4][4][4];

  static const Double_t f_m_lep=0.10566;
  //static const Double_t f_m_lep=0.000511;

  TH1F   *f_flux_histo;
  Int_t   f_N_rate_bins;
  Double_t *f_rate_regions;
  //Bool_t *f_rate_bins_on;

  Bool_t f_do_Pauli_blocking;
  Bool_t f_reject_q_tilde_lt_zero;
  Bool_t f_do_deForest_prescription;
  Bool_t f_do_zero_FP;

  Double_t f_numerical_threshold;

  static const Int_t f_N_Processes=5;
  Bool_t f_processes_on[f_N_Processes];

  Double_t f_M_target;
  Double_t f_M_recoil;

  Int_t f_N_Events;
  Int_t f_N_successes;
  Double_t f_rate_factor;

  TString f_filename;

  TT_params(char *config_file);
  virtual ~TT_params();

  virtual void Init(Int_t init_style);
  virtual void Init_val0();
  virtual void Init_val1();
  virtual void Init_val2();
  virtual void Init_val3();
  virtual void Init_val4();
  virtual void Init_val5();
  virtual void Init_val6();
  virtual void Init_val7();
  virtual void Init_valgrind();

  void Set_rand_type_and_seed(Int_t,Int_t);
  void Set_rate_bins(Int_t);
  void Set_flux_mode(Int_t);
  void Set_MA_mode(Int_t);
  void Set_masses(Double_t,Double_t);
  void Set_hacks_to_be_fixed();
  Bool_t Check_for_problems();

  virtual void Init_g();
  virtual void Init_epsilon_lower();
  Int_t Power(Int_t base,Int_t exp);
  Int_t Count_transpositions(const Int_t a[]);
  void  Upper(Double_t a_lower[4],Double_t a_upper[4]);

  //form factors are from Katori's Llewellyn-Smith BooNE memo
  Double_t G_D (Double_t q2);
  //  Double_t G_pE(Double_t q2);
  //  Double_t G_nE(Double_t q2);
  //  Double_t G_pM(Double_t q2);
  //  Double_t G_nM(Double_t q2);
  //  Double_t G_VE(Double_t q2);
  //  Double_t G_VM(Double_t q2);
  //  Double_t F1_old(Double_t q2);
  //  Double_t F2_old(Double_t q2);
  Double_t F1(Double_t q2);
  Double_t F2(Double_t q2);
  Double_t FA(Double_t q2,Double_t target_axial_mass);
  Double_t FP(Double_t q2,Double_t target_axial_mass);
  Double_t H1(Double_t q2,Double_t target_axial_mass);
  Double_t H2(Double_t q2,Double_t target_axial_mass);
  Double_t H3(Double_t q2,Double_t target_axial_mass);
  Double_t H4(Double_t q2,Double_t target_axial_mass);
  Double_t H5(Double_t q2,Double_t target_axial_mass);
};

#endif

#ifdef TT_params_cxx
TT_params::TT_params(char *config_file)
{
  f_rand=0;
  f_flux_histo=0;
  f_rate_regions=0;
  //f_rate_bins_on=0;
  Init_g();
  Init_epsilon_lower();
}

TT_params::~TT_params() {
  if (f_flux_histo)  delete f_flux_histo;
  if (f_rate_regions)  delete [] f_rate_regions;
  //if (f_rate_bins_on) delete [] f_rate_bins_on;
  if (f_MA)           delete [] f_MA;
  if (f_rand) delete f_rand;
}

#endif // #ifdef TT_params_cxx

