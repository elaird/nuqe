#ifndef TT_param_class_h
#define TT_param_class_h

#include <iostream>
#include "Rtypes.h"
#include "TRandom3.h"

using namespace std;

class TT_param_class {
 public :

  TRandom3 f_rand;

  //from PDG
  static const Double_t f_M_p=0.93827;
  static const Double_t f_M=0.93956;
  static const Double_t f_GF=1.1664e-5;
  static const Double_t f_coscab=0.974;

  //source?
  static const Double_t f_Mpi=0.13957;
  static const Double_t f_gA=-1.267;
  static const Double_t f_axial_mass=1.03;
  static const Double_t f_MV=0.84;
  static const Double_t f_xi=3.706;

  Double_t f_g[4][4];
  Double_t f_epsilon_lower[4][4][4][4];

  static const Double_t f_m_lep=0.10566;
  //static const Double_t f_m_lep=0.000511;

  Bool_t f_do_AS_SF;
  Bool_t f_do_Pauli_blocking;
  Bool_t f_do_Smith_Moniz;
  Bool_t f_do_Smith_Moniz_gen;
  Bool_t f_do_phony_gen;
  Bool_t f_do_deForest_prescription;
  Bool_t f_single_nucleon_mode;
  Bool_t f_use_M_instead_of_MplusW;

  TT_param_class(char *config_file);
  virtual ~TT_param_class();

  virtual void Init_g();
  virtual void Init_epsilon_lower();
  Int_t Power(Int_t base,Int_t exp);
  Int_t Count_transpositions(const Int_t a[]);
  void  Upper(Double_t g[4][4],Double_t a_lower[4],Double_t a_upper[4]);

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
  Double_t FA(Double_t q2);
  Double_t FP(Double_t q2);
  Double_t H1(Double_t q2);
  Double_t H2(Double_t q2);
  Double_t H3(Double_t q2);
  Double_t H4(Double_t q2);
  Double_t H5(Double_t q2);
};

#endif

#ifdef TT_param_class_cxx
TT_param_class::TT_param_class(char *config_file)
{
  f_rand.SetSeed(0);
  f_do_Pauli_blocking=1;

  f_do_AS_SF=0;
  f_do_Smith_Moniz=0;
  f_do_Smith_Moniz_gen=1;
  f_do_phony_gen=0;
  f_do_deForest_prescription=1;
  f_single_nucleon_mode=0;
  f_use_M_instead_of_MplusW=0;
  Init_g();
  Init_epsilon_lower();
}

TT_param_class::~TT_param_class() {
}

#endif // #ifdef TT_param_class_cxx

