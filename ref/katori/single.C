//--Katori BooNE memo

#include <iostream>
#include "Rtypes.h"
#include "TMath.h"
#include "TGraph.h"
#include "TFile.h"

Double_t d_sigma_d_Q2(Double_t Enu,Double_t Q2);
Double_t d_sigma_d_E_mu(Double_t E_nu,Double_t E_mu);
Double_t  A(Double_t Q2);
Double_t  B(Double_t Q2);
Double_t  C(Double_t Q2);

Double_t G_D(Double_t q2);
//Double_t G_pE(Double_t q2);
//Double_t G_nE(Double_t q2);
//Double_t G_pM(Double_t q2);
//Double_t G_nM(Double_t q2);
//Double_t G_VE(Double_t q2);
//Double_t G_VM(Double_t q2);
Double_t F1(Double_t q2);
Double_t F2(Double_t q2);
Double_t FA(Double_t q2);
Double_t FP(Double_t q2);

Bool_t do_FP_eq_0=1;
Double_t f_M_p=0.93827;
Double_t f_M_n=0.93956;
Double_t f_M=f_M_p;
Double_t f_GF=1.1664e-5;
Double_t f_coscab=0.974;
Double_t unit_conv=0.389379e-27;

//source?
Double_t f_Mpi=0.13957;
Double_t f_gA=1.267;
Double_t f_axial_mass=1.03;
//Double_t f_MV=0.71;
Double_t f_MV=0.84;
Double_t f_xi=3.706;

Bool_t f_do_Pauli_blocking;
Bool_t f_do_Smith_Moniz;
Bool_t f_do_Paschos_FF;
Bool_t f_do_Benhar_FF;
Bool_t f_do_deForest_prescription;
Bool_t f_single_nucleon_mode;

//Double_t G_VE(Double_t q2);
//Double_t G_VM(Double_t q2);
//Double_t F1(Double_t q2);
//Double_t F2(Double_t q2);
//Double_t FA(Double_t q2);
//Double_t FP(Double_t q2);

using namespace std;

///////////////////////////////////////////////////////////
Int_t main(Int_t argc,char **argv) {

  f_do_Pauli_blocking=0;
  f_do_Smith_Moniz=0;
  f_do_deForest_prescription=0;
  f_single_nucleon_mode=0;

  TFile f("single.root","RECREATE");

  const Int_t N_points=200000;
  Double_t E_nu=0.8;

  Double_t    Q2[N_points];
  Double_t xs_Q2[N_points];

  for (Int_t iPoint=0;iPoint<N_points;iPoint++) {
    Double_t Q2_min=0.0;
    Double_t Q2_max=2.0;

    Q2[iPoint]=Q2_min+(Q2_max-Q2_min)*iPoint/N_points;

    xs_Q2[iPoint]=unit_conv*d_sigma_d_Q2(E_nu,Q2[iPoint]);
  }
    
  TGraph gr_Q2(N_points,Q2,xs_Q2);
  gr_Q2.SetName("xs_Q2");
  gr_Q2.SetTitle(Form("E_{#nu}=%g GeV;Q^{2} (GeV^{2});d#sigma/dQ^{2} (cm^{2}/GeV^{2})",E_nu));
  gr_Q2.Write();


  f.Close();

  return 0;
}

///////////////////////////////////////////////////////////
Double_t d_sigma_d_E_mu(Double_t E_nu,Double_t E_mu) {
  const Double_t M_p=0.93827;
  
  Double_t Q2=2.0*M_p*(E_nu-E_mu);
  Double_t ans=d_sigma_d_Q2(E_nu,Q2)*(2.0*M_p);
  return ans;
}
///////////////////////////////////////////////////////////
Double_t d_sigma_d_Q2(Double_t Enu,Double_t Q2) {

  //from PDG
  const Double_t M_p=0.93827;
  const Double_t GF=1.1664e-5;
  const Double_t coscab=0.974;
  const Double_t mlep=0.10566;

  Double_t s_minus_u=4.0*Enu*M_p-Q2-mlep*mlep;
  Double_t ans=A(Q2)+B(Q2)*s_minus_u/(M_p*M_p)+C(Q2)*s_minus_u*s_minus_u/pow(M_p,4);

  ans*=pow(M_p*GF*coscab/Enu,2)/(8.0*TMath::Pi());
  return ans;
}

///////////////////////////////////////////////////////////
Double_t A(Double_t Q2) {
  const Double_t mlep=0.10566;
  const Double_t M_p=0.93827;

  Double_t tau=Q2/(4.0*M_p*M_p);

  Double_t ans=0.0;
  ans+=(1.0+tau)*pow(FA(-Q2),2);
  ans-=(1.0-tau)*pow(F1(-Q2),2);
  ans+=tau*(1.0-tau)*pow(F2(-Q2),2);
  ans+=4.0*tau*F1(-Q2)*F2(-Q2);

  Double_t ans2=pow(F1(-Q2)+F2(-Q2),2)+pow(FA(-Q2)+2.0*FP(-Q2),2)-4.0*(1.0+tau)*pow(FP(-Q2),2);
  ans2*=-mlep*mlep/(4.0*M_p*M_p);

  ans+=ans2;
  ans*=(mlep*mlep+Q2)/(M_p*M_p);

  return ans;
}

///////////////////////////////////////////////////////////
Double_t B(Double_t Q2) {
  const Double_t M_p=0.93827;
  return Q2*FA(-Q2)*(F1(-Q2)+F2(-Q2))/(M_p*M_p);
}

///////////////////////////////////////////////////////////
Double_t C(Double_t Q2) {
  const Double_t M_p=0.93827;
  Double_t ans=pow(FA(-Q2),2)+pow(F1(-Q2),2)+Q2*pow(F2(-Q2),2)/(4.0*M_p*M_p);
  return ans/4.0;
}

////////////////////////////////////////////////////////////////////////
Double_t G_D(Double_t q2)
{
  return 1.0/pow(1.0-q2/(f_MV*f_MV),2);
}

//////////////////////////////////////////////////////////////////////////
//Double_t G_pE(Double_t q2)
//{
//  return G_D(q2);
//}
//
//////////////////////////////////////////////////////////////////////////
//Double_t G_nE(Double_t q2)
//{
//  return 0.0;
//}
//
//////////////////////////////////////////////////////////////////////////
//Double_t G_pM(Double_t q2)
//{
//  return f_mu_p*G_D(q2);
//}
//
//////////////////////////////////////////////////////////////////////////
//Double_t G_nM(Double_t q2)
//{
//  return f_mu_n*G_D(q2);
//}
//
//////////////////////////////////////////////////////////////////////////
//Double_t G_VE(Double_t q2)
//{
//  return G_pE(q2)-G_nE(q2);
//}
//
//////////////////////////////////////////////////////////////////////////
//Double_t G_VM(Double_t q2)
//{
//  return G_pM(q2)-G_nM(q2);
//}
//
//////////////////////////////////////////////////////////////////////////
//Double_t F1_old(Double_t q2)
//{
//  return (G_VE(q2)-G_VM(q2)*q2/(4.0*f_M*f_M)) / (1.0-q2/(4.0*f_M*f_M));
//}
//
//////////////////////////////////////////////////////////////////////////
//Double_t F2_old(Double_t q2)
//{
//  return (G_VM(q2)-G_VE(q2)) / (1.0-q2/(4.0*f_M*f_M));
//}
//
////////////////////////////////////////////////////////////////////////
Double_t F1(Double_t q2)
{
  return (1.0-q2*(1.0+f_xi)/(4.0*f_M*f_M))*G_D(q2) / (1.0-q2/(4.0*f_M*f_M));
}

////////////////////////////////////////////////////////////////////////
Double_t F2(Double_t q2)
{
  return f_xi*G_D(q2) / (1.0-q2/(4.0*f_M*f_M));
}

////////////////////////////////////////////////////////////////////////
Double_t FA(Double_t q2)
{
  return f_gA/pow(1.0-q2/f_axial_mass/f_axial_mass,2);
}

////////////////////////////////////////////////////////////////////////
Double_t FP(Double_t q2)
{
  if (do_FP_eq_0) return 0.0;
  else return 2*f_M*f_M*FA(q2)/(f_Mpi*f_Mpi-q2);
}
