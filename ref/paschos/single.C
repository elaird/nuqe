//--implementing hep-ph/0107261v1

#include <iostream>
#include "Rtypes.h"
#include "TMath.h"
#include "TGraph.h"
#include "TFile.h"

Double_t d_sigma_d_q2(Double_t E_nu,Double_t q2);
Double_t d_sigma_d_E_mu(Double_t E_nu,Double_t E_mu);
Double_t F1(Double_t q2);
Double_t xi_F2(Double_t q2);
Double_t FA(Double_t q2);
Double_t FP(Double_t q2);
Double_t G_VE(Double_t q2);
Double_t G_VM(Double_t q2);

using namespace std;

///////////////////////////////////////////////////////////
Int_t main(Int_t argc,char **argv) {

  TFile f("single.root","RECREATE");

  const Double_t M=0.93956;
  const Double_t m=0.10566;

  const Int_t NE_nu=200;
  const Int_t N_points=600;

  Double_t E_nu_min=0.8;
  Double_t E_nu_max=5.8;
  Double_t E_nu[NE_nu];
  Double_t XS[NE_nu];

  for (Int_t iE_nu=0;iE_nu<NE_nu;iE_nu++) {

    E_nu[iE_nu]=E_nu_min + (E_nu_max-E_nu_min)*iE_nu/NE_nu;

    Double_t    Q2[N_points];
    Double_t  E_mu[N_points];
    Double_t xs_Q2  [N_points];
    Double_t xs_E_mu[N_points];

    XS[iE_nu]=0.0;
    for (Int_t iPoint=0;iPoint<N_points;iPoint++) {
      Double_t Q2_min=0.0;
      Double_t Q2_max=4.0*M*(E_nu[iE_nu]-m);

      Double_t E_mu_min=0.10566;//E_nu[iE_nu]*7.0/8.0;
      Double_t E_mu_max=E_nu[iE_nu];

      Q2  [iPoint]=Q2_min   + (Q2_max-Q2_min)*iPoint/N_points;
      E_mu[iPoint]=E_mu_min + (E_mu_max-E_mu_min)*iPoint/N_points;

      xs_Q2  [iPoint]=0.389379e-27*d_sigma_d_q2(E_nu[iE_nu],-Q2[iPoint]);
      xs_E_mu[iPoint]=0.389379e-27*d_sigma_d_E_mu(E_nu[iE_nu],E_mu[iPoint]);

      XS[iE_nu]+=xs_Q2[iPoint];
    }
    XS[iE_nu]*=(Q2[N_points-1]-Q2[0])/N_points;

    TGraph gr_Q2(N_points,Q2,xs_Q2);
    gr_Q2.SetName(Form("xs_Q2_%d",iE_nu));
    gr_Q2.SetTitle(Form("E_{#nu}=%g GeV;Q^{2} (GeV^{2});d#sigma/dQ^{2} (cm^{2}/GeV^{2})",E_nu[iE_nu]));
    gr_Q2.Write();

    TGraph gr_E_mu(N_points,E_mu,xs_E_mu);
    gr_E_mu.SetName(Form("xs_E_mu_%d",iE_nu));
    gr_E_mu.SetTitle(Form("E_{#nu}=%g GeV;E_{#mu} (GeV);d#sigma/dE_{#mu} (cm^{2}/GeV)",E_nu[iE_nu]));
    gr_E_mu.Write();

  }

  TGraph gr_XS(NE_nu,E_nu,XS);
  gr_XS.SetName("XS");
  gr_XS.SetTitle(";E_{#nu} (GeV);#sigma (cm^{2})");
  gr_XS.Write();

  f.Close();

  return 0;
}

///////////////////////////////////////////////////////////
Double_t d_sigma_d_E_mu(Double_t E_nu,Double_t E_mu) {
  const Double_t M_p=0.93827;
  
  Double_t Q2=2.0*M_p*(E_nu-E_mu);
  Double_t ans=d_sigma_d_q2(E_nu,-Q2)*(2.0*M_p);
  return ans;
}

///////////////////////////////////////////////////////////
Double_t d_sigma_d_q2(Double_t E_nu,Double_t q2) {

  //from PDG

  const Double_t M_p=0.93827;
  const Double_t M=0.93956;
  const Double_t GF=1.1664e-5;
  const Double_t coscab=0.974;
  const Double_t m=0.10566;

  Double_t s_minus_u=4.0*E_nu*M+q2-m*m;
  Double_t ans=0.0;

  ans += pow(F1(q2),2)*(q2*q2-4.0*M*M*(m*m-q2)-pow(m,4))/(4.0*M*M);
  ans += pow(xi_F2(q2),2)*(4.0*M*M*(q2*q2-pow(m,4))-q2*q2*(m*m-q2))/(16.0*pow(M,4));
  ans += pow(FA(q2),2)*(q2*q2-4.0*M*M*(m*m-q2)-pow(m,4))/(4.0*M*M);
  ans -= pow(FP(q2),2)*m*m*q2*(-q2+m*m)/(4.0*pow(M,4));
  ans += F1(q2)*xi_F2(q2)*(2.0*q2*q2+q2*m*m+pow(m,4))/(2.0*M*M);
  ans -= FA(q2)*FP(q2)*m*m*(-q2+m*m)/(2.0*M*M);
  ans += FA(q2)*(F1(q2)+xi_F2(q2))*q2*s_minus_u/(M*M);
  ans += ( pow(F1(q2),2) - pow(xi_F2(q2),2)*q2/(4.0*M*M) + pow(FA(q2),2) )*pow(s_minus_u,2)/(4.0*M*M);

  ans*=GF*GF*coscab*coscab/(8.0*TMath::Pi()*E_nu*E_nu);
  return ans;
}

///////////////////////////////////////////////////////////
Double_t G_VE(Double_t q2) {
  const Double_t MV=0.84;

  return 1.0/pow(1.0-q2/MV/MV,2);
}

///////////////////////////////////////////////////////////
Double_t G_VM(Double_t q2) {
  const Double_t MV=0.84;
  const Double_t xi=3.706;

  return (1.0+xi)/pow(1.0-q2/MV/MV,2);
}

///////////////////////////////////////////////////////////
Double_t F1(Double_t q2) {
  const Double_t M=0.93956;

  Double_t tau=-q2/(4.0*M*M);
  return (G_VE(q2)+tau*G_VM(q2))/(1.0+tau);
}

///////////////////////////////////////////////////////////
Double_t xi_F2(Double_t q2) {
  const Double_t M=0.93956;

  Double_t tau=-q2/(4.0*M*M);
  return (G_VM(q2)-G_VE(q2))/(1.0+tau);
}

///////////////////////////////////////////////////////////
Double_t FA(Double_t q2) {
  const Double_t MA=1.0;

  return -1.23/pow(1.0-q2/MA/MA,2);
}

///////////////////////////////////////////////////////////
Double_t FP(Double_t q2) {
  const Double_t M=0.93956;
  const Double_t M_pi=0.14;
  
  return 2.0*M*M*FA(q2)/(M_pi*M_pi-q2);
}

