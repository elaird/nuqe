//implementing nucl-th/0512004v4
// ./m to compile

////list of ambiguities/issues

//code
//--check complex numbers (complex<Double_t>??)
//--make stuff const?
//--classify?

//physics
//--initialize q in dsimgadwdq3?
//--de Forest prescription
//--ksi in F2?
//--Mpi+-0?
//--poles in Smith-Moniz form factors?
//--check+/- in form factor denoms
//--G,coscab values
//--M value fixed?
//--use delta function to eliminate the integral over E?
//--overall sign of epsilon
//--epsilon upper/lower -1 proof
//--figure out why no imaginary stuff is showing up
//--what if contraction is imaginary?
//--SM delta function
//--sign of e_bind
//--understand equation (21)
//--normalize SF to #neutrons

#include <iostream>
#include <complex>
#include "TMath.h"
#include "TRandom3.h"

using namespace std;

Double_t dsigmadwdq3       (Double_t E_nu,Double_t q3,Double_t w,Double_t theta,Double_t phi);
Double_t d3sigmad3kprime   (Double_t E_nu,Double_t kprime_lower[4]);
void     g_init        (Double_t a[][4]);
void     epsilon_lower_init(Double_t e[][4][4][4]);
Int_t    power         (Int_t base,Int_t exp);
Int_t    count_transpositions(const Int_t a[]);
Double_t dot           (Double_t g[4][4],Double_t a[4],Double_t b[4]);
void     upper         (Double_t g[4][4],Double_t a_lower[4],Double_t a_upper[4]);
Double_t integral      (Double_t k_lower[4],Double_t q_lower[4]);
Double_t integrand     (const complex<Double_t> L_lower[4][4],Double_t q_lower[4],Double_t p_lower[4]);
Double_t H1(Double_t q2,Double_t M);
Double_t H2(Double_t q2,Double_t M);
Double_t H3(Double_t q2,Double_t M);
Double_t H4(Double_t q2,Double_t M);
Double_t H5(Double_t q2,Double_t M);
Double_t SF_SM(Double_t p_lower[4]);

//make const?
Double_t g[4][4];
Double_t epsilon_lower[4][4][4][4];

const Double_t mass_nuc_A = 12.0*0.931;
const Double_t mass_nuc_Am1 = 11.0*0.931;

const Double_t SM_p_fermi=0.220;
const Double_t SM_e_bind=0.025;
const Double_t M=0.93956;

Int_t main(void) {
  g_init(g);
  epsilon_lower_init(epsilon_lower);

  TRandom3 b;
  b.SetSeed(0);
  for (Int_t i=0;i<20;i++) {
    Double_t theta=b.Uniform(0.0,2*TMath::Pi());
    Double_t phi  =b.Uniform(0.0,2*TMath::Pi());
    dsigmadwdq3(1.0,0.3,0.8,theta,phi);
  }
  return 0;
}

//equation (1)
Double_t d3sigmad3kprime(Double_t E_nu,Double_t kprime[4]) {
  //Double_t q_lower[4] = {w,0.0,0.0,0.0};
  //q_lower[1]=q3*sin(theta)*cos(phi);
  //q_lower[2]=q3*sin(theta)*sin(phi);
  //q_lower[3]=q3*cos(theta);
  //Double_t k_lower[4] = {E_nu,0.0,0.0,E_nu};
  //
  //Double_t ans=GF*GF*coscab*coscab*q3/4/TMath::Pi()/E_nu/E_nu;
  //ans *= integral(k_lower,q_lower);
  //return ans;
  return 1.0;
}


////equation (8)
//Double_t dsigmadwdq3(Double_t E_nu,Double_t q3,Double_t w,Double_t theta,Double_t phi) {
//  const Double_t GF=1.0;
//  const Double_t coscab=1.0;
//
//  Double_t q_lower[4] = {w,0.0,0.0,0.0};
//  q_lower[1]=q3*sin(theta)*cos(phi);
//  q_lower[2]=q3*sin(theta)*sin(phi);
//  q_lower[3]=q3*cos(theta);
//  Double_t k_lower[4] = {E_nu,0.0,0.0,E_nu};
//
//  Double_t ans=GF*GF*coscab*coscab*q3/4/TMath::Pi()/E_nu/E_nu;
//  ans *= integral(k_lower,q_lower);
//  return ans;
//}


Double_t integral(Double_t k_lower[4],Double_t q_lower[4]) {

  Double_t kprime_lower[4];
  for (Int_t i=0;i<4;i++) kprime_lower[i] = k_lower[i] - q_lower[i];

  Double_t k_upper[4];	    upper(g,k_lower,k_upper);
  Double_t kprime_upper[4]; upper(g,kprime_lower,kprime_upper);

  Double_t kk=dot(g,k_lower,kprime_lower);

  complex<Double_t> L_lower[4][4];
  for (Int_t mu=0;mu<4;mu++) {
    for (Int_t nu=0;nu<4;nu++) {
      L_lower[mu][nu] = complex<Double_t>(2.0*(k_lower[mu]*kprime_lower[nu]+kprime_lower[mu]*k_lower[nu] - kk*g[mu][nu]),0.0);
      for (Int_t rho=0;rho<4;rho++) {
        for (Int_t sigma=0;sigma<4;sigma++) {
          L_lower[mu][nu]-= complex<Double_t>(0.0,2.0*epsilon_lower[mu][nu][rho][sigma]*k_upper[rho]*kprime_upper[sigma]);
        }
      }
      //cout << L_lower[mu][nu] << endl;
    }
  }


  Bool_t doSM=1;
  if (!doSM) return 0.0;


  TRandom3 r;
  r.SetSeed(0);
  Int_t N=400;
  Double_t p_lower[4];

  Double_t results=0.0;
  Double_t results2=0.0;
  for (Int_t i=0;i<N;i++) {
    for (Int_t j=1;j<4;j++) p_lower[j]=r.Uniform(0,SM_p_fermi);
    //according to equation (18)
    p_lower[0]=sqrt(pow(p_lower[1],2)+pow(p_lower[2],2)+pow(p_lower[3],2)+M*M);
    //and E is not needed
    Double_t result=integrand(L_lower,q_lower,p_lower);
    //cout << "integrand: " << result << endl;
    results+=result;
    results2+=result*result;
  }

  Double_t V=4.0/3.0*TMath::Pi()*pow(SM_p_fermi,3);
  Double_t ans=V*results/N;
  Double_t err=V*pow(results2/N-pow(results/N,2),0.5);

  cout << "ans, errs = " << ans << "," << err << endl;
  return ans;

}

Double_t integrand(const complex<Double_t> L_lower[4][4], Double_t q_lower[4],Double_t p_lower[4]) {

  Double_t q_upper[4];      upper(g,q_lower,q_upper);
  Double_t p_upper[4];	    upper(g,p_lower,p_upper);

  Double_t q2=0.0;
  for (Int_t mu=0;mu<4;mu++) {
    q2+=q_lower[mu]*q_upper[mu];
  }

  complex<Double_t> H_upper[4][4];
  for (Int_t mu=0;mu<4;mu++) {
    for (Int_t nu=0;nu<4;nu++) {
      H_upper[mu][nu] = complex<Double_t>(-g[mu][nu]*H1(q2,M),0.0);
      H_upper[mu][nu]+= complex<Double_t>(p_upper[mu]*p_upper[nu]*H2(q2,M)/M/M,0.0);
      H_upper[mu][nu]+= complex<Double_t>(-q_upper[mu]*q_upper[nu]*H4(q2,M)/M/M,0.0);
      H_upper[mu][nu]+= complex<Double_t>((p_upper[mu]*q_upper[nu]+q_upper[mu]*p_upper[nu])*H5(q2,M)/2.0/M/M,0.0);
      for (Int_t kappa=0;kappa<4;kappa++) {
        for (Int_t lambda=0;lambda<4;lambda++) {
	  //epsilon_upper=-epsilon_lower (according to Mathworld)
          H_upper[mu][nu]+= complex<Double_t>(0.0,-epsilon_lower[mu][nu][kappa][lambda]*p_lower[kappa]*q_lower[lambda]*H3(q2,M)/(2.0*M*M));
        }
      }
      H_upper[mu][nu]*=M*M;
      //cout << H_upper[mu][nu] << endl;
    }
  }

  complex<Double_t> contraction(0.0,0.0);
  for (Int_t mu=0;mu<4;mu++) {
    for (Int_t nu=0;nu<4;nu++) {
      contraction+=L_lower[mu][nu]*H_upper[mu][nu];
      //cout << "L: " << L_lower[mu][nu] << endl;
      //cout << "H: " << H_upper[mu][nu] << endl;
      //cout << "contraction: " << contraction << endl;
    }
  }
  ////from equation 16
  //Double_t E_pprime=mass_nuc_A + q_lower[0];
  //E_pprime -= sqrt(pow(mass_nuc_Am1,2)+pow(p_lower[1],2)+pow(p_lower[2],2)+pow(p_lower[3],2));

  //from equation (19)
  Double_t E_pprime=q_lower[0]-SM_e_bind+pow(M*M+pow(p_lower[1],2)+pow(p_lower[2],2)+pow(p_lower[3],2),0.5);
  Double_t ans=real(contraction)*SF_SM(p_lower)/p_lower[0]/E_pprime;
  return ans;
}

//(equation 21)
Double_t SF_SM(Double_t p_lower[4]) {
  
  Double_t p2=pow(p_lower[1],2)+pow(p_lower[2],2)+pow(p_lower[3],2);
  if (p2>SM_p_fermi*SM_p_fermi) return 0.0;
  else return 3.0*pow(SM_p_fermi,3)/(4.0*TMath::Pi());
}
