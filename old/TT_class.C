#define TT_class_cxx
#include "TT_class.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TCanvas.h"

////////////////////////////////////////////////////////////////
void TT_class::Compute_XS()
{
  Integrate();

  Double_t factor=pow(f_params->f_GF*f_params->f_coscab,2)/8.0;
  factor=factor/pow(TMath::Pi(),2)/f_k_lower[0]/f_kprime_lower[0];
  f_d3sigma_d3kprime     *=factor;
  f_d3sigma_d3kprime_err *=factor;

  Double_t p_mag=pow(pow(f_kprime_lower[1],2)+pow(f_kprime_lower[2],2)+pow(f_kprime_lower[3],2),0.5);
  Double_t another_factor=2.0*TMath::Pi()*p_mag*f_kprime_lower[0];
  f_d2sigma_de_dcostheta     = f_d3sigma_d3kprime    *another_factor;
  f_d2sigma_de_dcostheta_err = f_d3sigma_d3kprime_err*another_factor;
}

////////////////////////////////////////////////////////////////
void TT_class::Integrate()
{
  Int_t N=300;

  Double_t results=0.0;
  Double_t results2=0.0;

  for (Int_t i=0;i<N;i++) {
    Double_t p_lower[4];
    Double_t E=0.0;
    Double_t result=0.0;
    Bool_t blocked=0;

    if (f_params->f_do_Smith_Moniz) {
      Double_t p_mag=f_params->f_rand.Uniform(0.0,f_params->f_SM_p_fermi);
      Double_t costh=f_params->f_rand.Uniform(-1.0,1.0);
      Double_t sinth=pow(1.0-costh*costh,0.5);
      Double_t phi=f_params->f_rand.Uniform(0.0,2*TMath::Pi());
      p_lower[1]=p_mag*sinth*cos(phi);
      p_lower[2]=p_mag*sinth*sin(phi);
      p_lower[3]=p_mag*costh;

      //according to equation (18)
      p_lower[0]=pow(pow(p_lower[1],2)+pow(p_lower[2],2)+pow(p_lower[3],2)+pow(f_params->f_M,2),0.5);
      //E not needed

      Double_t p2_fin=0.0;
      for (Int_t i=1;i<4;i++) p2_fin+=pow(p_lower[i]+f_q_lower[i],2);

      for (Int_t i=0;i<4;i++) {
	cout << "p_final_lower[" << i << "]: " << p_lower[i]+f_q_lower[i] << endl;
      }

      blocked=p2_fin<pow(f_params->f_SM_p_fermi,2);
      //cout << "p2_fin: " << p2_fin << endl;
      //cout << "blocked: " << blocked << endl << endl;
    }

    result=Evaluate_integrand(E,p_lower);
    
    if (f_params->f_do_Pauli_blocking && blocked) result=0.0;
    
    results+=result;
    results2+=result*result;
  }

  Double_t V=0.0;
  if (f_params->f_do_Smith_Moniz) V=4.0/3.0*TMath::Pi()*pow(f_params->f_SM_p_fermi,3);

  f_d3sigma_d3kprime     = V*results/N;
  f_d3sigma_d3kprime_err = V*pow(results2/N-pow(results/N,2),0.5)/pow(N,0.5);
  //cout << "ans, errs = " << f_d3sigma_d3kprime << "," << f_d3sigma_d3kprime_err << endl;
}

////////////////////////////////////////////////////////////////
Double_t TT_class::Evaluate_integrand(Double_t E,Double_t p_lower[4])
{

  Double_t p_upper[4];
  f_params->Upper(f_params->f_g,p_lower,p_upper);
  Double_t q2=0.0;
  for (Int_t mu=0;mu<4;mu++) q2+=f_q_lower[mu]*f_q_upper[mu];

  Double_t Re_H_upper[4][4];
  Double_t Im_H_upper[4][4];
  Double_t M=f_params->f_M;

  for (Int_t mu=0;mu<4;mu++) {
    for (Int_t nu=0;nu<4;nu++) {
      Re_H_upper[mu][nu] = -f_params->f_g[mu][nu]*f_params->H1(q2);
      Re_H_upper[mu][nu]+= p_upper[mu]*p_upper[nu]*f_params->H2(q2)/M/M;
      Re_H_upper[mu][nu]+= -f_q_upper[mu]*f_q_upper[nu]*f_params->H4(q2)/M/M;
      Re_H_upper[mu][nu]+= (p_upper[mu]*f_q_upper[nu]+f_q_upper[mu]*p_upper[nu])*f_params->H5(q2)/2.0/M/M;
      for (Int_t kappa=0;kappa<4;kappa++) {
        for (Int_t lambda=0;lambda<4;lambda++) {
	  //epsilon_upper=-epsilon_lower (according to Mathworld)
          Im_H_upper[mu][nu] = -f_params->f_epsilon_lower[mu][nu][kappa][lambda]*p_lower[kappa]*f_q_lower[lambda]*f_params->H3(q2)/(2.0*M*M);
        }
      }
      Re_H_upper[mu][nu]*=M*M;
      Im_H_upper[mu][nu]*=M*M;
    }
  }

  Double_t Re_contraction=0.0;
  Double_t Im_contraction=0.0;
  for (Int_t mu=0;mu<4;mu++) {
    for (Int_t nu=0;nu<4;nu++) {
      Re_contraction+=f_Re_L_lower[mu][nu]*Re_H_upper[mu][nu];
      Re_contraction-=f_Im_L_lower[mu][nu]*Im_H_upper[mu][nu];
      Im_contraction+=f_Re_L_lower[mu][nu]*Im_H_upper[mu][nu];
      Im_contraction+=f_Im_L_lower[mu][nu]*Re_H_upper[mu][nu];
    }
  }
  //cout << "Re_contraction = " << Re_contraction << endl;
  //cout << "Im_contraction = " << Im_contraction << endl;

  ////from equation 16
  //Double_t E_pprime=mass_nuc_A + q_lower[0];
  //E_pprime -= sqrt(pow(mass_nuc_Am1,2)+pow(p_lower[1],2)+pow(p_lower[2],2)+pow(p_lower[3],2));

  Double_t ans=0.0;
  if (f_params->f_do_Smith_Moniz) {
    //from equation (19)
    Double_t E_pprime=f_q_lower[0]-f_params->f_SM_e_bind+pow(M*M+pow(p_lower[1],2)+pow(p_lower[2],2)+pow(p_lower[3],2),0.5);
    ans=Re_contraction*SF_SM(E,p_lower)/p_lower[0]/E_pprime;
  }

  return ans;
}

////////////////////////////////////////////////////////////////
Double_t TT_class::SF_SM(Double_t E,Double_t p_lower[4]) {
  //(equation 21)
  Double_t p2=pow(p_lower[1],2)+pow(p_lower[2],2)+pow(p_lower[3],2);
  if ( p2>pow(f_params->f_SM_p_fermi,2) ) return 0.0;
  else return 3.0*pow(f_params->f_SM_p_fermi,3)/(4.0*TMath::Pi());
}
