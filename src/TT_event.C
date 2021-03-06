#define TT_event_cxx
#include "TT_event.h"
#include "TMath.h"
#include "TRandom2.h"
#include "TCanvas.h"

////////////////////////////////////////////////////////////////////////
void TT_event::Print_Stuff()
{
  //debug info here
  std::cout << "M_i=" << sqrt(pow(f_p_lower[0],2)-pow(f_p_lower[1],2)-pow(f_p_lower[2],2)-pow(f_p_lower[3],2)) << std::endl;
  std::cout << "M_f=" << sqrt(pow(f_pprime_lower[0],2)-pow(f_pprime_lower[1],2)-pow(f_pprime_lower[2],2)-pow(f_pprime_lower[3],2)) << std::endl;
  std::cout << "M_lep=" << sqrt(pow(f_kprime_lower[0],2)-pow(f_kprime_lower[1],2)-pow(f_kprime_lower[2],2)-pow(f_kprime_lower[3],2)) << std::endl;
  SM_compute_enuqe_and_Q2qe();
  std::cout << "EnuQE=" << f_enuqe << std::endl;
}

////////////////////////////////////////////////////////////////////////
Bool_t TT_event::Init(Double_t Enu,Double_t w)
{
  Double_t qbold=w*w+2.0*f_params->f_M_target*w;
  qbold+=f_params->f_M_target*f_params->f_M_target;
  qbold-=f_params->f_M_recoil*f_params->f_M_recoil;
  qbold=sqrt(qbold);

  return Init(Enu,w,qbold,0.0,0.0);
}

////////////////////////////////////////////////////////////////////////
Bool_t TT_event::Init(Double_t Enu,Double_t w,Double_t qbold,Double_t mag_p,Double_t phi_p)
{
  f_mag_q=qbold;
  f_mag_p=mag_p;
  f_phi_p=phi_p;

  f_k_lower[0]=Enu;
  f_q_lower[0]=w;
  if (!Setup_kinematics1()) return kFALSE;

  Double_t p_lep_sq=(Enu-w)*(Enu-w)-f_params->f_m_lepton*f_params->f_m_lepton;
  Double_t cos_theta_q  =(Enu*Enu-p_lep_sq+qbold*qbold)/(2.0*Enu*qbold);
  if (cos_theta_q>1.0 || cos_theta_q<-1.0) return kFALSE;
  //Double_t cos_theta_lep=(Enu*Enu+p_lep_sq-qbold*qbold)/(2.0*Enu*sqrt(p_lep_sq));
  //if (cos_theta_lep>1.0 || cos_theta_lep<-1.0) return;
  f_cos_theta_q=cos_theta_q;

  Init_k(Enu);
  Double_t q_lower[4];
  q_lower[0]=w;
  q_lower[1]=0.0;
  q_lower[2]=-qbold*sqrt(1.0-cos_theta_q*cos_theta_q);
  q_lower[3]=qbold*cos_theta_q;
  Init_q(q_lower);

  if (!Setup_kinematics2()) return kFALSE;

  return kTRUE;
}

////////////////////////////////////////////////////////////////////////
Bool_t TT_event::Init(Double_t Enu,Double_t w,Double_t qbold,Double_t mag_p,Double_t cos_theta_pq,Double_t phi_p)
{
  //check for early exit
  Double_t pos=w+f_params->f_M_target;
  pos-=sqrt(f_params->f_M_recoil*f_params->f_M_recoil + mag_p*mag_p + qbold*qbold + 2.0*mag_p*qbold*cos_theta_pq);
  pos-=f_nucleus->f_E2+mag_p*mag_p/(2.0*f_nucleus->f_M_residual_nucleus);
  if (pos<0.0) {
    //std::cout << "pos=" << pos << std::endl;
    return kFALSE;
  }

  f_cos_theta_pq=cos_theta_pq;
  return Init(Enu,w,qbold,mag_p,phi_p);
}

////////////////////////////////////////////////////////////////////////
void TT_event::Init_k(Double_t Enu)
{
  f_k_lower[0]=Enu;
  f_k_lower[1]=0.0;
  f_k_lower[2]=0.0;
  f_k_lower[3]=Enu;
}

////////////////////////////////////////////////////////////////////////
void TT_event::Init_q(Double_t q[4])
{
  for (Int_t i=0;i<4;i++) {
    f_q_lower[i]=q[i];
    f_kprime_lower[i]=f_k_lower[i]-f_q_lower[i];
  }

  f_params->Upper(f_k_lower,f_k_upper);
  f_params->Upper(f_kprime_lower,f_kprime_upper);
  f_params->Upper(f_q_lower,f_q_upper);
}

////////////////////////////////////////////////////////////////
void TT_event::Evaluate_dsigma_dall()
{
  f_dsigma_dall=Evaluate_integrand();

  Double_t factor=pow(f_params->f_GF*f_params->f_coscab,2)/(4.0*TMath::Pi());
  factor=factor*f_mag_q/pow(f_k_lower[0],2);
  f_dsigma_dall*=factor*f_params->f_unit_conv;

  Double_t contractions[f_params->f_N_MA];
  for (Int_t iMA=0;iMA<f_params->f_N_MA;iMA++) {
    Double_t target_axial_mass=f_params->f_MA[iMA];
    contractions[iMA]=Get_contraction(target_axial_mass);
  }
  f_dsigma_dall*=contractions[f_params->f_MA_base_index];
  for (Int_t iMA=0;iMA<f_params->f_N_MA;iMA++) {
    f_MA_weights[iMA]=contractions[iMA]/contractions[f_params->f_MA_base_index];
  }
  SM_compute_enuqe_and_Q2qe();
}

////////////////////////////////////////////////////////////////
Bool_t TT_event::Setup_kinematics1()
{
  f_p_lower[0]=sqrt(f_mag_p*f_mag_p+f_params->f_M_target*f_params->f_M_target);
  f_block_momentum=-1.0;

  switch (f_process) {

  case 0:   //single nucleon
    f_cos_theta_pq=0.0;
    break;

  case 1:   //AS_SM
    f_block_momentum=f_nucleus->f_SM_p_fermi;
    AS_SM_determine_cos_theta_pq_and_E();
    if (TMath::Abs(f_cos_theta_pq)>1.0) return kFALSE;
    if (TMath::Abs(SM_argument_of_delta())>f_params->f_numerical_threshold) return kFALSE;
    break;

  case 2:   //AS_MF
    f_block_momentum=f_nucleus->f_AS_SF_p_fermi;
    AS_MF_determine_cos_theta_pq_and_E();
    if (TMath::Abs(f_cos_theta_pq)>1.0) return kFALSE;
    if (TMath::Abs(AS_MF_argument_of_delta())>f_params->f_numerical_threshold) return kFALSE;
    break;

  case 3:   //AS_corr
    f_block_momentum=f_nucleus->f_AS_SF_p_fermi;
    break;

  case 4:   //SM_SM
    f_block_momentum=f_nucleus->f_SM_p_fermi;
    f_p_lower[0]-=f_nucleus->f_SM_e_bind;
    SM_SM_determine_cos_theta_pq_and_E();
    if (TMath::Abs(f_cos_theta_pq)>1.0) return kFALSE;
    //if (TMath::Abs(SM_argument_of_delta())>f_params->f_numerical_threshold) return kFALSE;
    break;
  }
  
  return kTRUE;
}

////////////////////////////////////////////////////////////////
Bool_t TT_event::Setup_kinematics2()
{
  Compute_p_and_pprime();

  //fix this
  if (f_process==0) {
    f_E=f_q_lower[0]+f_params->f_M_target-f_pprime_lower[0];
  }

  if (f_process==3) {
    f_E=f_q_lower[0]+f_params->f_M_target-f_pprime_lower[0];
    //Double_t pos2=f_E-f_nucleus->f_E2-f_mag_p*f_mag_p/(2.0*f_nucleus->f_M_target_residual_nucleus);
    //if (pos2<0.0) return;
  }

  //Pauling blocking here
  if (f_params->f_do_Pauli_blocking) {
    Double_t mag_pprime=0.0;
    for (Int_t iComp=1;iComp<4;iComp++) {
      mag_pprime+=f_pprime_lower[iComp]*f_pprime_lower[iComp];
    }
    mag_pprime=sqrt(mag_pprime);
    if (mag_pprime<f_block_momentum) return kFALSE;
  }

  return kTRUE;
}

////////////////////////////////////////////////////////////////
void TT_event::Compute_p_and_pprime()
{
  Double_t sinth_q=sqrt(1.0-f_cos_theta_q*f_cos_theta_q);
  Double_t sinth_pq=sqrt(1.0-f_cos_theta_pq*f_cos_theta_pq);

  f_p_lower[1]=f_mag_p*( sinth_pq*cos(f_phi_p) );
  f_p_lower[2]=f_mag_p*( -f_cos_theta_pq*sinth_q + f_cos_theta_q*sinth_pq*sin(f_phi_p) );
  f_p_lower[3]=f_mag_p*(  f_cos_theta_pq*f_cos_theta_q + sinth_q*sinth_pq*sin(f_phi_p) );

  f_pprime_lower[0]=f_params->f_M_recoil*f_params->f_M_recoil;
  for (Int_t iComp=1;iComp<4;iComp++) {
    f_pprime_lower[iComp]=f_p_lower[iComp]+f_q_lower[iComp];
    f_pprime_lower[0]+=pow(f_pprime_lower[iComp],2);
  }
  f_pprime_lower[0]=sqrt(f_pprime_lower[0]);
}

////////////////////////////////////////////////////////////////
Double_t TT_event::Evaluate_integrand()
{
  Double_t SF_norm=0.0;
  Double_t delta_denom=0.0;

  //single nucleon
  if (f_process==0) {
    return 1.0/(f_p_lower[0]*f_mag_q);    
  }

  //AS_SM or SM_SM
  if (f_process==1 || f_process==4) {
    SF_norm=f_nucleus->n_sm(f_mag_p);
    delta_denom=TMath::Abs(SM_argument_of_delta_prime());
  }

  //AS_MF
  if (f_process==2) {
    SF_norm=f_nucleus->n_mf(f_mag_p);
    delta_denom=TMath::Abs(AS_MF_argument_of_delta_prime());
  }
  
  //AS_corr
  if (f_process==3) {
    SF_norm=f_nucleus->P_corr(f_mag_p,f_E);
    delta_denom=1.0;
  }

  Double_t ans=SF_norm/(f_p_lower[0]*f_pprime_lower[0]);
  if (delta_denom==0.0) std::cout << "vanishing denominator a" << std::endl;
  
  //-1-we converted to spherical coordinates for the initial momentum
  //-2-for some processes, the delta function in the spectral function
  //   killed the p_mag integral
  //---so we need two more factors
  return ans*f_mag_p*f_mag_p/delta_denom;
}

////////////////////////////////////////////////////////////////
Double_t TT_event::Get_contraction(Double_t target_axial_mass)
{
  Double_t Re_contraction=0.0;
  Double_t Im_contraction=0.0;

  Double_t kk=0.0;

  Double_t q_tilde_lower[4],q_tilde_upper[4];
  for (Int_t mu=0;mu<4;mu++) {
    kk+= f_k_lower[mu]*f_kprime_upper[mu];
    q_tilde_lower[mu]=f_q_lower[mu];
  }

  if (f_params->f_do_deForest_prescription) {
    q_tilde_lower[0]=f_pprime_lower[0]-f_p_lower[0];
    //Double_t bind=f_p_lower[0]+f_q_lower[0]-f_pprime_lower[0];
  }
  if (f_params->f_reject_q_tilde_lt_zero && q_tilde_lower[0]<0.0) return 0.0;

  f_params->Upper(q_tilde_lower,q_tilde_upper);
  f_params->Upper(f_p_lower,f_p_upper);

  Double_t q2=0.0;
  Double_t q2_tilde=0.0;
  for (Int_t mu=0;mu<4;mu++) {
    q2      +=f_q_lower[mu]*f_q_upper[mu];
    q2_tilde+=q_tilde_lower[mu]*q_tilde_upper[mu];
  }
  if (q2>0.0) std::cout <<       "q2=      " << q2 << std::endl;
  if (q2_tilde>0.0) std::cout << "q2_tilde=" << q2_tilde << std::endl;

  Double_t h1_factor=f_params->H1(q2_tilde,target_axial_mass)*pow(f_params->f_M_target,2);
  Double_t h2_factor=f_params->H2(q2_tilde,target_axial_mass);
  Double_t h3_factor=f_params->H3(q2_tilde,target_axial_mass)/2.0;
  Double_t h4_factor=f_params->H4(q2_tilde,target_axial_mass);
  Double_t h5_factor=f_params->H5(q2_tilde,target_axial_mass)/2.0;

  for (Int_t mu=0;mu<4;mu++) {
    for (Int_t nu=0;nu<4;nu++) {
      Double_t Re_H_upper_mu_nu=-f_params->f_g[mu][nu]*h1_factor;
      Re_H_upper_mu_nu+= f_p_upper[mu]*f_p_upper[nu]*h2_factor;
      Re_H_upper_mu_nu+= -q_tilde_upper[mu]*q_tilde_upper[nu]*h4_factor;
      Re_H_upper_mu_nu+= (f_p_upper[mu]*q_tilde_upper[nu]+q_tilde_upper[mu]*f_p_upper[nu])*h5_factor;
      Double_t Im_H_upper_mu_nu=0.0;

      Double_t Re_L_lower_mu_nu = 2.0*(f_k_lower[mu]*f_kprime_lower[nu]+f_kprime_lower[mu]*f_k_lower[nu] - kk*f_params->f_g[mu][nu]);
      Double_t Im_L_lower_mu_nu = 0.0;

      if (mu!=nu) {
	for (Int_t kappa=0;kappa<4;kappa++) {
	  if (kappa==mu || kappa==nu) continue;
	  //only one value of lambda contributes
	  Int_t lambda=6-mu-nu-kappa;
	  //epsilon_upper=-epsilon_lower (according to Mathworld)
	  Im_L_lower_mu_nu += -2.0*f_params->f_epsilon_lower[mu][nu][kappa][lambda]*f_k_upper[kappa]*f_kprime_upper[lambda];
	  Im_H_upper_mu_nu += -f_params->f_epsilon_lower[mu][nu][kappa][lambda]*f_p_lower[kappa]*q_tilde_lower[lambda]*h3_factor;
	}
	//std::cout << "Im_H_upper[" << mu << "][" << nu << "]: " << Im_H_upper_mu_nu << std::endl;
      }

      Re_contraction+=Re_L_lower_mu_nu*Re_H_upper_mu_nu;
      Re_contraction-=Im_L_lower_mu_nu*Im_H_upper_mu_nu;
      Im_contraction+=Re_L_lower_mu_nu*Im_H_upper_mu_nu;
      Im_contraction+=Im_L_lower_mu_nu*Re_H_upper_mu_nu;
    }
  }

  //std::cout << "Re:" << f_Re_contraction << std::endl;
  //std::cout << "Im:" << f_Im_contraction << std::endl;
  return Re_contraction;
}

////////////////////////////////////////////////////////////////
void TT_event::AS_SM_determine_cos_theta_pq_and_E()
{
  f_E=f_params->f_M_target+f_nucleus->f_SM_e_bind-f_p_lower[0];
  Double_t temp=f_q_lower[0]+f_params->f_M_target-f_E;
  f_cos_theta_pq=temp*temp-f_params->f_M_recoil*f_params->f_M_recoil - f_mag_p*f_mag_p - f_mag_q*f_mag_q;
  f_cos_theta_pq /=(2.0*f_mag_p*f_mag_q);
}

////////////////////////////////////////////////////////////////
void TT_event::SM_SM_determine_cos_theta_pq_and_E()
{
  f_E=f_params->f_M_target-f_p_lower[0];
  Double_t temp=f_q_lower[0]+f_params->f_M_target-f_E;
  f_cos_theta_pq=temp*temp-f_params->f_M_recoil*f_params->f_M_recoil - f_mag_p*f_mag_p - f_mag_q*f_mag_q;
  f_cos_theta_pq /=(2.0*f_mag_p*f_mag_q);
}

////////////////////////////////////////////////////////////////
Double_t TT_event::SM_argument_of_delta()
{
  Double_t var=f_params->f_M_recoil*f_params->f_M_recoil;
  var+=f_mag_p*f_mag_p;
  var+=f_mag_q*f_mag_q;
  var+=2.0*f_mag_p*f_mag_q*f_cos_theta_pq;
  return f_q_lower[0]+f_params->f_M_target-f_E-sqrt(var);
}

////////////////////////////////////////////////////////////////
Double_t TT_event::SM_argument_of_delta_prime()
{
  Double_t var=f_params->f_M_recoil*f_params->f_M_recoil;
  var+=f_mag_p*f_mag_p;
  var+=f_mag_q*f_mag_q;
  var+=2.0*f_mag_p*f_mag_q*f_cos_theta_pq;
  return -f_mag_p*f_mag_q/sqrt(var);
}

////////////////////////////////////////////////////////////////
void TT_event::AS_MF_determine_cos_theta_pq_and_E()
{
  f_E=f_nucleus->f_E1+f_mag_p*f_mag_p/(2.0*f_nucleus->f_M_residual_nucleus);
  Double_t temp=f_q_lower[0]+f_params->f_M_target-f_E;
  f_cos_theta_pq=temp*temp-f_params->f_M_recoil*f_params->f_M_recoil - f_mag_p*f_mag_p - f_mag_q*f_mag_q;
  f_cos_theta_pq /=(2.0*f_mag_p*f_mag_q);
}

////////////////////////////////////////////////////////////////
Double_t TT_event::AS_MF_argument_of_delta()
{
  Double_t var=f_params->f_M_recoil*f_params->f_M_recoil;
  var+=f_mag_p*f_mag_p;
  var+=f_mag_q*f_mag_q;
  var+=2.0*f_mag_p*f_mag_q*f_cos_theta_pq;
  return f_q_lower[0]+f_params->f_M_target-f_E-sqrt(var);
}

////////////////////////////////////////////////////////////////
Double_t TT_event::AS_MF_argument_of_delta_prime()
{
  return SM_argument_of_delta_prime();
}

////////////////////////////////////////////////////////////////
void TT_event::SM_compute_enuqe_and_Q2qe()
{
  Double_t M=f_params->f_M_target-f_nucleus->f_SM_e_bind;
  Double_t num,denom;

  num=pow(f_params->f_M_recoil,2)-pow(f_params->f_m_lepton,2)-M*M+2.0*M*f_kprime_lower[0];
  denom=M-f_kprime_lower[0]+f_kprime_lower[3];

  f_enuqe=0.5*num/denom;

  f_Q2qe=2.0*f_enuqe*(f_kprime_lower[0]-f_kprime_lower[3]);
  f_Q2qe+=f_params->f_m_lepton*f_params->f_m_lepton;

  f_kappa=sqrt(f_q_lower[1]*f_q_lower[1]+f_q_lower[2]*f_q_lower[2]+f_q_lower[3]*f_q_lower[3])/(2.0*f_params->f_M_target);
  f_lambda=(f_q_lower[0]-f_nucleus->f_SM_e_bind)/(2.0*f_params->f_M_target);
  f_tau=f_kappa*f_kappa-f_lambda*f_lambda;
  f_psi=(f_lambda-f_tau) / sqrt( (1.0+f_lambda)*f_tau + f_kappa*sqrt(f_tau*(1.0+f_tau)) );
  Double_t eta=f_nucleus->f_SM_p_fermi/f_params->f_M_target;
  Double_t ksi=sqrt(1.0+eta*eta)-1.0;
  f_psi/=sqrt(ksi);
}
