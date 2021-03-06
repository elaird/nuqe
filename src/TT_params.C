#define TT_params_cxx
#include "TT_params.h"
#include "TFile.h"
#include <math.h>

////////////////////////////////////////////////////////////////////////
void TT_params::Init(Int_t init_style)
{
  Set_hacks_to_be_fixed();

  if (init_style<0) {
    Set_rate_bins(1);
    Set_flux_mode("numu",1200);
    //Set_Enu_max(2.0);

    //Set_MA_mode(1);
    //Set_MA_mode(10);
    Set_MA_mode(120);
    //Set_MA_mode(100);

    Set_masses(f_m_muon,f_M_proton,f_M_neutron);

    f_do_Pauli_blocking       =true;
    f_reject_q_tilde_lt_zero  =true;
    f_do_deForest_prescription=false;
    f_do_zero_FP              =true;
    f_do_antinu               =false;

    f_filename="~ted/work/boone/ccqe/events/may09_events_2M_proc4_ma120_nfp_mu_1200.root";
    //f_filename="~ted/work/boone/ccqe/events/test.root";

    f_N_Events=2000000;

    f_processes_on[0]=0;
    f_processes_on[1]=0;
    f_processes_on[2]=0;
    f_processes_on[3]=0;
    f_processes_on[4]=1;
  }
  else {
    f_N_Events=200000;
  }

  if (init_style==0) Init_val0();
  if (init_style==1) Init_val1();
  if (init_style==2) Init_val2();
  if (init_style==3) Init_val3();
  if (init_style==4) Init_val4();
  if (init_style==5) Init_val5();
  if (init_style==6) Init_val6();
  if (init_style==7) Init_val7();
  if (init_style==8) Init_val8();
  if (init_style==9) Init_val9();
}


////////////////////////////////////////////////////////////////////////
void TT_params::Init_valgrind()
{
  Set_rate_bins(1);
  Set_flux_mode("numu",800);
  Set_MA_mode(1);
  Set_masses(f_m_muon,f_M_neutron,f_M_proton);
  
  f_do_Pauli_blocking       =kTRUE;
  f_reject_q_tilde_lt_zero  =kTRUE;
  f_do_deForest_prescription=kFALSE;
  f_do_zero_FP              =kFALSE;
  f_do_antinu               =kFALSE;
  
  f_filename="../events/events_memtest.root";
  
  f_N_Events=4;
  
  f_processes_on[0]=1;
  f_processes_on[1]=0;
  f_processes_on[2]=0;
  f_processes_on[3]=1;
  f_processes_on[4]=0;
}

////////////////////////////////////////////////////////////////////////
void TT_params::Init_val0()
{
  Set_rate_bins(1);
  Set_flux_mode("numu",800);
  Set_MA_mode(1);
  Set_masses(f_m_muon,f_M_neutron,f_M_proton);

  f_do_Pauli_blocking       =kTRUE;
  f_reject_q_tilde_lt_zero  =kTRUE;
  f_do_deForest_prescription=kFALSE;
  f_do_zero_FP              =kFALSE;
  f_do_antinu               =kFALSE;

  f_filename="../events/events_val0.root";

  f_processes_on[0]=1;
  f_processes_on[1]=0;
  f_processes_on[2]=0;
  f_processes_on[3]=0;
  f_processes_on[4]=0;
}

////////////////////////////////////////////////////////////////////////
void TT_params::Init_val1()
{
  Set_rate_bins(1);
  Set_flux_mode("numu",800);
  Set_MA_mode(1);
  Set_masses(f_m_muon,f_M_neutron,f_M_neutron);

  f_do_Pauli_blocking       =kTRUE;
  f_reject_q_tilde_lt_zero  =kTRUE;
  f_do_deForest_prescription=kTRUE;
  f_do_zero_FP              =kFALSE;
  f_do_antinu               =kFALSE;

  f_filename="../events/events_val1.root";

  f_processes_on[0]=0;
  f_processes_on[1]=1;
  f_processes_on[2]=0;
  f_processes_on[3]=0;
  f_processes_on[4]=0;
}

////////////////////////////////////////////////////////////////////////
void TT_params::Init_val2()
{
  Set_rate_bins(1);
  Set_flux_mode("numu",800);
  Set_MA_mode(1);
  Set_masses(f_m_muon,f_M_neutron,f_M_neutron);

  f_do_Pauli_blocking       =kFALSE;
  f_reject_q_tilde_lt_zero  =kTRUE;
  f_do_deForest_prescription=kTRUE;
  f_do_zero_FP              =kFALSE;
  f_do_antinu               =kFALSE;

  f_filename="../events/events_val2.root";

  f_processes_on[0]=0;
  f_processes_on[1]=0;
  f_processes_on[2]=1;
  f_processes_on[3]=0;
  f_processes_on[4]=0;
}

////////////////////////////////////////////////////////////////////////
void TT_params::Init_val3()
{
  Set_rate_bins(1);
  Set_flux_mode("numu",800);
  Set_MA_mode(1);
  Set_masses(f_m_muon,f_M_neutron,f_M_neutron);

  f_do_Pauli_blocking       =kFALSE;
  f_reject_q_tilde_lt_zero  =kTRUE;
  f_do_deForest_prescription=kTRUE;
  f_do_zero_FP              =kFALSE;
  f_do_antinu               =kFALSE;

  f_filename="../events/events_val3.root";

  f_processes_on[0]=0;
  f_processes_on[1]=0;
  f_processes_on[2]=0;
  f_processes_on[3]=1;
  f_processes_on[4]=0;
}

////////////////////////////////////////////////////////////////////////
void TT_params::Init_val4()
{
  Set_rate_bins(1);
  Set_flux_mode("numu",800);
  Set_MA_mode(1);
  Set_masses(f_m_muon,f_M_neutron,f_M_neutron);

  f_do_Pauli_blocking       =kTRUE;
  f_reject_q_tilde_lt_zero  =kTRUE;
  f_do_deForest_prescription=kTRUE;
  f_do_zero_FP              =kFALSE;
  f_do_antinu               =kFALSE;

  f_filename="../events/events_val4.root";

  f_processes_on[0]=0;
  f_processes_on[1]=0;
  f_processes_on[2]=1;
  f_processes_on[3]=1;
  f_processes_on[4]=0;
}

////////////////////////////////////////////////////////////////////////
void TT_params::Init_val5()
{
  Set_rate_bins(1);
  Set_flux_mode("numu",800);
  Set_MA_mode(1);
  Set_masses(f_m_muon,f_M_neutron,f_M_proton);

  f_do_Pauli_blocking       =kTRUE;
  f_reject_q_tilde_lt_zero  =kTRUE;
  f_do_deForest_prescription=kFALSE;
  f_do_zero_FP              =kTRUE;
  f_do_antinu               =kFALSE;

  f_filename="../events/events_val5.root";

  f_processes_on[0]=0;
  f_processes_on[1]=0;
  f_processes_on[2]=0;
  f_processes_on[3]=0;
  f_processes_on[4]=1;
}

////////////////////////////////////////////////////////////////////////
void TT_params::Init_val6()
{
  Set_rate_bins(1);
  Set_flux_mode("numu",300);
  Set_MA_mode(1);
  Set_masses(f_m_muon,f_M_neutron,f_M_proton);

  f_do_Pauli_blocking       =kTRUE;
  f_reject_q_tilde_lt_zero  =kTRUE;
  f_do_deForest_prescription=kFALSE;
  f_do_zero_FP              =kTRUE;
  f_do_antinu               =kFALSE;

  f_filename="../events/events_val6.root";

  f_processes_on[0]=0;
  f_processes_on[1]=0;
  f_processes_on[2]=0;
  f_processes_on[3]=0;
  f_processes_on[4]=1;
}

////////////////////////////////////////////////////////////////////////
void TT_params::Init_val7()
{
  Set_rate_bins(20);
  Set_flux_mode("numu",0);
  //hack
  f_rate_regions[f_N_rate_bins]=5.0;
  Set_MA_mode(1);
  Set_masses(f_m_muon,f_M_neutron,f_M_proton);

  f_do_Pauli_blocking       =kTRUE;
  f_reject_q_tilde_lt_zero  =kTRUE;
  f_do_deForest_prescription=kFALSE;
  f_do_zero_FP              =kTRUE;
  f_do_antinu               =kFALSE;

  f_filename="../events/events_val7.root";

  f_processes_on[0]=0;
  f_processes_on[1]=0;
  f_processes_on[2]=0;
  f_processes_on[3]=0;
  f_processes_on[4]=1;
}

////////////////////////////////////////////////////////////////////////
void TT_params::Init_val8()
{
  Set_rate_bins(1);
  Set_flux_mode("numu",800);
  Set_MA_mode(1);
  Set_masses(f_m_muon,f_M_proton,f_M_neutron);

  f_do_Pauli_blocking       =kTRUE;
  f_reject_q_tilde_lt_zero  =kTRUE;
  f_do_deForest_prescription=kFALSE;
  f_do_zero_FP              =kTRUE;
  f_do_antinu               =kTRUE;

  f_filename="../events/events_val8.root";

  f_processes_on[0]=0;
  f_processes_on[1]=0;
  f_processes_on[2]=0;
  f_processes_on[3]=0;
  f_processes_on[4]=1;
}

////////////////////////////////////////////////////////////////////////
void TT_params::Init_val9()
{
  Set_rate_bins(20);
  Set_flux_mode("numu",0);
  //hack
  f_rate_regions[f_N_rate_bins]=5.0;
  Set_MA_mode(1);
  Set_masses(f_m_muon,f_M_proton,f_M_neutron);

  f_do_Pauli_blocking       =kTRUE;
  f_reject_q_tilde_lt_zero  =kTRUE;
  f_do_deForest_prescription=kFALSE;
  f_do_zero_FP              =kTRUE;
  f_do_antinu               =kTRUE;

  f_filename="../events/events_val9.root";

  f_processes_on[0]=0;
  f_processes_on[1]=0;
  f_processes_on[2]=0;
  f_processes_on[3]=0;
  f_processes_on[4]=1;
}

////////////////////////////////////////////////////////////////////////
void TT_params::Set_rand_type_and_seed(Int_t rand_type,Int_t seed)
{
  if (rand_type==2) {
    f_rand=new TRandom2;
  }
  if (rand_type==3) {
    f_rand=new TRandom3;
  }
  f_rand->SetSeed(seed);
}

////////////////////////////////////////////////////////////////////////
void TT_params::Set_rate_bins(Int_t bins)
{
  f_N_rate_bins=bins;
}

////////////////////////////////////////////////////////////////////////
void TT_params::Set_Enu_max(Double_t Enu_max)
{
  f_rate_regions[f_N_rate_bins]=Enu_max;
  for (Int_t j_rate_bin=0;j_rate_bin<=f_N_rate_bins;j_rate_bin++) {
    printf("y[%d]=%g\n",j_rate_bin,f_rate_regions[j_rate_bin]);
  }
}

//////////////////////////////////////////////////////////////////////////
//void TT_params::Set_flux_mode(Int_t flux_mode)
//{
//  TString file_name;
//  TString hist_name;
//  if (flux_mode==0) {
//    file_name="helper/april07_baseline_rgen610.6_flux_8gev.root";
//    hist_name="h7003";
//  }
//  else {
//    file_name=Form("helper/flux_%d.root",flux_mode);
//    hist_name="flux";
//  }
//  TFile flux_file(file_name);
//  f_flux_histo=(TH1F*)flux_file.Get(hist_name)->Clone();
//  f_flux_histo->SetDirectory(0);
//  flux_file.Close();
//
//  Double_t *x=new Double_t[f_N_rate_bins];
//  f_rate_regions=new Double_t[f_N_rate_bins+1];
//
//  for (Int_t i_rate_bin=0;i_rate_bin<f_N_rate_bins;i_rate_bin++) {
//    x[i_rate_bin]=(i_rate_bin+0.0)/f_N_rate_bins;
//  }
//  f_flux_histo->GetQuantiles(f_N_rate_bins,f_rate_regions,x);
//  f_rate_regions[f_N_rate_bins]=f_flux_histo->GetXaxis()->GetXmax();
//  delete [] x;
//
//  for (Int_t j_rate_bin=0;j_rate_bin<=f_N_rate_bins;j_rate_bin++) {
//    printf("y[%d]=%g\n",j_rate_bin,f_rate_regions[j_rate_bin]);
//  }
//
//  //f_N_rate_bins=rate_bins;
//  //f_rate_bins_on=new Bool_t[f_N_rate_bins];
//  //for (Int_t iBin=0;iBin<f_N_rate_bins;iBin++) {
//  //  //turn on by default
//  //  f_flux_bins_on[iBin]=1;
//  //
//  //  //turn off some
//  //  if (iBin==0) f_flux_bins_on[iBin]=0;
//  //  if (iBin==f_N_flux_bins-1) f_flux_bins_on[iBin]=0;
//  //
//  //  if (flux_mode==0 && iBin<3) f_flux_bins_on[iBin]=0;
//  //  if (flux_mode==0 && iBin>9) f_flux_bins_on[iBin]=0;
//  //
//  //  //this is not perfect if one interpolates
//  //  if (f_flux_histo->GetBinContent(iBin)==0.0) f_flux_bins_on[iBin]=0;
//  //}
//}

////////////////////////////////////////////////////////////////////////
void TT_params::Set_flux_mode(TString neutrino_string,Int_t energy)
{
  TString file_name;
  TString hist_name;
  if (energy==0) {
    file_name="helper/april07_baseline_rgen610.6_flux_8gev.root";
    if (neutrino_string=="numu") {
      hist_name="h7003";
    }
    else if (neutrino_string=="nue") {
      hist_name="";
    }
  }
  else {
    file_name=Form("helper/flux_%d.root",energy);
    hist_name="flux";
  }
  TFile flux_file(file_name);
  f_flux_histo=(TH1F*)flux_file.Get(hist_name)->Clone();
  f_flux_histo->SetDirectory(0);
  flux_file.Close();

  Double_t *x=new Double_t[f_N_rate_bins];
  f_rate_regions=new Double_t[f_N_rate_bins+1];

  for (Int_t i_rate_bin=0;i_rate_bin<f_N_rate_bins;i_rate_bin++) {
    x[i_rate_bin]=(i_rate_bin+0.0)/f_N_rate_bins;
  }
  f_flux_histo->GetQuantiles(f_N_rate_bins,f_rate_regions,x);
  f_rate_regions[f_N_rate_bins]=f_flux_histo->GetXaxis()->GetXmax();
  delete [] x;

  for (Int_t j_rate_bin=0;j_rate_bin<=f_N_rate_bins;j_rate_bin++) {
    printf("y[%d]=%g\n",j_rate_bin,f_rate_regions[j_rate_bin]);
  }

  //f_N_rate_bins=rate_bins;
  //f_rate_bins_on=new Bool_t[f_N_rate_bins];
  //for (Int_t iBin=0;iBin<f_N_rate_bins;iBin++) {
  //  //turn on by default
  //  f_flux_bins_on[iBin]=1;
  //
  //  //turn off some
  //  if (iBin==0) f_flux_bins_on[iBin]=0;
  //  if (iBin==f_N_flux_bins-1) f_flux_bins_on[iBin]=0;
  //
  //  if (flux_mode==0 && iBin<3) f_flux_bins_on[iBin]=0;
  //  if (flux_mode==0 && iBin>9) f_flux_bins_on[iBin]=0;
  //
  //  //this is not perfect if one interpolates
  //  if (f_flux_histo->GetBinContent(iBin)==0.0) f_flux_bins_on[iBin]=0;
  //}
}

////////////////////////////////////////////////////////////////////////
void TT_params::Set_MA_mode(Int_t mode)
{
  f_MA_base_index=0;

  if (mode==1) {
    f_N_MA=1;
    f_MA=new Double_t[f_N_MA];
    f_MA[0]=1.03;
  }

  if (mode==10) {
    f_N_MA=1;
    f_MA=new Double_t[f_N_MA];
    f_MA[0]=1.23;
  }

  if (mode==100) {
    f_N_MA=1;
    f_MA=new Double_t[f_N_MA];
    f_MA[0]=1.00;
  }

  if (mode==120) {
    f_N_MA=1;
    f_MA=new Double_t[f_N_MA];
    f_MA[0]=1.20;
  }

  if (mode==4) {
    f_N_MA=4;
    f_MA=new Double_t[f_N_MA];
    for (Int_t iMA=0;iMA<f_N_MA;iMA++) {
      Double_t base=1.03;
      f_MA[iMA]=base+0.1*(iMA-1);
    }
  }

}

////////////////////////////////////////////////////////////////////////
void TT_params::Set_masses(Double_t m_lep,Double_t m_target, Double_t m_recoil)
{
  f_m_lepton=m_lep;
  f_M_target=m_target;
  f_M_recoil=m_recoil;
}

////////////////////////////////////////////////////////////////////////
void TT_params::Set_hacks_to_be_fixed()
{
  std::cerr << "improve these hacks someday" << std::endl;
  f_numerical_threshold=2.0e-6;
  //f_numerical_threshold=2.0e-13;
  f_N_successes=200;
  f_rate_factor=2.0;
}

////////////////////////////////////////////////////////////////////////
Bool_t TT_params::Check_for_problems()
{
  if (f_do_deForest_prescription && f_processes_on[4]) {
    std::cout << "yikes-- double off-shell correction!" << std::endl;
    return kFALSE;
  }
  if (f_MA_base_index>=f_N_MA) {
    std::cout << "fix your base index!" << std::endl;
    return kFALSE;
  }
  return kTRUE;
}

////////////////////////////////////////////////////////////////////////
void TT_params::Init_g()
{
  //if you use another, convention, don't forget to change TT_params::Upper
  //(coded for speed, not f_g-changability)
  for (Int_t mu=0;mu<4;mu++) {
    for (Int_t nu=0;nu<4;nu++) {
      f_g[mu][nu]=0.0;
    }
  }
  f_g[0][0]= 1.0;
  f_g[1][1]=-1.0;
  f_g[2][2]=-1.0;
  f_g[3][3]=-1.0;
}

////////////////////////////////////////////////////////////////////////
void TT_params::Init_epsilon_lower()
{
  Int_t i[4];
  for (i[0]=0;i[0]<4;i[0]++) {
    for (i[1]=0;i[1]<4;i[1]++) {
      for (i[2]=0;i[2]<4;i[2]++) {
	for (i[3]=0;i[3]<4;i[3]++) {
	  Int_t c=0;
	  for (Int_t j=0;j<4;j++) c+=1<<i[j];
	  
	  if (c==15) {
	    f_epsilon_lower[i[0]] [i[1]] [i[2]] [i[3]]=Power(-1,Count_transpositions(i));
	    //printf("%d,%d,%d,%d: %d,%g",i[0],i[1],i[2],i[3],Count_transpositions(i),f_epsilon_lower[i[0]] [i[1]] [i[2]] [i[3]]);
	  }
	  else f_epsilon_lower[i[0]] [i[1]] [i[2]] [i[3]]=0.0;
	}
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////
Int_t TT_params::Power(Int_t base,Int_t exp) {
  Int_t ans=1;
  for (Int_t i=0;i<exp;i++) ans*=base;
  return ans;
}

////////////////////////////////////////////////////////////////////////
Int_t TT_params::Count_transpositions(const Int_t a[]) {
  //bubble sort!
  Int_t b[4];
  for (Int_t j=0;j<4;j++) b[j]=a[j];

  Int_t ans=0;
  for (Int_t inst=0;inst<3;inst++) {
    for (Int_t index=0;index<3;index++) {
      if (b[index+1]<b[index]) {
	Int_t temp=b[index+1];
	b[index+1]=b[index];
	b[index]=temp;
	ans++;
      }
    }
  }
  return ans;
}

////////////////////////////////////////////////////////////////////////
void TT_params::Upper(Double_t a_lower[4],Double_t a_upper[4])
{
  a_upper[0]=a_lower[0];
  a_upper[1]=-a_lower[1];
  a_upper[2]=-a_lower[2];
  a_upper[3]=-a_lower[3];
}

////////////////////////////////////////////////////////////////////////
Double_t TT_params::G_D(Double_t q2)
{
  return 1.0/pow(1.0-q2/(f_target_vector_mass*f_target_vector_mass),2);
}

//////////////////////////////////////////////////////////////////////////
//Double_t TT_params::G_pE(Double_t q2)
//{
//  return G_D(q2);
//}
//
//////////////////////////////////////////////////////////////////////////
//Double_t TT_params::G_nE(Double_t q2)
//{
//  return 0.0;
//}
//
//////////////////////////////////////////////////////////////////////////
//Double_t TT_params::G_pM(Double_t q2)
//{
//  return f_mu_p_factor*f_mu_N*G_D(q2);
//}
//
//////////////////////////////////////////////////////////////////////////
//Double_t TT_params::G_nM(Double_t q2)
//{
//  return f_mu_n_factor*f_mu_N*G_D(q2);
//}
//
//////////////////////////////////////////////////////////////////////////
//Double_t TT_params::G_VE(Double_t q2)
//{
//  if (f_do_Paschos_FF || f_do_Benhar_FF) return 1.0/pow(1.0-q2/(f_vector_mass*f_vector_mass),2);
//  return G_pE(q2)-G_nE(q2);
//}
//
//////////////////////////////////////////////////////////////////////////
//Double_t TT_params::G_VM(Double_t q2)
//{
//  if (f_do_Paschos_FF || f_do_Benhar_FF) return 4.706*G_VE(q2);
//  return G_pM(q2)-G_nM(q2);
//}
//
//////////////////////////////////////////////////////////////////////////
//Double_t TT_params::F1_old(Double_t q2)
//{
//  return (G_VE(q2)-G_VM(q2)*q2/(4.0*f_M*f_M)) / (1.0-q2/(4.0*f_M*f_M));
//}
//
//////////////////////////////////////////////////////////////////////////
//Double_t TT_params::F2_old(Double_t q2)
//{
//  if (f_do_Benhar_FF)  return (G_VM(q2)-G_VE(q2)) / (1.0-q2/(4.0*f_M*f_M));
//  if (f_do_Paschos_FF) return (G_VM(q2)-G_VE(q2)) / (1.0-q2/(4.0*f_M*f_M));
//  //or perhaps this?
//  //if (f_do_Paschos_FF) return (G_VM(q2)-G_VE(q2)) / (1.0-q2/(4.0*f_M*f_M)) / f_xi;
//  return (G_VM(q2)-G_VE(q2)) / (1.0-q2/(4.0*f_M*f_M));
//}

////////////////////////////////////////////////////////////////////////
Double_t TT_params::F1(Double_t q2)
{
  return (1.0-q2*(1.0+f_xi)/(4.0*f_M_target*f_M_target))*G_D(q2) / (1.0-q2/(4.0*f_M_target*f_M_target));
}

////////////////////////////////////////////////////////////////////////
Double_t TT_params::F2(Double_t q2)
{
  return f_xi*G_D(q2) / (1.0-q2/(4.0*f_M_target*f_M_target));
}

////////////////////////////////////////////////////////////////////////
Double_t TT_params::FA(Double_t q2,Double_t target_axial_mass)
{
  return f_target_gA/pow(1.0-q2/(target_axial_mass*target_axial_mass),2);
}

////////////////////////////////////////////////////////////////////////
Double_t TT_params::FP(Double_t q2,Double_t target_axial_mass)
{
  if (f_do_zero_FP) return 0.0;
  else return 2*f_M_target*f_M_target*FA(q2,target_axial_mass)/(f_Mpi*f_Mpi-q2);
}

////////////////////////////////////////////////////////////////////////
Double_t TT_params::H1(Double_t q2,Double_t target_axial_mass)
{
  Double_t tau=-q2/(4.0*f_M_target*f_M_target);
  return pow(FA(q2,target_axial_mass),2)*(1+tau) + tau*pow(F1(q2)+F2(q2),2);
}

////////////////////////////////////////////////////////////////////////
Double_t TT_params::H2(Double_t q2,Double_t target_axial_mass)
{
  Double_t tau=-q2/(4.0*f_M_target*f_M_target);
  return pow(FA(q2,target_axial_mass),2) + pow(F1(q2),2) + tau*pow(F2(q2),2);
}

////////////////////////////////////////////////////////////////////////
Double_t TT_params::H3(Double_t q2,Double_t target_axial_mass)
{
  Double_t value=2.0*FA(q2,target_axial_mass)*(F1(q2)+F2(q2));
  if (f_do_antinu) return -value;
  return value;
}

////////////////////////////////////////////////////////////////////////
Double_t TT_params::H4(Double_t q2,Double_t target_axial_mass)
{
  Double_t tau=-q2/(4.0*f_M_target*f_M_target);
  return pow(F2(q2),2)*(1-tau)/4.0 + F1(q2)*F2(q2)/2.0 + FA(q2,target_axial_mass)*FP(q2,target_axial_mass) - tau*pow(FP(q2,target_axial_mass),2);
}

////////////////////////////////////////////////////////////////////////
Double_t TT_params::H5(Double_t q2,Double_t target_axial_mass)
{
  return H2(q2,target_axial_mass);
}

