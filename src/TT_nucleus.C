#define TT_nucleus_cxx
#include "TT_nucleus.h"
#include "TMath.h"
#include "TGraph.h"
#include "TFile.h"

////////////////////////////////////////////////////////////////////////
void TT_nucleus::compute_p2_avg()
{
  double numerator=0.0;
  double denominator=0.0;
  int NP=1000000;

  for (int iP=0;iP<NP;iP++) {
    double p=f_AS_MF_gen_p_max*(iP+0.0)/NP;
    double nmf=n_mf(p);
    numerator+=pow(p,4)*nmf;
    denominator+=pow(p,2)*nmf;
  }
  f_avg_p2=numerator/denominator;
  std::cout << "f_avg_p2=" << f_avg_p2 << std::endl;
}

////////////////////////////////////////////////////////////////////////
Double_t TT_nucleus::n_mf(Double_t p)
{
  Double_t ans=1.0;
  ans+=f_n_mf_C*pow(p,2);
  ans+=f_n_mf_D*pow(p,4);
  ans+=f_n_mf_E*pow(p,6);
  ans+=f_n_mf_F*pow(p,8);

  ans*=f_n_mf_A*TMath::Exp(-f_n_mf_B*p*p);
  return ans*(f_N+f_Z)/(8.0*TMath::Pi());
}

////////////////////////////////////////////////////////////////////////
Double_t TT_nucleus::n_corr(Double_t p)
{
  Double_t ans=0.0;
  ans+=f_n_corr_A*TMath::Exp(-f_n_corr_B*p*p);
  ans+=f_n_corr_C*TMath::Exp(-f_n_corr_D*p*p);
  return ans*(f_N+f_Z)/(8.0*TMath::Pi());
}

////////////////////////////////////////////////////////////////////////
Double_t TT_nucleus::n_sm(Double_t p)
{
  if (p>f_SM_p_fermi) return 0.0;
  if (f_SM_p_fermi==0.0) {
    std::cout << "f_SM_p_fermi=" << f_SM_p_fermi << std::endl;
    return 0.0;
  }
  
  return 3.0*(f_N+f_Z)/(8.0*TMath::Pi()*pow(f_SM_p_fermi,3.0));
}

////////////////////////////////////////////////////////////////////////
void TT_nucleus::n_plot()
{
  Int_t n_points=1000;
  Double_t p_min=0.0;
  Double_t p_max=1.0;
  
  Double_t    p[n_points];
  Double_t   mf[n_points];
  Double_t corr[n_points];
  Double_t  tot[n_points];
  Double_t   sm[n_points];

  for (Int_t i=0;i<n_points;i++) {
    p[i]=i*(p_max-p_min)/n_points + p_min;
    sm[i]=n_sm(p[i]);
    mf[i]=n_mf(p[i]);
    corr[i]=n_corr(p[i]);
    tot[i]=mf[i]+corr[i];
  }

  TGraph sm_graph(n_points,p,sm);
  sm_graph.SetName("sm");
  sm_graph.SetTitle(";p (GeV);n_{SM} (GeV^{-3})");

  TGraph mf_graph(n_points,p,mf);
  mf_graph.SetName("mf");
  mf_graph.SetTitle(";p (GeV);n_{MF} (GeV^{-3})");

  TGraph corr_graph(n_points,p,corr);
  corr_graph.SetName("corr");
  corr_graph.SetTitle(";p (GeV);n_{corr} (GeV^{-3})");
  
  TGraph tot_graph(n_points,p,tot);
  tot_graph.SetName("tot");
  tot_graph.SetTitle(";p (GeV);n_{tot} (GeV^{-3})");
  
  TFile f("../helper/p_dists.root","RECREATE");
  sm_graph.Write();
  mf_graph.Write();
  corr_graph.Write();
  tot_graph.Write();
  f.Close();

}

////////////////////////////////////////////////////////////////////////
Double_t TT_nucleus::P_corr(Double_t p_mag,Double_t E)
{
  if (p_mag==0) {
    std::cout << "warning: p_mag=0" << std::endl;
    return 0.0;
  }

  Double_t beta=(f_Z+f_N-2.0)/(f_Z+f_N-1.0);
  Double_t alpha=3.0/(4.0*f_avg_p2*beta);

  Double_t ans=n_corr(p_mag)*f_M_n/p_mag*sqrt(alpha/TMath::Pi());

  Double_t arg=2.0*f_M_n*beta*(E-f_E2-p_mag*p_mag/(2.0*f_M_residual_nucleus));
  if (arg<=0) {
    std::cout << "warning: arg<=0" << std::endl;
    return 0.0;
  }

  Double_t p2_min=beta*p_mag-sqrt(arg);
  Double_t p2_max=beta*p_mag+sqrt(arg);
  p2_min=p2_min*p2_min;
  p2_max=p2_max*p2_max;

  ans*=(TMath::Exp(-alpha*p2_min)-TMath::Exp(-alpha*p2_max));
  return ans;
}
////////////////////////////////////////////////////////////////////////
void TT_nucleus::init_neutron() {
  f_M_target_nucleus=f_M_n;
  f_M_residual_nucleus=0.0;
  
  f_AS_SF_p_fermi=0.209;
  //actually compute these
  f_AS_MF_gen_p_max=0.4;
  f_AS_corr_gen_p_max=1.0;
  
  f_SM_p_fermi=0.005;
  f_SM_e_bind=0.000;
  f_N=1;
  f_Z=0;
  
  f_n_mf_A=pow(f_fm_conv,-3)*2.74;
  f_n_mf_B=pow(f_fm_conv,-2)*3.33;
  f_n_mf_C=pow(f_fm_conv,-2)*6.66;
  f_n_mf_D=pow(f_fm_conv,-4)*0.00;
  f_n_mf_E=pow(f_fm_conv,-6)*0.00;
  f_n_mf_F=pow(f_fm_conv,-8)*0.00;
  
  f_n_corr_A=pow(f_fm_conv,-3)*0.326;
  f_n_corr_B=pow(f_fm_conv,-2)*1.40;
  f_n_corr_C=pow(f_fm_conv,-3)*0.0263;
  f_n_corr_D=pow(f_fm_conv,-2)*0.22;

  f_E1_n=0.0;
  f_E1_p=0.0;
  f_E2=0.0;
  f_avg_p2=0.0;
}

////////////////////////////////////////////////////////////////////////
void TT_nucleus::compute_e2_oxygen16() {
  f_E2=f_amu_conv*(14.0032420+14.0030740+14.0085953)/3.0;//M_{A-2}
  //f_E2+=f_M_p+f_M_n; //2M
  f_E2+=2.0*f_M_n; //2M
  f_E2-=f_amu_conv*15.9949146;//-M_{A}

  //f_E2=0.02633;
  std::cout << "f_E2=" << f_E2 << std::endl;
}

////////////////////////////////////////////////////////////////////////
void TT_nucleus::compute_e1_oxygen16() {

  double totalOcc=0.0;
  totalOcc+= 2*0.06;//2 s
  totalOcc+=10*0.14;//1 d
  //alpha_F
  totalOcc+= 6*0.80;//1 p
  totalOcc+= 2*0.87;//1 s
  std::cout << "occ=" << totalOcc << std::endl;

  f_E1_p=0.0;
  //        N   occ.   E
  f_E1_p+=( 4 * 0.14 *(-4.65) );//1 d 3/2
  f_E1_p+=( 2 * 0.06 *  0.08  );//2 s 1/2
  f_E1_p+=( 6 * 0.14 *  0.59  );//1 d 5/2
  //alpha_F
  f_E1_p+=( 2 * 0.80 * 12.11  );//1 p 1/2
  f_E1_p+=( 4 * 0.80 * 18.44  );//1 p 3/2
  f_E1_p+=( 2 * 0.87 * 45     );//1 s 1/2

  f_E1_p/=totalOcc;
  std::cout << "f_E1_p=" << f_E1_p << std::endl;

  f_E1_n=0.0;
  //        N   occ.   E
  f_E1_n+=( 4 * 0.14 *(-0.93) );//1 d 3/2
  f_E1_n+=( 2 * 0.06 *  3.27  );//2 s 1/2
  f_E1_n+=( 6 * 0.14 *  4.15  );//1 d 5/2
  //alpha_F
  f_E1_n+=( 2 * 0.80 * 15.65  );//1 p 1/2
  f_E1_n+=( 4 * 0.80 * 21.80  );//1 p 3/2
  f_E1_n+=( 2 * 0.87 * 47     );//1 s 1/2
  f_E1_n/=totalOcc;
  std::cout << "f_E1_n=" << f_E1_n << std::endl;
}

////////////////////////////////////////////////////////////////////////
void TT_nucleus::compute_e1_carbon12() {
  double totalOcc=0.0;
  totalOcc+= 2*0.00;//2 s
  totalOcc+=10*0.09;//1 d
  //alpha_F
  totalOcc+= 6*0.58;//1 p
  totalOcc+= 2*0.79;//1 s
  std::cout << "occ=" << totalOcc << std::endl;

  f_E1_p=0.0;
  //      N   occ.   E
  f_E1_p+=( 4 * 0.09 *(-3.39) );//1 d 3/2
  f_E1_p+=( 2 * 0.00 *  1.10  );//2 s 1/2
  f_E1_p+=( 6 * 0.09 *  1.86  );//1 d 5/2
  //alpha_F
  f_E1_p+=( 2 * 0.58 *  4.95  );//1 p 1/2
  f_E1_p+=( 4 * 0.58 * 18.72  );//1 p 3/2
  f_E1_p+=( 2 * 0.79 * 35     );//1 s 1/2
  f_E1_p/=totalOcc;
  
  //use same numbers for neutrons
  f_E1_n=f_E1_p;

  std::cout << "f_E1_n=" << f_E1_n << std::endl;
}

////////////////////////////////////////////////////////////////////////
void TT_nucleus::compute_e2_carbon12() {
  f_E2=f_amu_conv*(10.0129370+10.0168531+10.0135337)/3.0;//M_{A-2}
  f_E2+=f_M_p+f_M_n; //2M
  f_E2-=f_amu_conv*12.0000000;//-M_{A}
  std::cout << "f_E2=" << f_E2 << std::endl;
}

////////////////////////////////////////////////////////////////////////
void TT_nucleus::init_oxygen16() {
  f_M_target_nucleus = 15.9949146*f_amu_conv;
  f_M_residual_nucleus = 15.0030654*f_amu_conv;
  //f_M_residual_nucleus = 15.0001089*f_amu_conv;//for anti-nu reactions
  
  f_AS_SF_p_fermi=0.209;
  //actually compute these
  //f_AS_MF_gen_p_max=0.4;
  //f_AS_corr_gen_p_max=1.0;
  f_AS_MF_gen_p_max=0.6;
  f_AS_corr_gen_p_max=1.5;
  
  f_SM_p_fermi=0.225;
  f_SM_e_bind=0.027;
  //f_SM_e_bind=0.0;
  //f_SM_p_fermi=0.001;
  //f_SM_e_bind=0.000;
  f_N=8;
  f_Z=8;
  
  f_n_mf_A=pow(f_fm_conv,-3)*2.74;
  f_n_mf_B=pow(f_fm_conv,-2)*3.33;
  f_n_mf_C=pow(f_fm_conv,-2)*6.66;
  f_n_mf_D=pow(f_fm_conv,-4)*0.00;
  f_n_mf_E=pow(f_fm_conv,-6)*0.00;
  f_n_mf_F=pow(f_fm_conv,-8)*0.00;
  
  f_n_corr_A=pow(f_fm_conv,-3)*0.326;
  f_n_corr_B=pow(f_fm_conv,-2)*1.40;
  f_n_corr_C=pow(f_fm_conv,-3)*0.0263;
  f_n_corr_D=pow(f_fm_conv,-2)*0.22;
  
  //f_E1=0.01918;//for protons
  //f_E1=0.02232;//for neutrons
  //f_avg_p2=pow(0.1622,2);
  compute_e1_oxygen16();
  compute_e2_oxygen16();
  compute_p2_avg();
}

////////////////////////////////////////////////////////////////////////
void TT_nucleus::init_carbon12() {

  f_M_target_nucleus   = 12.0000000*f_amu_conv;
  f_M_residual_nucleus = 11.0114338*f_amu_conv;
  //f_M_residual_nucleus = 11.0093055*f_amu_conv;//for anti-nu reactions
  
  f_AS_SF_p_fermi=0.209;
  //actually compute these
  //f_AS_MF_gen_p_max=0.4;
  //f_AS_corr_gen_p_max=1.0;
  f_AS_MF_gen_p_max=0.6;
  f_AS_corr_gen_p_max=1.5;
  
  f_SM_p_fermi=0.225;
  f_SM_e_bind=0.027;
  //f_SM_e_bind=0.0;
  //f_SM_p_fermi=0.001;
  //f_SM_e_bind=0.000;
  f_N=6;
  f_Z=6;
  
  f_n_mf_A=pow(f_fm_conv,-3)*2.61;
  f_n_mf_B=pow(f_fm_conv,-2)*2.66;
  f_n_mf_C=pow(f_fm_conv,-2)*3.54;
  f_n_mf_D=pow(f_fm_conv,-4)*0.00;
  f_n_mf_E=pow(f_fm_conv,-6)*0.00;
  f_n_mf_F=pow(f_fm_conv,-8)*0.00;
  
  f_n_corr_A=pow(f_fm_conv,-3)*0.426;
  f_n_corr_B=pow(f_fm_conv,-2)*1.60;
  f_n_corr_C=pow(f_fm_conv,-3)*0.0237;
  f_n_corr_D=pow(f_fm_conv,-2)*0.22;
  
  compute_e1_carbon12();
  compute_e2_carbon12();
  compute_p2_avg();
}

