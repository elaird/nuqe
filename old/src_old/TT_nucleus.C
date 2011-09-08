#define TT_nucleus_cxx
#include "TT_nucleus.h"
#include "TMath.h"
#include "TGraph.h"
#include "TFile.h"

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
    cout << "f_SM_p_fermi=" << f_SM_p_fermi << endl;
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
    cout << "warning: p_mag=0" << endl;
    return 0.0;
  }

  Double_t beta=(f_Z+f_N-2.0)/(f_Z+f_N-1.0);
  Double_t alpha=3.0/(4.0*f_avg_p2*beta);

  Double_t ans=n_corr(p_mag)*f_M_n/p_mag*sqrt(alpha/TMath::Pi());

  Double_t arg=2.0*f_M_n*beta*(E-f_E2-p_mag*p_mag/(2.0*f_M_residual_nucleus));
  if (arg<=0) {
    cout << "warning: arg<=0" << endl;
    return 0.0;
  }

  Double_t p2_min=beta*p_mag-sqrt(arg);
  Double_t p2_max=beta*p_mag+sqrt(arg);
  p2_min=p2_min*p2_min;
  p2_max=p2_max*p2_max;

  ans*=(TMath::Exp(-alpha*p2_min)-TMath::Exp(-alpha*p2_max));
  return ans;
}

