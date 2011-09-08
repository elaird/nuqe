#define TT_drawer_cxx
#include "TT_drawer.h"
#include "TT_event.h"

////////////////////////////////////////////////////////////////////////
Double_t TT_drawer::Flux_histo_height()
{
  TH1F *histo=f_params->f_flux_histo;
  Double_t Enu=f_event->f_k_lower[0];

  if (0) return histo->GetBinContent(f_params->f_flux_histo->FindBin(Enu));

  Int_t fluxbins=histo->GetNbinsX();
  Int_t bin=histo->FindBin(Enu);
  if (bin==0 || bin==fluxbins+1) return 0.0;
  
  Int_t left_bin=bin;
  Int_t right_bin=bin;
  Double_t Enu_center=histo->GetBinCenter(bin);
  if (Enu>Enu_center) {
    left_bin=bin;
    right_bin=bin+1;
  }
  else {
    left_bin=bin-1;
    right_bin=bin;
  }
  if (left_bin<1) left_bin=1;
  if (right_bin>fluxbins) right_bin=fluxbins;
  Double_t Enu_left=histo->GetBinCenter(left_bin);
  Double_t Enu_right=histo->GetBinCenter(right_bin);
  Double_t fraction=0.0;
  if (left_bin!=right_bin) fraction=(Enu-Enu_left)/(Enu_right-Enu_left);
  
  Double_t height_left=histo->GetBinContent(left_bin);
  Double_t height_right=histo->GetBinContent(right_bin);
  return height_left+fraction*(height_right-height_left);
}

////////////////////////////////////////////////////////////////////////
void TT_drawer::Init_randomly()
{
  f_initializing=kTRUE;
  f_height=0.0;

  Double_t heights[f_params->f_N_successes];

  Int_t iSuccess=0;
  while (iSuccess<f_params->f_N_successes) {
    if (!Draw_point()) continue;
    Double_t result=f_event->f_dsigma_dall*Flux_histo_height();

    if (result>=f_height/2.0) {
      heights[iSuccess]=result;
      iSuccess++;
      if (iSuccess%100==0) std::cout << iSuccess << std::endl;
    }
    if (result>f_height) {
      f_height=result;
    }
  }
  
  TString name=Form("init_histo_proc%d_bin%d",f_process,f_bin);
  f_init_histo=new TH1D(name,name,100,0.0,2.0*f_height);
  f_init_histo->SetDirectory(0);
  for (Int_t iSuccess=0;iSuccess<f_params->f_N_successes;iSuccess++) {
    f_init_histo->Fill(heights[iSuccess]);
  }
  
  f_height *=f_params->f_rate_factor;
  f_initializing=kFALSE;
}

////////////////////////////////////////////////////////////////////////
void TT_drawer::Compute_integral()
{
  f_integral=f_height;
  f_integral*=(f_Enu_max  -f_Enu_min);
  f_integral*=(f_w_max    -f_w_min);
  if (f_process>0) {
    f_integral*=(f_qbold_max-f_qbold_min);
    f_integral*=(f_mag_p_max-f_mag_p_min);
    f_integral*=(f_phi_p_max-f_phi_p_min);
  }
  if (f_process==3) f_integral*=(f_cos_theta_pq_max-f_cos_theta_pq_min);
}

////////////////////////////////////////////////////////////////////////
Bool_t TT_drawer::Event_goodness()
{
  Double_t Enu=f_params->f_rand->Uniform(f_Enu_min,f_Enu_max);
  Double_t w=f_params->f_rand->Uniform(f_w_min,f_w_max);
  if (w>Enu-f_params->f_m_lepton) return false1();

  if (f_process==0) return f_event->Init(Enu,w);
  else {
    //Double_t p_lep=sqrt((Enu-f_w_min)*(Enu-f_w_min)-(f_params->f_m_lepton)*(f_params->f_m_lepton));
    Double_t p_lep=sqrt((Enu-w)*(Enu-w)-(f_params->f_m_lepton)*(f_params->f_m_lepton));
    Double_t qbold_min=Enu-p_lep;
    Double_t qbold_max=Enu+p_lep;
    Double_t qbold=f_params->f_rand->Uniform(f_qbold_min,f_qbold_max);
    if (qbold<qbold_min || qbold>qbold_max) return false2();

    Double_t mag_p=f_params->f_rand->Uniform(f_mag_p_min,f_mag_p_max);
    Double_t phi_p=f_params->f_rand->Uniform(f_phi_p_min,f_phi_p_max);

    if (f_process==3) {
      Double_t cos_theta_pq=f_params->f_rand->Uniform(-1.0,1.0);
      return f_event->Init(Enu,w,qbold,mag_p,cos_theta_pq,phi_p);
    }
    else return f_event->Init(Enu,w,qbold,mag_p,phi_p);
  }

  return kTRUE;
}

////////////////////////////////////////////////////////////////////////
Bool_t TT_drawer::Draw_point()
{
  if (!Event_goodness()) return false4();
  f_event->Evaluate_dsigma_dall();

  Double_t rate=f_event->f_dsigma_dall*Flux_histo_height();
  if (f_height<rate) {
    TString nums=Form("proc=%d; bin=%d; ratio=%8.6g; h=%8.6g; rate=%8.6g;",f_process,f_bin,rate/f_height,f_height,rate);
    if (f_initializing) {
      std::cout << "initializing: "+nums << std::endl;
    }
    else {
      f_DIE_NOW=kTRUE;
      std::cerr << "FATAL ERROR: "+nums << std::endl;
    }
  }
  
  //rejection method
  Double_t y=f_params->f_rand->Uniform(0.0,f_height);
  if (y<rate) return kTRUE;
  else return false5();
}

////////////////////////////////////////////////////////////////////////
Bool_t TT_drawer::false1()
{
  return kFALSE;
}

////////////////////////////////////////////////////////////////////////
Bool_t TT_drawer::false2()
{
  return kFALSE;
}

////////////////////////////////////////////////////////////////////////
Bool_t TT_drawer::false3()
{
  return kFALSE;
}

////////////////////////////////////////////////////////////////////////
Bool_t TT_drawer::false4()
{
  return kFALSE;
}

////////////////////////////////////////////////////////////////////////
Bool_t TT_drawer::false5()
{
  return kFALSE;
}
