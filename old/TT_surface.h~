#ifndef TT_surface_h
#define TT_surface_h

#include <iostream>
#include "Rtypes.h"
#include <math.h>
#include "TT_params.h"
#include "TT_nucleus.h"
#include "TV_class.h"
#include "TRandom3.h"
#include "TH1D.h"

using namespace std;

class TT_surface {
 public :

  TT_params  *f_params;
  TT_nucleus *f_nucleus;
  Double_t f_Enu;

  TRandom3 f_rand;
  Double_t f_w_min;
  Double_t f_w_max;
  Double_t f_cth_p_min;
  Double_t f_cth_p_max;
  Double_t f_phi_p_min;
  Double_t f_phi_p_max;

  Int_t f_N_dims;
  Int_t f_N_var;
  Int_t f_N_w;
  Int_t f_N_qbold;
  Int_t f_N_cth_p;
  Int_t f_N_phi_p;

  //Double_t f_surface[f_N_w][f_N_qbold][f_N_cth_p][f_N_phi_p];
  //Double_t f_integral[f_N_w*f_N_qbold*f_N_cth_p*f_N_phi_p];
  Double_t ****f_surface;
  Double_t *f_integral;
  Double_t f_total_integral;
  Double_t f_total_volume;

  TH1D *f_integral_histo;
  TH1D *f_integral_check_histo;

  TT_surface(Double_t Enu,TT_params *params,TT_nucleus *nucleus,Int_t N_dims,Int_t N_var);
  virtual ~TT_surface();
  void    Init_surface();
  Double_t Get_volume_factor(Int_t i_w,Int_t i_qbold,Int_t i_cth_p,Int_t i_phi_p);
  void    Compute_volume_and_integral();
  void  Adjust_surface_with_neighbors();
  void  Adjust_surface_with_neighbors2();
  void  Scale_surface(Double_t);
  void Replace_low_bins();
  Double_t Test_surface(Int_t);
  Double_t Test_surface_randomly(Int_t);
  Double_t Fix_surface(Int_t);
  Double_t Count_surface_zeroes();
  Bool_t Draw_point(TV_class*,Double_t&);
  Int_t Index(Int_t,Int_t,Int_t,Int_t);

  Double_t Use_Minuit(Int_t);

};

#endif

#ifdef TT_surface_cxx
TT_surface::TT_surface(Double_t Enu,TT_params *params,TT_nucleus *nucleus,Int_t N_dims,Int_t N_var)
{
  f_params=params;
  f_nucleus=nucleus;
  f_Enu=Enu;
  f_N_dims=N_dims;
  f_N_var=N_var;
  f_N_w=f_N_var;
  f_N_qbold=f_N_var;
  f_N_cth_p=f_N_var;
  f_N_phi_p=f_N_var;

  f_integral=new Double_t[f_N_w*f_N_qbold*f_N_cth_p*f_N_phi_p];

  f_surface=new Double_t***[f_N_w];
  for (Int_t i_w=0;i_w<f_N_w;i_w++) {
    f_surface[i_w]=new Double_t**[f_N_qbold];
    for (Int_t i_qbold=0;i_qbold<f_N_qbold;i_qbold++) {
      f_surface[i_w][i_qbold]=new Double_t*[f_N_cth_p];
      for (Int_t i_cth_p=0;i_cth_p<f_N_cth_p;i_cth_p++) {
	f_surface[i_w][i_qbold][i_cth_p]=new Double_t[f_N_phi_p];
      }
    }
  }
  
  f_rand.SetSeed(0);
  f_w_min=0.0;
  f_w_max=Enu-params->f_m_lep;
  f_cth_p_min=-1.0;
  f_cth_p_max=1.0;
  f_phi_p_min=-TMath::Pi();
  f_phi_p_max=TMath::Pi();

  Int_t nbins=f_N_var*f_N_var*f_N_var*f_N_var;
  f_integral_histo      =new TH1D("integral_histo"      ,"",nbins,-0.5,-0.5+nbins);
  f_integral_check_histo=new TH1D("integral_check_histo","",nbins,-0.5,-0.5+nbins);

}

TT_surface::~TT_surface() {
  for (Int_t i_w=0;i_w<f_N_w;i_w++) {
    for (Int_t i_qbold=0;i_qbold<f_N_qbold;i_qbold++) {
      for (Int_t i_cth_p=0;i_cth_p<f_N_cth_p;i_cth_p++) {
	delete f_surface[i_w][i_qbold][i_cth_p];
      }
      delete f_surface[i_w][i_qbold];
    }
    delete f_surface[i_w];
  }
  delete f_surface;
  
  delete f_integral;

  delete f_integral_histo;
  delete f_integral_check_histo;
}

#endif // #ifdef TT_surface_cxx

