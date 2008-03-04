#ifndef TT_surface_h
#define TT_surface_h

#include <iostream>
#include "Rtypes.h"
#include <math.h>
#include "TT_param_class.h"
#include "TT_nucleus.h"
#include "TRandom3.h"
#include "TH1D.h"

using namespace std;

class TT_surface {
 public :

  TT_param_class *f_params;
  TT_nucleus      *f_nucleus;
  Double_t f_Enu;

  TRandom3 f_rand;
  Double_t f_w_min;
  Double_t f_w_max;
  Double_t f_q_bold_min;
  Double_t f_q_bold_max;
  Double_t f_cth_min;
  Double_t f_cth_max;
  Double_t f_phi_min;
  Double_t f_phi_max;

  static const Int_t f_N_dims=4;
  static const Int_t f_N_var=10;
  static const Int_t f_N_w=f_N_var;
  static const Int_t f_N_q_bold=f_N_var;
  static const Int_t f_N_cth=f_N_var;
  static const Int_t f_N_phi=f_N_var;
  Int_t f_N_total;

  //perhaps use the heap
  Double_t f_surface[f_N_w][f_N_q_bold][f_N_cth][f_N_phi];
  Double_t f_volume[f_N_w*f_N_q_bold*f_N_cth*f_N_phi];
  Double_t f_total_volume;

  TH1D *f_volume_histo;
  TH1D *f_volume_check_histo;

  TT_surface(Double_t Enu,TT_param_class *params,TT_nucleus *nucleus,Int_t N_dims);
  virtual ~TT_surface();
  void    Init_surface();
  void    Compute_volume();
  void  Adjust_surface_with_neighbors();
  void  Scale_surface(Double_t);
  Double_t Test_surface(Int_t);
  Double_t Test_surface_randomly(Int_t);
  Double_t Fix_surface(Int_t);
  Double_t Count_surface_zeroes();
  void Draw_point(Double_t&,Double_t&,Double_t&,Double_t&,Double_t&,Double_t&);
  Int_t Index(Int_t,Int_t,Int_t,Int_t);

};

#endif

#ifdef TT_surface_cxx
TT_surface::TT_surface(Double_t Enu,TT_param_class *params,TT_nucleus *nucleus,Int_t N_dims)
{
  f_params=params;
  f_nucleus=nucleus;
  f_Enu=Enu;

  f_rand.SetSeed(0);
  f_w_min=0.0;
  f_w_max=Enu-params->f_m_lep;
  f_q_bold_min=0.0;
  f_q_bold_max=2.0;//4.0*params->f_M*w;//do a better job here
  f_cth_min=-1.0;
  f_cth_max=1.0;
  f_phi_min=-TMath::Pi();
  f_phi_max=TMath::Pi();

  f_volume_histo      =new TH1D("volume_histo"      ,"",f_N_var*f_N_var*f_N_var*f_N_var,-0.5,-0.5+f_N_var*f_N_var*f_N_var*f_N_var);
  f_volume_check_histo=new TH1D("volume_check_histo","",f_N_var*f_N_var*f_N_var*f_N_var,-0.5,-0.5+f_N_var*f_N_var*f_N_var*f_N_var);

}

TT_surface::~TT_surface() {
  delete f_volume_histo;
  delete f_volume_check_histo;
}

#endif // #ifdef TT_surface_cxx

