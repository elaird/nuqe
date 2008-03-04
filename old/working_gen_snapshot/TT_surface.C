#define TT_surface_cxx
#include "TT_surface.h"
#include "TV_class.h"

////////////////////////////////////////////////////////////////////////
void TT_surface::Init_surface()
{
  //sample function for an initial surface
  for (Int_t i_w=0;i_w<f_N_w;i_w++) {
    Double_t w=(i_w+0.5)*(f_w_max-f_w_min)/f_N_w + f_w_min;

    for (Int_t i_q_bold=0;i_q_bold<f_N_q_bold;i_q_bold++) {
      Double_t q_bold=(i_q_bold+0.5)*(f_q_bold_max-f_q_bold_min)/f_N_q_bold + f_q_bold_min;

      for (Int_t i_cth=0;i_cth<f_N_cth;i_cth++) {
	Double_t cth=(i_cth+0.5)*(f_cth_max-f_cth_min)/f_N_cth + f_cth_min;

	for (Int_t i_phi=0;i_phi<f_N_phi;i_phi++) {
	  Double_t phi=(i_phi+0.5)*(f_phi_max-f_phi_min)/f_N_phi + f_phi_min;

	  f_surface[i_w][i_q_bold][i_cth][i_phi]=0.0;	  

	  TV_class foo(f_params,f_nucleus,f_Enu,w,q_bold,cth,phi);
	  foo.Evaluate_d4sigma_dw_dq_bold_domega_p();
	  if (foo.f_d4sigma_dw_dq_bold_domega_p>0.0) {
	    f_surface[i_w][i_q_bold][i_cth][i_phi]=foo.f_d4sigma_dw_dq_bold_domega_p;
	    //cout << "f_s=" << f_surface[i_w][i_q_bold][i_cth][i_phi] << "; foo=" << foo.f_d4sigma_dw_dq_bold_domega_p << endl;
	  }

	}
      }
    }
  }

}

////////////////////////////////////////////////////////////////////////
Double_t TT_surface::Count_surface_zeroes()
{
  Int_t zeroes=0;
  for (Int_t i_w=0;i_w<f_N_w;i_w++) {
    for (Int_t i_q_bold=0;i_q_bold<f_N_q_bold;i_q_bold++) {
      for (Int_t i_cth=0;i_cth<f_N_cth;i_cth++) {
	for (Int_t i_phi=0;i_phi<f_N_phi;i_phi++) {
	  if (f_surface[i_w][i_q_bold][i_cth][i_phi]==0.0) zeroes++;
	  //else cout << f_surface[i_w][i_q_bold][i_cth][i_phi] << ", z:" << zeroes << endl;
	}
      }
    }
  }

  return (zeroes+0.0)/(f_N_w*f_N_q_bold*f_N_cth*f_N_phi);
}

////////////////////////////////////////////////////////////////////////
void TT_surface::Compute_volume()
{

  Double_t factor=1.0;
  factor*=(f_w_max-f_w_min)/f_N_var;
  factor*=(f_q_bold_max-f_q_bold_min)/f_N_var;
  factor*=(f_cth_max-f_cth_min)/f_N_var;
  factor*=(f_phi_max-f_phi_min)/f_N_var;

  Double_t previous_value=0.0;
  for (Int_t i_w=0;i_w<f_N_w;i_w++) {
    for (Int_t i_q_bold=0;i_q_bold<f_N_q_bold;i_q_bold++) {
      for (Int_t i_cth=0;i_cth<f_N_cth;i_cth++) {
	for (Int_t i_phi=0;i_phi<f_N_phi;i_phi++) {

	  Int_t index=Index(i_w,i_q_bold,i_cth,i_phi);
	  f_volume[index]=previous_value;

	  previous_value=f_volume[index]+factor*f_surface[i_w][i_q_bold][i_cth][i_phi];
	  
	  f_volume_histo->SetBinContent(index+1,factor*f_surface[i_w][i_q_bold][i_cth][i_phi]);
	  f_volume_histo->SetBinError(index+1,0.0);

	}
      }
    }
  }

  f_total_volume=previous_value;
}

////////////////////////////////////////////////////////////////////////
void TT_surface::Adjust_surface_with_neighbors()
{

  //ugh--copy surface
  Double_t ini_surface[f_N_w][f_N_q_bold][f_N_cth][f_N_phi];
  for (Int_t i_w=0;i_w<f_N_w;i_w++) {
    for (Int_t i_q_bold=0;i_q_bold<f_N_q_bold;i_q_bold++) {
      for (Int_t i_cth=0;i_cth<f_N_cth;i_cth++) {
	for (Int_t i_phi=0;i_phi<f_N_phi;i_phi++) {
	  ini_surface[i_w][i_q_bold][i_cth][i_phi]=f_surface[i_w][i_q_bold][i_cth][i_phi];
	}
      }
    }
  }

  //now adjust surface
  for (Int_t i_w=0;i_w<f_N_w;i_w++) {
    for (Int_t i_q_bold=0;i_q_bold<f_N_q_bold;i_q_bold++) {
      for (Int_t i_cth=0;i_cth<f_N_cth;i_cth++) {
	for (Int_t i_phi=0;i_phi<f_N_phi;i_phi++) {
	  Double_t neighbors[4][2];
	  
	  if (i_w==0    )             neighbors[0][0]=0.0;
	  else                        neighbors[0][0]=ini_surface[i_w-1][i_q_bold][i_cth][i_phi];
	  if (i_w==f_N_w-1)           neighbors[0][1]=0.0;
	  else                        neighbors[0][1]=ini_surface[i_w+1][i_q_bold][i_cth][i_phi];
	  if (i_q_bold==0         )   neighbors[1][0]=0.0;
	  else                        neighbors[1][0]=ini_surface[i_w][i_q_bold-1][i_cth][i_phi];
	  if (i_q_bold==f_N_q_bold-1) neighbors[1][1]=0.0;
	  else                        neighbors[1][1]=ini_surface[i_w][i_q_bold+1][i_cth][i_phi];
	  if (i_cth==0      )         neighbors[2][0]=0.0;
	  else                        neighbors[2][0]=ini_surface[i_w][i_q_bold][i_cth-1][i_phi];
	  if (i_cth==f_N_cth-1)       neighbors[2][1]=0.0;
	  else                        neighbors[2][1]=ini_surface[i_w][i_q_bold][i_cth+1][i_phi];
	  if (i_phi==0      )         neighbors[3][0]=0.0;
	  else                        neighbors[3][0]=ini_surface[i_w][i_q_bold][i_cth][i_phi-1];
	  if (i_phi==f_N_phi-1)       neighbors[3][1]=0.0;
	  else                        neighbors[3][1]=ini_surface[i_w][i_q_bold][i_cth][i_phi+1];

	  Double_t max=ini_surface[i_w][i_q_bold][i_cth][i_phi];
	  //pick the biggest neighbor
	  for (Int_t i_dim=0;i_dim<f_N_dims;i_dim++) {
	    for (Int_t i_side=0;i_side<2;i_side++) {
	      if (neighbors[i_dim][i_side]>max) max=neighbors[i_dim][i_side];
	    }
	  }
	  f_surface[i_w][i_q_bold][i_cth][i_phi]=max;
	}
      }
    }
  }

}

////////////////////////////////////////////////////////////////////////
Double_t TT_surface::Test_surface(Int_t factor)
{
  Int_t problems=0;
  for (Int_t i_w=0;i_w<factor*f_N_w;i_w++) {
    Double_t w=(i_w+0.5)*(f_w_max-f_w_min)/(factor*f_N_w) + f_w_min;

    for (Int_t i_q_bold=0;i_q_bold<factor*f_N_q_bold;i_q_bold++) {
      Double_t q_bold=(i_q_bold+0.5)*(f_q_bold_max-f_q_bold_min)/(factor*f_N_q_bold) + f_q_bold_min;

      for (Int_t i_cth=0;i_cth<factor*f_N_cth;i_cth++) {
	Double_t cth=(i_cth+0.5)*(f_cth_max-f_cth_min)/(factor*f_N_cth) + f_cth_min;

	for (Int_t i_phi=0;i_phi<factor*f_N_phi;i_phi++) {
	  Double_t phi=(i_phi+0.5)*(f_phi_max-f_phi_min)/(factor*f_N_phi) + f_phi_min;

	  TV_class foo(f_params,f_nucleus,f_Enu,w,q_bold,cth,phi);
	  foo.Evaluate_d4sigma_dw_dq_bold_domega_p();

	  Int_t ind_w=i_w/factor;
	  Int_t ind_q_bold=i_q_bold/factor;
	  Int_t ind_cth=i_cth/factor;
	  Int_t ind_phi=i_phi/factor;

	  if (foo.f_d4sigma_dw_dq_bold_domega_p>f_surface[ind_w][ind_q_bold][ind_cth][ind_phi]) {
	    //cout << "f_s=" << f_surface[ind_w][ind_q_bold][ind_cth][ind_phi] << "; foo=" << foo.f_d4sigma_dw_dq_bold_domega_p << endl;
	    problems++;
	  }

	}
      }
    }
  }
  return (problems+0.0)/pow(f_N_var*factor,f_N_dims);
}

////////////////////////////////////////////////////////////////////////
Double_t TT_surface::Test_surface_randomly(Int_t factor)
{
  Int_t problems=0;

  Int_t max=factor*f_N_var*f_N_var*f_N_var*f_N_var;
  for (Int_t i=0;i<max;i++) {

    Int_t ind_w=f_rand.Integer(f_N_var);
    Int_t ind_q_bold=f_rand.Integer(f_N_var);
    Int_t ind_cth=f_rand.Integer(f_N_var);
    Int_t ind_phi=f_rand.Integer(f_N_var);

    Double_t w_min=f_w_min+ind_w*(f_w_max-f_w_min)/f_N_w;
    Double_t w_max=f_w_min+(ind_w+1)*(f_w_max-f_w_min)/f_N_w;

    Double_t q_bold_min=f_q_bold_min+ ind_q_bold   *(f_q_bold_max-f_q_bold_min)/f_N_q_bold;
    Double_t q_bold_max=f_q_bold_min+(ind_q_bold+1)*(f_q_bold_max-f_q_bold_min)/f_N_q_bold;

    Double_t cth_min=f_cth_min+ ind_cth   *(f_cth_max-f_cth_min)/f_N_cth;
    Double_t cth_max=f_cth_min+(ind_cth+1)*(f_cth_max-f_cth_min)/f_N_cth;

    Double_t phi_min=f_phi_min+ ind_phi   *(f_phi_max-f_phi_min)/f_N_phi;
    Double_t phi_max=f_phi_min+(ind_phi+1)*(f_phi_max-f_phi_min)/f_N_phi;


    Double_t w=f_rand.Uniform(w_min,w_max);
    Double_t q_bold=f_rand.Uniform(q_bold_min,q_bold_max);
    Double_t cth=f_rand.Uniform(cth_min,cth_max);
    Double_t phi=f_rand.Uniform(phi_min,phi_max);

    TV_class foo(f_params,f_nucleus,f_Enu,w,q_bold,cth,phi);
    foo.Evaluate_d4sigma_dw_dq_bold_domega_p();

    if (foo.f_d4sigma_dw_dq_bold_domega_p>f_surface[ind_w][ind_q_bold][ind_cth][ind_phi]) {
      cout << "random test: f_s=" << f_surface[ind_w][ind_q_bold][ind_cth][ind_phi] << "; foo=" << foo.f_d4sigma_dw_dq_bold_domega_p << endl;
      problems++;
    }

  }
  return (problems+0.0)/max;
}

////////////////////////////////////////////////////////////////////////
Double_t TT_surface::Fix_surface(Int_t factor)
{
  Int_t problems=0;
  for (Int_t i_w=0;i_w<factor*f_N_w;i_w++) {
    Double_t w=(i_w+0.5)*(f_w_max-f_w_min)/(factor*f_N_w) + f_w_min;

    for (Int_t i_q_bold=0;i_q_bold<factor*f_N_q_bold;i_q_bold++) {
      Double_t q_bold=(i_q_bold+0.5)*(f_q_bold_max-f_q_bold_min)/(factor*f_N_q_bold) + f_q_bold_min;

      for (Int_t i_cth=0;i_cth<factor*f_N_cth;i_cth++) {
	Double_t cth=(i_cth+0.5)*(f_cth_max-f_cth_min)/(factor*f_N_cth) + f_cth_min;

	for (Int_t i_phi=0;i_phi<factor*f_N_phi;i_phi++) {
	  Double_t phi=(i_phi+0.5)*(f_phi_max-f_phi_min)/(factor*f_N_phi) + f_phi_min;

	  TV_class foo(f_params,f_nucleus,f_Enu,w,q_bold,cth,phi);
	  foo.Evaluate_d4sigma_dw_dq_bold_domega_p();

	  Int_t ind_w=i_w/factor;
	  Int_t ind_q_bold=i_q_bold/factor;
	  Int_t ind_cth=i_cth/factor;
	  Int_t ind_phi=i_phi/factor;

	  if (foo.f_d4sigma_dw_dq_bold_domega_p>f_surface[ind_w][ind_q_bold][ind_cth][ind_phi]) {
	    f_surface[ind_w][ind_q_bold][ind_cth][ind_phi]=2.0*foo.f_d4sigma_dw_dq_bold_domega_p;
	    problems++;
	  }

	}
      }
    }
  }
  return (problems+0.0)/pow(f_N_var*factor,f_N_dims);
}

////////////////////////////////////////////////////////////////////////
void TT_surface::Scale_surface(Double_t scale_factor)
{
  for (Int_t i_w=0;i_w<f_N_w;i_w++) {
    for (Int_t i_q_bold=0;i_q_bold<f_N_q_bold;i_q_bold++) {
      for (Int_t i_cth=0;i_cth<f_N_cth;i_cth++) {
	for (Int_t i_phi=0;i_phi<f_N_phi;i_phi++) {
	  f_surface[i_w][i_q_bold][i_cth][i_phi]*=scale_factor;
	  Int_t index=Index(i_w,i_q_bold,i_cth,i_phi);
	  f_volume[index]*=scale_factor;
	}
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////
Int_t TT_surface::Index(Int_t i_w,Int_t i_q_bold,Int_t i_cth,Int_t i_phi)
{
  Int_t ans=0;
  ans+=i_w     *f_N_var*f_N_var*f_N_var;
  ans+=i_q_bold*f_N_var*f_N_var;
  ans+=i_cth   *f_N_var;
  ans+=i_phi;

  return ans;
}

////////////////////////////////////////////////////////////////////////
void TT_surface::Draw_point(Double_t& w,Double_t& q_bold,Double_t& cth,Double_t& phi,Double_t& xs,Double_t& inv_weight)
{
  Int_t num=f_N_var*f_N_var*f_N_var*f_N_var;
  Double_t x=f_rand.Uniform(0.0,f_total_volume);
  Int_t index=TMath::BinarySearch(num,f_volume,x);
  //printf("index=%8.6d; x=%8.6g; vl=%8.6g; vr=%8.6g\n",index,x,f_volume[index],f_volume[index+1]);
  f_volume_check_histo->Fill(index,1.0);
  Int_t index_orig=index;

  Int_t i_phi=index%f_N_var;
  index-=i_phi;

  Int_t i_cth=(index/f_N_var)%f_N_var;
  index-=f_N_var*i_cth;

  Int_t i_q_bold=(index/(f_N_var*f_N_var))%f_N_var;
  index-=f_N_var*f_N_var*i_q_bold;

  Int_t i_w=(index/(f_N_var*f_N_var*f_N_var))%f_N_var;
  index-=f_N_var*f_N_var*f_N_var*i_w;

  if (index_orig!=Index(i_w,i_q_bold,i_cth,i_phi)) printf("uh oh: index=%d; (i_w;i_q_bold;i_cth;i_phi)=(%d;%d;%d;%d)\n",index_orig,i_w,i_q_bold,i_cth,i_phi);

  Double_t w_min= i_w   *(f_w_max-f_w_min)/f_N_w+f_w_min;
  Double_t w_max=(i_w+1)*(f_w_max-f_w_min)/f_N_w+f_w_min;

  Double_t q_bold_min= i_q_bold   *(f_q_bold_max-f_q_bold_min)/f_N_q_bold+f_q_bold_min;
  Double_t q_bold_max=(i_q_bold+1)*(f_q_bold_max-f_q_bold_min)/f_N_q_bold+f_q_bold_min;

  Double_t cth_min= i_cth   *(f_cth_max-f_cth_min)/f_N_cth+f_cth_min;
  Double_t cth_max=(i_cth+1)*(f_cth_max-f_cth_min)/f_N_cth+f_cth_min;

  Double_t phi_min= i_phi   *(f_phi_max-f_phi_min)/f_N_phi+f_phi_min;
  Double_t phi_max=(i_phi+1)*(f_phi_max-f_phi_min)/f_N_phi+f_phi_min;

  w=f_rand.Uniform(w_min,w_max);
  q_bold=f_rand.Uniform(q_bold_min,q_bold_max);
  cth=f_rand.Uniform(cth_min,cth_max);
  phi=f_rand.Uniform(phi_min,phi_max);

  TV_class foo(f_params,f_nucleus,f_Enu,w,q_bold,cth,phi);
  foo.Evaluate_d4sigma_dw_dq_bold_domega_p();

  //rejection method
  Double_t y=f_rand.Uniform(0.0,f_surface[i_w][i_q_bold][i_cth][i_phi]);

  inv_weight=1.0;
  if (f_surface[i_w][i_q_bold][i_cth][i_phi]<foo.f_d4sigma_dw_dq_bold_domega_p) {
    //cout << "broken surface=" << Index(i_w,i_q_bold,i_cth,i_phi) << "," << f_surface[i_w][i_q_bold][i_cth][i_phi]/foo.f_d4sigma_dw_dq_bold_domega_p << endl;
    inv_weight=f_surface[i_w][i_q_bold][i_cth][i_phi]/foo.f_d4sigma_dw_dq_bold_domega_p;
  }

  xs=0.0;
  if (y<foo.f_d4sigma_dw_dq_bold_domega_p) xs=foo.f_d4sigma_dw_dq_bold_domega_p;
}

