#define TT_surface_cxx
#include "TT_surface.h"
#include "TT_event.h"
#include "TMinuit.h"

void Minuit_FCN(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

////////////////////////////////////////////////////////////////////////
void TT_surface::Init_surface()
{
  //sample function for an initial surface
  for (Int_t i_w=0;i_w<f_N_w;i_w++) {
    Double_t w =(i_w+0.5)*(f_w_max-f_w_min)/f_N_w + f_w_min;
    Double_t ws=(i_w)    *(f_w_max-f_w_min)/f_N_w + f_w_min;

    Double_t p_lep=sqrt(pow(f_Enu-ws,2)-pow(f_params->f_m_lep,2));
    Double_t ws_qbold_min=f_Enu-p_lep;
    Double_t ws_qbold_max=f_Enu+p_lep;

    for (Int_t i_qbold=0;i_qbold<f_N_qbold;i_qbold++) {
      Double_t qbold=(i_qbold+0.5)*(ws_qbold_max-ws_qbold_min)/f_N_qbold + ws_qbold_min;

      for (Int_t i_mag_p=0;i_mag_p<f_N_mag_p;i_mag_p++) {
	Double_t mag_p=(i_mag_p+0.5)*(f_mag_p_max-f_mag_p_min)/f_N_mag_p + f_mag_p_min;

	for (Int_t i_phi_p=0;i_phi_p<f_N_phi_p;i_phi_p++) {
	  Double_t phi_p=(i_phi_p+0.5)*(f_phi_p_max-f_phi_p_min)/f_N_phi_p + f_phi_p_min;

	  f_surface[i_w][i_qbold][i_mag_p][i_phi_p]=0.0;	  

	  TT_event foo;
	  foo.Init(f_params,f_nucleus,f_Enu,w,qbold,mag_p,phi_p);
	  foo.Evaluate_dsigma_dall();
	  if (foo.f_dsigma_dall>0.0) {
	    f_surface[i_w][i_qbold][i_mag_p][i_phi_p]=foo.f_dsigma_dall;
	  }
	  //cout << "index=" << Index(i_w,i_qbold,i_mag_p,i_phi_p) << "; xs=" << foo.f_dsigma_dall << endl;

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
    for (Int_t i_qbold=0;i_qbold<f_N_qbold;i_qbold++) {
      for (Int_t i_mag_p=0;i_mag_p<f_N_mag_p;i_mag_p++) {
	for (Int_t i_phi_p=0;i_phi_p<f_N_phi_p;i_phi_p++) {
	  if (f_surface[i_w][i_qbold][i_mag_p][i_phi_p]==0.0) zeroes++;
	  //else cout << f_surface[i_w][i_qbold][i_mag_p][i_phi_p] << ", z:" << zeroes << endl;
	}
      }
    }
  }

  return (zeroes+0.0)/(f_N_w*f_N_qbold*f_N_mag_p*f_N_phi_p);
}

////////////////////////////////////////////////////////////////////////
Double_t TT_surface::Get_volume_factor(Int_t i_w,Int_t i_qbold,Int_t i_mag_p,Int_t i_phi_p)
{
  Double_t ws=(i_w)    *(f_w_max-f_w_min)/f_N_w + f_w_min;
  Double_t p_lep=sqrt(pow(f_Enu-ws,2)-pow(f_params->f_m_lep,2));
  Double_t ws_qbold_min=f_Enu-p_lep;
  Double_t ws_qbold_max=f_Enu+p_lep;

  Double_t factor=1.0;
  factor*=(f_w_max-f_w_min)/f_N_var;
  factor*=(ws_qbold_max-ws_qbold_min)/f_N_var;
  factor*=(f_mag_p_max-f_mag_p_min)/f_N_var;
  factor*=(f_phi_p_max-f_phi_p_min)/f_N_var;
  return factor;
}

////////////////////////////////////////////////////////////////////////
void TT_surface::Compute_volume_and_integral()
{
  Double_t previous_value=0.0;
  f_total_volume=0.0;
  for (Int_t i_w=0;i_w<f_N_w;i_w++) {
    for (Int_t i_qbold=0;i_qbold<f_N_qbold;i_qbold++) {
      for (Int_t i_mag_p=0;i_mag_p<f_N_mag_p;i_mag_p++) {
	for (Int_t i_phi_p=0;i_phi_p<f_N_phi_p;i_phi_p++) {

	  Int_t index=Index(i_w,i_qbold,i_mag_p,i_phi_p);
	  f_integral[index]=previous_value;

	  Double_t factor=Get_volume_factor(i_w,i_qbold,i_mag_p,i_phi_p);
	  f_total_volume+=factor;

	  previous_value=f_integral[index]+factor*f_surface[i_w][i_qbold][i_mag_p][i_phi_p];
	  
	  f_integral_histo->SetBinContent(index+1,factor*f_surface[i_w][i_qbold][i_mag_p][i_phi_p]);
	  f_integral_histo->SetBinError(index+1,0.0);

	}
      }
    }
  }

  f_total_integral=previous_value;
}

////////////////////////////////////////////////////////////////////////
void TT_surface::Adjust_surface_with_neighbors()
{

  //ugh--copy surface
  Double_t ini_surface[f_N_w][f_N_qbold][f_N_mag_p][f_N_phi_p];
  for (Int_t i_w=0;i_w<f_N_w;i_w++) {
    for (Int_t i_qbold=0;i_qbold<f_N_qbold;i_qbold++) {
      for (Int_t i_mag_p=0;i_mag_p<f_N_mag_p;i_mag_p++) {
	for (Int_t i_phi_p=0;i_phi_p<f_N_phi_p;i_phi_p++) {
	  ini_surface[i_w][i_qbold][i_mag_p][i_phi_p]=f_surface[i_w][i_qbold][i_mag_p][i_phi_p];
	}
      }
    }
  }

  //now adjust surface
  for (Int_t i_w=0;i_w<f_N_w;i_w++) {
    for (Int_t i_qbold=0;i_qbold<f_N_qbold;i_qbold++) {
      for (Int_t i_mag_p=0;i_mag_p<f_N_mag_p;i_mag_p++) {
	for (Int_t i_phi_p=0;i_phi_p<f_N_phi_p;i_phi_p++) {
	  Double_t neighbors[4][2];
	  
	  if (i_w==0    )             neighbors[0][0]=0.0;
	  else                        neighbors[0][0]=ini_surface[i_w-1][i_qbold][i_mag_p][i_phi_p];
	  if (i_w==f_N_w-1)           neighbors[0][1]=0.0;
	  else                        neighbors[0][1]=ini_surface[i_w+1][i_qbold][i_mag_p][i_phi_p];
	  if (i_qbold==0         )    neighbors[1][0]=0.0;
	  else                        neighbors[1][0]=ini_surface[i_w][i_qbold-1][i_mag_p][i_phi_p];
	  if (i_qbold==f_N_qbold-1)   neighbors[1][1]=0.0;
	  else                        neighbors[1][1]=ini_surface[i_w][i_qbold+1][i_mag_p][i_phi_p];
	  if (i_mag_p==0      )       neighbors[2][0]=0.0;
	  else                        neighbors[2][0]=ini_surface[i_w][i_qbold][i_mag_p-1][i_phi_p];
	  if (i_mag_p==f_N_mag_p-1)   neighbors[2][1]=0.0;
	  else                        neighbors[2][1]=ini_surface[i_w][i_qbold][i_mag_p+1][i_phi_p];
	  if (i_phi_p==0      )       neighbors[3][0]=0.0;
	  else                        neighbors[3][0]=ini_surface[i_w][i_qbold][i_mag_p][i_phi_p-1];
	  if (i_phi_p==f_N_phi_p-1)   neighbors[3][1]=0.0;
	  else                        neighbors[3][1]=ini_surface[i_w][i_qbold][i_mag_p][i_phi_p+1];

	  Double_t max=ini_surface[i_w][i_qbold][i_mag_p][i_phi_p];
	  //pick the biggest neighbor
	  for (Int_t i_dim=0;i_dim<f_N_dims;i_dim++) {
	    for (Int_t i_side=0;i_side<2;i_side++) {
	      if (neighbors[i_dim][i_side]>max) max=neighbors[i_dim][i_side];
	    }
	  }
	  f_surface[i_w][i_qbold][i_mag_p][i_phi_p]=max;
	}
      }
    }
  }

}

////////////////////////////////////////////////////////////////////////
void TT_surface::Adjust_surface_with_neighbors2()
{

  //ugh--copy surface
  Double_t ini_surface[f_N_w][f_N_qbold][f_N_mag_p][f_N_phi_p];
  for (Int_t i_w=0;i_w<f_N_w;i_w++) {
    for (Int_t i_qbold=0;i_qbold<f_N_qbold;i_qbold++) {
      for (Int_t i_mag_p=0;i_mag_p<f_N_mag_p;i_mag_p++) {
	for (Int_t i_phi_p=0;i_phi_p<f_N_phi_p;i_phi_p++) {
	  ini_surface[i_w][i_qbold][i_mag_p][i_phi_p]=f_surface[i_w][i_qbold][i_mag_p][i_phi_p];
	}
      }
    }
  }

  //now adjust surface
  for (Int_t i_w=0;i_w<f_N_w;i_w++) {
    for (Int_t i_qbold=0;i_qbold<f_N_qbold;i_qbold++) {
      for (Int_t i_mag_p=0;i_mag_p<f_N_mag_p;i_mag_p++) {
	for (Int_t i_phi_p=0;i_phi_p<f_N_phi_p;i_phi_p++) {

	  //ni=neighboring_indices-- not necessarily physical neighbors
	  Int_t ni[f_N_dims][3];
	  Int_t potential_index;
	  for (Int_t side=0;side<3;side++) {
	    potential_index=i_phi_p-1+side;
	    ni[3][side]=i_phi_p;
	    if (potential_index>=0 && potential_index<f_N_phi_p) ni[3][side]=potential_index;

	    potential_index=i_mag_p-1+side;
	    ni[2][side]=i_mag_p;
	    if (potential_index>=0 && potential_index<f_N_mag_p) ni[2][side]=potential_index;

	    potential_index=i_qbold-1+side;
	    ni[1][side]=i_qbold;
	    if (potential_index>=0 && potential_index<f_N_qbold) ni[1][side]=potential_index;

	    potential_index=i_w-1+side;
	    ni[0][side]=i_w;
	    if (potential_index>=0 && potential_index<f_N_w) ni[0][side]=potential_index;
	  }

	  //pick the biggest neighbor
	  Double_t max=ini_surface[i_w][i_qbold][i_mag_p][i_phi_p];
	  for (Int_t side0=0;side0<3;side0++) {
	    for (Int_t side1=0;side1<3;side1++) {
	      for (Int_t side2=0;side2<3;side2++) {
		for (Int_t side3=0;side3<3;side3++) {
		  Int_t ind0=ni[0][side0];
		  Int_t ind1=ni[1][side1];
		  Int_t ind2=ni[2][side2];
		  Int_t ind3=ni[3][side3];
		  Double_t xs=ini_surface[ind0][ind1][ind2][ind3];
		  if (xs>max) max=xs;
		}
	      }
	    }
	  } //done picking the biggest neighbor

	  f_surface[i_w][i_qbold][i_mag_p][i_phi_p]=max;
	}
      }
    }
  }

}

////////////////////////////////////////////////////////////////////////
Double_t TT_surface::Find_max()
{
  Double_t ans=0.0;
  for (Int_t i_w=0;i_w<f_N_w;i_w++) {
    for (Int_t i_qbold=0;i_qbold<f_N_qbold;i_qbold++) {
      for (Int_t i_mag_p=0;i_mag_p<f_N_mag_p;i_mag_p++) {
  	for (Int_t i_phi_p=0;i_phi_p<f_N_phi_p;i_phi_p++) {
	  if (f_surface[i_w][i_qbold][i_mag_p][i_phi_p]>ans) ans=f_surface[i_w][i_qbold][i_mag_p][i_phi_p];
  	}
      }
    }
  }
  return ans;
}

////////////////////////////////////////////////////////////////////////
void TT_surface::Install_value(Double_t value)
{
  for (Int_t i_w=0;i_w<f_N_w;i_w++) {
    for (Int_t i_qbold=0;i_qbold<f_N_qbold;i_qbold++) {
      for (Int_t i_mag_p=0;i_mag_p<f_N_mag_p;i_mag_p++) {
  	for (Int_t i_phi_p=0;i_phi_p<f_N_phi_p;i_phi_p++) {
	  f_surface[i_w][i_qbold][i_mag_p][i_phi_p]=value;
  	}
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////
void TT_surface::Replace_low_bins()
{
  Double_t content=f_total_integral/f_total_volume;
  for (Int_t i_w=0;i_w<f_N_w;i_w++) {
    for (Int_t i_qbold=0;i_qbold<f_N_qbold;i_qbold++) {
      for (Int_t i_mag_p=0;i_mag_p<f_N_mag_p;i_mag_p++) {
  	for (Int_t i_phi_p=0;i_phi_p<f_N_phi_p;i_phi_p++) {
  	  if (f_surface[i_w][i_qbold][i_mag_p][i_phi_p]<content) {
  	    f_surface[i_w][i_qbold][i_mag_p][i_phi_p]=content;
  	  }
  	}
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////
Double_t TT_surface::Use_Minuit(Int_t factor)
{
  Int_t problems=0;
  for (Int_t i_w=0;i_w<f_N_w;i_w++) {
    Double_t w_min=f_w_min+i_w*(f_w_max-f_w_min)/f_N_w;
    Double_t w_max=f_w_min+(i_w+1)*(f_w_max-f_w_min)/f_N_w;

    Double_t ws=(i_w)    *(f_w_max-f_w_min)/f_N_w + f_w_min;
    Double_t p_lep=sqrt(pow(f_Enu-ws,2)-pow(f_params->f_m_lep,2));
    Double_t ws_qbold_min=f_Enu-p_lep;
    Double_t ws_qbold_max=f_Enu+p_lep;

    for (Int_t i_qbold=0;i_qbold<f_N_qbold;i_qbold++) {
      for (Int_t i_mag_p=0;i_mag_p<f_N_mag_p;i_mag_p++) {
	for (Int_t i_phi_p=0;i_phi_p<f_N_phi_p;i_phi_p++) {

	  Double_t qbold_min=ws_qbold_min+ i_qbold   *(ws_qbold_max-ws_qbold_min)/f_N_qbold;
	  Double_t qbold_max=ws_qbold_min+(i_qbold+1)*(ws_qbold_max-ws_qbold_min)/f_N_qbold;

	  Double_t mag_p_min=f_mag_p_min+ i_mag_p   *(f_mag_p_max-f_mag_p_min)/f_N_mag_p;
	  Double_t mag_p_max=f_mag_p_min+(i_mag_p+1)*(f_mag_p_max-f_mag_p_min)/f_N_mag_p;

	  Double_t phi_p_min=f_phi_p_min+ i_phi_p   *(f_phi_p_max-f_phi_p_min)/f_N_phi_p;
	  Double_t phi_p_max=f_phi_p_min+(i_phi_p+1)*(f_phi_p_max-f_phi_p_min)/f_N_phi_p;

	  for (Int_t trial=0;trial<factor;trial++) {

	    //static?
	    Double_t vstart[4];
	    Double_t step[4];
	    Double_t lower[4];
	    Double_t upper[4];

	    vstart[0]=0.0922963; 
	    vstart[1]=0.3128; 
	    vstart[2]=0.470622; 
	    vstart[3]=0.715178;
	    
	    //vstart[0]=f_rand.Uniform(w_min,w_max);
	    //vstart[1]=f_rand.Uniform(qbold_min,qbold_max);
	    //vstart[2]=f_rand.Uniform(mag_p_min,mag_p_max);
	    //vstart[3]=f_rand.Uniform(phi_p_min,phi_p_max);
	    
	    lower[0]=w_min;
	    upper[0]=w_max;
	    lower[1]=qbold_min;
	    upper[1]=qbold_max;
	    lower[2]=mag_p_min;
	    upper[2]=mag_p_max;
	    lower[3]=phi_p_min;
	    upper[3]=phi_p_max;

	    //vstart[0]=0.5; 
	    //vstart[1]=0.5;
	    //vstart[2]=0.5; 
	    //vstart[3]=0.5;
	    //
	    //lower[0]=0.0;
	    //upper[0]=1.0;
	    //lower[1]=0.0;
	    //upper[1]=1.0;
	    //lower[2]=0.0;
	    //upper[2]=1.0;
	    //lower[3]=0.0;
	    //upper[3]=1.0;

	    for (Int_t iStep=0;iStep<4;iStep++) {
	      step[iStep]=(upper[iStep]-lower[iStep])/2.0;
	      //printf("step[%d]=%8.6g\n",iStep,step[iStep]);
	    }
	    
	    TMinuit minuit(4);
	    minuit.SetFCN(Minuit_FCN);

	    Double_t arglist[10];
	    Int_t ierflg = 0;
	    
	    arglist[0]=-1;
	    minuit.mnexcm("SET PRI",arglist,1,ierflg);

	    arglist[0] = 1.0;
	    minuit.mnexcm("SET ERR", arglist ,1,ierflg);

	    minuit.mnparm(0, "w",      vstart[0],step[0],lower[0],upper[0],ierflg);
	    minuit.mnparm(1, "q_bold", vstart[1],step[1],lower[1],upper[1],ierflg);
	    minuit.mnparm(2, "mag_p",  vstart[2],step[2],lower[2],upper[2],ierflg);
	    minuit.mnparm(3, "phi_p",  vstart[3],step[3],lower[3],upper[3],ierflg);
	    
	    // Now ready for minimization step
	    arglist[0] = 500;
	    arglist[1] = 1.0e-1;
	    minuit.mnexcm("MIGRAD", arglist ,2,ierflg);
	    
	    // Print results
	    Double_t amin,edm,errdef;
	    Int_t nvpar,nparx,icstat;
	    minuit.mnstat(amin,edm,errdef,nvpar,nparx,icstat);
	    minuit.mnprin(3,amin);


	  }
	}
      }
    }
  }
  return (problems+0.0)/pow(f_N_var*factor,f_N_dims);
}

////////////////////////////////////////////////////////////////////////
Double_t TT_surface::Test_surface(Int_t factor)
{
  Int_t problems=0;
  for (Int_t i_w=0;i_w<factor*f_N_w;i_w++) {
    Double_t w =(i_w+0.5)*(f_w_max-f_w_min)/(factor*f_N_w) + f_w_min;
    Double_t ws=(i_w)    *(f_w_max-f_w_min)/(factor*f_N_w) + f_w_min;

    Double_t p_lep=sqrt(pow(f_Enu-ws,2)-pow(f_params->f_m_lep,2));
    Double_t ws_qbold_min=f_Enu-p_lep;
    Double_t ws_qbold_max=f_Enu+p_lep;

    for (Int_t i_qbold=0;i_qbold<factor*f_N_qbold;i_qbold++) {
      Double_t qbold=(i_qbold+0.5)*(ws_qbold_max-ws_qbold_min)/(factor*f_N_qbold) + ws_qbold_min;

      for (Int_t i_mag_p=0;i_mag_p<factor*f_N_mag_p;i_mag_p++) {
	Double_t mag_p=(i_mag_p+0.5)*(f_mag_p_max-f_mag_p_min)/(factor*f_N_mag_p) + f_mag_p_min;

	for (Int_t i_phi_p=0;i_phi_p<factor*f_N_phi_p;i_phi_p++) {
	  Double_t phi_p=(i_phi_p+0.5)*(f_phi_p_max-f_phi_p_min)/(factor*f_N_phi_p) + f_phi_p_min;

	  TT_event foo;
	  foo.Init(f_params,f_nucleus,f_Enu,w,qbold,mag_p,phi_p);
	  foo.Evaluate_dsigma_dall();

	  Int_t ind_w=i_w/factor;
	  Int_t ind_qbold=i_qbold/factor;
	  Int_t ind_mag_p=i_mag_p/factor;
	  Int_t ind_phi_p=i_phi_p/factor;

	  if (foo.f_dsigma_dall>f_surface[ind_w][ind_qbold][ind_mag_p][ind_phi_p]) {
	    //cout << "f_s=" << f_surface[ind_w][ind_qbold][ind_mag_p][ind_phi_p] << "; foo=" << foo.f_dsigma_dall << endl;
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
    Int_t ind_qbold=f_rand.Integer(f_N_var);
    Int_t ind_mag_p=f_rand.Integer(f_N_var);
    Int_t ind_phi_p=f_rand.Integer(f_N_var);

    Double_t w_min=f_w_min+ind_w*(f_w_max-f_w_min)/f_N_w;
    Double_t w_max=f_w_min+(ind_w+1)*(f_w_max-f_w_min)/f_N_w;

    Double_t ws=(ind_w)    *(f_w_max-f_w_min)/(factor*f_N_w) + f_w_min;

    Double_t p_lep=sqrt(pow(f_Enu-ws,2)-pow(f_params->f_m_lep,2));
    Double_t ws_qbold_min=f_Enu-p_lep;
    Double_t ws_qbold_max=f_Enu+p_lep;

    Double_t qbold_min=ws_qbold_min+ ind_qbold   *(ws_qbold_max-ws_qbold_min)/f_N_qbold;
    Double_t qbold_max=ws_qbold_min+(ind_qbold+1)*(ws_qbold_max-ws_qbold_min)/f_N_qbold;

    Double_t mag_p_min=f_mag_p_min+ ind_mag_p   *(f_mag_p_max-f_mag_p_min)/f_N_mag_p;
    Double_t mag_p_max=f_mag_p_min+(ind_mag_p+1)*(f_mag_p_max-f_mag_p_min)/f_N_mag_p;

    Double_t phi_p_min=f_phi_p_min+ ind_phi_p   *(f_phi_p_max-f_phi_p_min)/f_N_phi_p;
    Double_t phi_p_max=f_phi_p_min+(ind_phi_p+1)*(f_phi_p_max-f_phi_p_min)/f_N_phi_p;


    Double_t w=f_rand.Uniform(w_min,w_max);
    Double_t qbold=f_rand.Uniform(qbold_min,qbold_max);
    Double_t mag_p=f_rand.Uniform(mag_p_min,mag_p_max);
    Double_t phi_p=f_rand.Uniform(phi_p_min,phi_p_max);

    TT_event foo;
    foo.Init(f_params,f_nucleus,f_Enu,w,qbold,mag_p,phi_p);
    foo.Evaluate_dsigma_dall();

    if (foo.f_dsigma_dall>f_surface[ind_w][ind_qbold][ind_mag_p][ind_phi_p]) {
      printf("random test: w,q,c,p=%d%d%d%d; surf=%8.6g; xs=%8.6g\n",ind_w,ind_qbold,ind_mag_p,ind_phi_p,f_surface[ind_w][ind_qbold][ind_mag_p][ind_phi_p],foo.f_dsigma_dall);
      f_surface[ind_w][ind_qbold][ind_mag_p][ind_phi_p]=foo.f_dsigma_dall;
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
    Double_t w =(i_w+0.5)*(f_w_max-f_w_min)/(factor*f_N_w) + f_w_min;
    Double_t ws=(i_w)    *(f_w_max-f_w_min)/(factor*f_N_w) + f_w_min;

    Double_t p_lep=sqrt(pow(f_Enu-ws,2)-pow(f_params->f_m_lep,2));
    Double_t ws_qbold_min=f_Enu-p_lep;
    Double_t ws_qbold_max=f_Enu+p_lep;

    for (Int_t i_qbold=0;i_qbold<factor*f_N_qbold;i_qbold++) {
      Double_t qbold=(i_qbold+0.5)*(ws_qbold_max-ws_qbold_min)/(factor*f_N_qbold) + ws_qbold_min;

      for (Int_t i_mag_p=0;i_mag_p<factor*f_N_mag_p;i_mag_p++) {
	Double_t mag_p=(i_mag_p+0.5)*(f_mag_p_max-f_mag_p_min)/(factor*f_N_mag_p) + f_mag_p_min;

	for (Int_t i_phi_p=0;i_phi_p<factor*f_N_phi_p;i_phi_p++) {
	  Double_t phi_p=(i_phi_p+0.5)*(f_phi_p_max-f_phi_p_min)/(factor*f_N_phi_p) + f_phi_p_min;

	  //printf("%8.6g; %8.6g; %8.6g; %8.6g\n",w,qbold,mag_p,phi_p);
	  TT_event foo;
	  foo.Init(f_params,f_nucleus,f_Enu,w,qbold,mag_p,phi_p);
	  foo.Evaluate_dsigma_dall();

	  Int_t ind_w=i_w/factor;
	  Int_t ind_qbold=i_qbold/factor;
	  Int_t ind_mag_p=i_mag_p/factor;
	  Int_t ind_phi_p=i_phi_p/factor;

	  if (foo.f_dsigma_dall>f_surface[ind_w][ind_qbold][ind_mag_p][ind_phi_p]) {
	    f_surface[ind_w][ind_qbold][ind_mag_p][ind_phi_p]=2.0*foo.f_dsigma_dall;
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
    for (Int_t i_qbold=0;i_qbold<f_N_qbold;i_qbold++) {
      for (Int_t i_mag_p=0;i_mag_p<f_N_mag_p;i_mag_p++) {
	for (Int_t i_phi_p=0;i_phi_p<f_N_phi_p;i_phi_p++) {
	  f_surface[i_w][i_qbold][i_mag_p][i_phi_p]*=scale_factor;
	  Int_t index=Index(i_w,i_qbold,i_mag_p,i_phi_p);
	  f_integral[index]*=scale_factor;
	}
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////
Int_t TT_surface::Index(Int_t i_w,Int_t i_qbold,Int_t i_mag_p,Int_t i_phi_p)
{
  Int_t ans=0;
  ans+=i_w     *f_N_var*f_N_var*f_N_var;
  ans+=i_qbold*f_N_var*f_N_var;
  ans+=i_mag_p   *f_N_var;
  ans+=i_phi_p;

  return ans;
}

////////////////////////////////////////////////////////////////////////
Bool_t TT_surface::Draw_point(TT_event *event,Double_t& weight)
{
  Int_t num=f_N_var*f_N_var*f_N_var*f_N_var;
  Double_t x=f_rand.Uniform(0.0,f_total_integral);
  Int_t index=TMath::BinarySearch(num,f_integral,x);
  //printf("index=%8.6d; x=%8.6g; vl=%8.6g; vr=%8.6g\n",index,x,f_integral[index],f_integral[index+1]);
  f_integral_check_histo->Fill(index,1.0);
  Int_t index_orig=index;

  Int_t i_phi_p=index%f_N_var;
  index-=i_phi_p;

  Int_t i_mag_p=(index/f_N_var)%f_N_var;
  index-=f_N_var*i_mag_p;

  Int_t i_qbold=(index/(f_N_var*f_N_var))%f_N_var;
  index-=f_N_var*f_N_var*i_qbold;

  Int_t i_w=(index/(f_N_var*f_N_var*f_N_var))%f_N_var;
  index-=f_N_var*f_N_var*f_N_var*i_w;

  if (index_orig!=Index(i_w,i_qbold,i_mag_p,i_phi_p)) printf("uh oh: index=%d; (i_w;i_qbold;i_mag_p;i_phi_p)=(%d;%d;%d;%d)\n",index_orig,i_w,i_qbold,i_mag_p,i_phi_p);

  Double_t w_min= i_w   *(f_w_max-f_w_min)/f_N_w+f_w_min;
  Double_t w_max=(i_w+1)*(f_w_max-f_w_min)/f_N_w+f_w_min;

  Double_t p_lep=sqrt(pow(f_Enu-w_min,2)-pow(f_params->f_m_lep,2));
  Double_t qbold_min=f_Enu-p_lep;
  Double_t qbold_max=f_Enu+p_lep;

  //Double_t qbold_min= i_qbold   *(ws_qbold_max-ws_qbold_min)/f_N_qbold+ws_qbold_min;
  //Double_t qbold_max=(i_qbold+1)*(ws_qbold_max-ws_qbold_min)/f_N_qbold+ws_qbold_min;

  Double_t mag_p_min= i_mag_p   *(f_mag_p_max-f_mag_p_min)/f_N_mag_p+f_mag_p_min;
  Double_t mag_p_max=(i_mag_p+1)*(f_mag_p_max-f_mag_p_min)/f_N_mag_p+f_mag_p_min;

  Double_t phi_p_min= i_phi_p   *(f_phi_p_max-f_phi_p_min)/f_N_phi_p+f_phi_p_min;
  Double_t phi_p_max=(i_phi_p+1)*(f_phi_p_max-f_phi_p_min)/f_N_phi_p+f_phi_p_min;

  Double_t w=f_rand.Uniform(w_min,w_max);
  Double_t qbold=f_rand.Uniform(qbold_min,qbold_max);
  Double_t mag_p=f_rand.Uniform(mag_p_min,mag_p_max);
  Double_t phi_p=f_rand.Uniform(phi_p_min,phi_p_max);

  event->Init(f_params,f_nucleus,f_Enu,w,qbold,mag_p,phi_p);
  event->Evaluate_dsigma_dall();

  //rejection method
  Double_t y=f_rand.Uniform(0.0,f_surface[i_w][i_qbold][i_mag_p][i_phi_p]);

  weight=1.0;
  if (f_surface[i_w][i_qbold][i_mag_p][i_phi_p]<event->f_dsigma_dall) {
    //cout << "broken surface=" << Index(i_w,i_qbold,i_mag_p,i_phi_p) << "," << f_surface[i_w][i_qbold][i_mag_p][i_phi_p]/event->f_dsigma_dall << endl;
    weight=event->f_dsigma_dall/f_surface[i_w][i_qbold][i_mag_p][i_phi_p];
  }

  if (y<event->f_dsigma_dall) return kTRUE;
  else return kFALSE;
}

////////////////////////////////////////////////////////////////////////
Bool_t TT_surface::Simple_draw_point(Double_t height,TT_event *event,Double_t& weight)
{
  Double_t p_lep=sqrt(pow(f_Enu-f_w_min,2)-pow(f_params->f_m_lep,2));
  Double_t qbold_min=f_Enu-p_lep;
  Double_t qbold_max=f_Enu+p_lep;

  Double_t w=f_rand.Uniform(f_w_min,f_w_max);
  Double_t qbold=f_rand.Uniform(qbold_min,qbold_max);
  Double_t mag_p=f_rand.Uniform(f_mag_p_min,f_mag_p_max);
  Double_t phi_p=f_rand.Uniform(f_phi_p_min,f_phi_p_max);

  if (f_params->f_do_AS_corr) {
    Double_t cos_theta_pq=f_rand.Uniform(-1.0,1.0);
    event->Init(f_params,f_nucleus,f_Enu,w,qbold,mag_p,cos_theta_pq,phi_p);
  }
  else event->Init(f_params,f_nucleus,f_Enu,w,qbold,mag_p,phi_p);

  event->Evaluate_dsigma_dall();

  //rejection method
  Double_t y=f_rand.Uniform(0.0,height);

  weight=1.0;
  if (height<event->f_dsigma_dall) {
    weight=event->f_dsigma_dall/height;
    printf("broken surface; r=%8.6g; h=%8.6g; xs=%8.6g; w=%8.6g; qbold=%8.6g; mag_p=%8.6g; phi_p=%8.6g\n",event->f_dsigma_dall/height,height,event->f_dsigma_dall,w,qbold,mag_p,phi_p);
  }

  if (y<event->f_dsigma_dall) return kTRUE;
  else return kFALSE;
}

