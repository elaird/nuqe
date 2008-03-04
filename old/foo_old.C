//implementing nucl-th/0512004v4
// ./m to compile

////list of ambiguities/issues
//--de Forest prescription
//--ksi in F2?
//--Mpi+-0?
//--poles in Smith-Moniz form factors?
//--M value fixed?
//--use delta function to eliminate the integral over E?
//--overall sign of epsilon
//--epsilon upper/lower -1 proof
//--what if contraction is complex?
//--sign of e_bind
//--understand equation (21)
//--normalize SF to #neutrons
//--args, options, seed
//--form factors
//--integrate over a flux?
//--nuclear magneton
//--hadronic tensor (M!=MA?)
//--proton mass!=neutron mass
//--generalize target and recoil masses (call M_prime) for e scattering?
//--average e_bind for SM
//--Double_t ok?
//--direction of q?
//--negative E/limits of integration for E?
//----has implications for norm. of SM spectral function
//----E lower threshold for binding e near 0
//--does having two p_mag solutions make sense?
//--virtual functions?
//--check cos_theta_pq uniformity
//--single nucleon lower energy thershold
//--overall factor of M/(M+w) between Katori and polish limit?
//--Llewellyn-Smith does not make sense at high Q2
//
//--for (pf,eb)=(0.001,0.0), nothing      >
//--for (pf,eb)=(0.03,0.0), looks good    >  understand
//--for (pf,eb)=(0.225,0.0), very little  >
//--for (pf,eb)=(0.225,0.025), seems ok   >
//--seems caused by a N_W that's too small

#include <iostream>
#include "Rtypes.h"
#include "TRandom3.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TFile.h"
#include "TTree.h"
#include "TT_param_class.h"
#include "TV_class.h"

using namespace std;

Int_t xs_e(TT_param_class*);
Int_t xs_w_q(TT_param_class*);
Int_t xs2_w_q(TT_param_class*);
Int_t xs_w_Q2(TT_param_class*);
Int_t xs_Q2(TT_param_class*);

////////////////////////////////////////////////////////////////////////
Int_t main(Int_t argc,char **argv) {
  Bool_t help=0;

  if (help) {
    cout << "--USAGE: " << argv[0] << " CONFIG_FILE" << endl;
    cout << "--for further help, see readme" << endl;
    return 0;
  }
  
  //read in config file
  TT_param_class *params = new TT_param_class(argv[1]);
  cout << params << " " << endl;

  //xs_w_Q2(params);
  cout << "1" << endl;
  Double_t Enu=0.8;

  const Double_t unit_conv=0.389379e-27;  

  const Int_t N_W=2000;
  const Int_t N_Q2=200;
  
  Double_t  w[N_W*N_Q2];
  Double_t Q2[N_W*N_Q2];
//  Double_t  E[N_W*N_Q2];
//  Double_t  C[N_W*N_Q2];
  Double_t XS2d[N_W*N_Q2];
  Double_t XS2d_err[N_W*N_Q2];
  cout << "2" << endl;
  Double_t Q21d[N_Q2];
  Double_t XS1d[N_Q2];
  Double_t XS1d_err[N_Q2];
  cout << "3" << endl;

  delete params;
  return 0;
}

////////////////////////////////////////////////////////////////////////
Int_t xs_w_q(TT_param_class *params) {

  Double_t Enu=0.8;

  const Double_t unit_conv=0.389379e-27;  

  Int_t N_E=200;
  Int_t N_C=200;
  
  Double_t      w[N_E*N_C];
  Double_t q_bold[N_E*N_C];
  Double_t  E[N_E*N_C];
  Double_t  C[N_E*N_C];
  Double_t XS2d[N_E*N_C];
  Double_t XS2d_err[N_E*N_C];

  Double_t  E1d[N_E];
  Double_t XS1d[N_E];
  Double_t XS1d_err[N_E];
  
  Double_t E_min=params->f_m_lep;
  Double_t E_max=Enu;
  Double_t C_min=-1.0;
  Double_t C_max=1.0;

  for (Int_t iE=0;iE<N_E;iE++) {
    XS1d[iE]=0.0;
    for (Int_t iC=0;iC<N_C;iC++) {
  
      Int_t index=iE*N_C+iC;
      E[index]=iE*(E_max-E_min)/N_E+E_min;
      C[index]=iC*(C_max-C_min)/N_C+C_min;

      Double_t P=sqrt(pow(E[index],2)-pow(params->f_m_lep,2));
      Double_t q_lower[4];
      q_lower[0]=Enu-E[index];
      q_lower[1]=0.0;
      q_lower[2]=-P*sqrt(1.0-pow(C[index],2));
      q_lower[3]=Enu-P*C[index];

      w[index]=q_lower[0];
      q_bold[index]=0.0;
      for (Int_t i=1;i<4;i++) q_bold[index]+=pow(q_lower[i],2);
      q_bold[index]=sqrt(q_bold[index]);

      XS2d    [index]=0.0;
      XS2d_err[index]=0.0;

      TV_class foo(params,Enu,q_lower);
      Double_t var=TMath::Abs(Enu-foo.SM_compute_enuqe())/Enu;
      if (var<0.05) {
	foo.Compute_XS();
	
	XS2d    [index]=unit_conv*foo.f_d2sigma_dw_dq;
	XS2d_err[index]=unit_conv*foo.f_d2sigma_dw_dq_err;

	Double_t jacobian=Enu*P/sqrt(Enu*Enu-2.0*Enu*P*C[index]+P*P);
	XS1d[iE] += jacobian*XS2d[index]*(C_max-C_min)/N_C;
      }
      //if (XS2d[index]>0.0) printf("enuqe=%5.3g; E=%5.3g; C=%5.3g; XS=%5.3g; err=%5.3g; \n",foo.SM_compute_enuqe(),E[index],C[index],XS2d[index],XS2d_err[index]);
      //printf("enuqe=%5.3g; E=%5.3g; C=%5.3g; XS=%5.3g; err=%5.3g; \n",foo.SM_compute_enuqe(),E[index],C[index],XS2d[index],XS2d_err[index]);
    }
    E1d[iE]=E[iE*N_C];
    cout << "E=" << E1d[iE] << endl;
  }

  TGraph xs1d_graph(N_E,E1d,XS1d);
  xs1d_graph.SetName("xs1d");
  xs1d_graph.SetTitle(";E_{#mu} (GeV);d#sigma/dE_{#mu} (cm^{2}/GeV)");
  
  TGraph2D xs2d_graph(N_E*N_C,E,C,XS2d);
  xs2d_graph.SetName("xs2d");
  xs2d_graph.SetTitle(";E_{#mu} (GeV);cos(#theta_{#mu#nu});d#sigma/dE_{#mu}dcos(#theta) (cm^{2}/GeV)");
  TGraph2D xs2d_err_graph(N_E*N_C,E,C,XS2d_err);
  xs2d_err_graph.SetName("xs2d_errs");
  xs2d_err_graph.SetTitle(";x-title;y-title;z-title");
  TFile f("~elaird/ccqe/xs_ted/xs.root","RECREATE");
  xs1d_graph.Write();
  xs2d_graph.Write();
  xs2d_err_graph.Write();
  f.Close();

  return 0;
}

////////////////////////////////////////////////////////////////////////
Int_t xs2_w_q(TT_param_class *params) {

  Double_t Enu=0.8;

  const Double_t unit_conv=0.389379e-27;  

  Int_t N_W=200;
  Int_t N_Q_BOLD=200;
  
  Double_t      w[N_W*N_Q_BOLD];
  Double_t q_bold[N_W*N_Q_BOLD];
  Double_t  E[N_W*N_Q_BOLD];
  Double_t  C[N_W*N_Q_BOLD];
  Double_t XS2d[N_W*N_Q_BOLD];
  Double_t XS2d_err[N_W*N_Q_BOLD];

  Double_t  w1d[N_W];
  Double_t XS1d[N_W];
  Double_t XS1d_err[N_W];
  
  Double_t w_min=0.0;
  Double_t w_max=Enu-params->f_m_lep;

  for (Int_t iW=0;iW<N_W;iW++) {
    XS1d[iW]=0.0;
    for (Int_t iQ_BOLD=0;iQ_BOLD<N_Q_BOLD;iQ_BOLD++) {
  
      Int_t index=iW*N_Q_BOLD+iQ_BOLD;
      XS2d    [index]=0.0;
      XS2d_err[index]=0.0;

      w[index]=iW*(w_max-w_min)/N_W+w_min;

      Double_t q_bold_min=0.0;
      Double_t q_bold_max=2.0*sqrt(pow(w[index],2)+2.0*params->f_M*w[index]);

      q_bold[index]=(iQ_BOLD+1)*(q_bold_max-q_bold_min)/N_Q_BOLD+q_bold_min;

      Double_t p_lep_sq=pow(Enu-w[index],2)-pow(params->f_m_lep,2);
      Double_t cos_theta_q=(pow(Enu,2)-p_lep_sq+pow(q_bold[index],2))/(2.0*Enu*q_bold[index]);
      if (cos_theta_q>1.0) continue;
      Double_t sin_theta_q=sqrt(1.0-pow(cos_theta_q,2));
      Double_t sin_theta_lep=-q_bold[index]*sin_theta_q/sqrt(p_lep_sq);
      //printf("w=%5.3g; q=%5.3g; c_q=%5.3g; s_l=%5.3g\n",w[index],q_bold[index],cos_theta_q,sin_theta_q);

      Double_t q_lower[4];
      q_lower[0]=w[index];
      q_lower[1]=0.0;
      q_lower[2]=q_bold[index]*sin_theta_q;
      q_lower[3]=q_bold[index]*cos_theta_q;


      TV_class foo(params,Enu,q_lower);
      {
	foo.Compute_XS();
	
	XS2d    [index]=unit_conv*foo.f_d2sigma_dw_dq;
	XS2d_err[index]=unit_conv*foo.f_d2sigma_dw_dq_err;

	XS1d[iW] += XS2d[index]*(q_bold_max-q_bold_min)/N_Q_BOLD;
      }
//      //if (XS2d[index]>0.0) printf("enuqe=%5.3g; E=%5.3g; C=%5.3g; XS=%5.3g; err=%5.3g; \n",foo.SM_compute_enuqe(),w[index],q_bold[index],XS2d[index],XS2d_err[index]);
//      //printf("enuqe=%5.3g; E=%5.3g; C=%5.3g; XS=%5.3g; err=%5.3g; \n",foo.SM_compute_enuqe(),w[index],q_bold[index],XS2d[index],XS2d_err[index]);
    }
    w1d[iW]=w[iW*N_Q_BOLD];
    cout << "w=" << w1d[iW] << endl;
  }

  TGraph xs1d_graph(N_W,w1d,XS1d);
  xs1d_graph.SetName("xs1d");
  xs1d_graph.SetTitle(";#omega (GeV);d#sigma/d#omega} (cm^{2}/GeV)");
  
  TGraph2D xs2d_graph(N_W*N_Q_BOLD,E,C,XS2d);
  xs2d_graph.SetName("xs2d");
  xs2d_graph.SetTitle(";E_{#mu} (GeV);cos(#theta_{#mu#nu});d#sigma/dE_{#mu}dcos(#theta) (cm^{2}/GeV)");
  TGraph2D xs2d_err_graph(N_W*N_Q_BOLD,E,C,XS2d_err);
  xs2d_err_graph.SetName("xs2d_errs");
  xs2d_err_graph.SetTitle(";x-title;y-title;z-title");
  TFile f("~elaird/ccqe/xs_ted/xs.root","RECREATE");
  xs1d_graph.Write();
  xs2d_graph.Write();
  xs2d_err_graph.Write();
  f.Close();

  return 0;
}

////////////////////////////////////////////////////////////////////////
Int_t xs_w_Q2(TT_param_class *params) {

  Double_t Enu=0.8;

  const Double_t unit_conv=0.389379e-27;  

  const Int_t N_W=2000;
  const Int_t N_Q2=200;
  
  Double_t Q21d[N_Q2];
  Double_t XS1d[N_Q2];
  Double_t XS1d_err[N_Q2];

  Double_t w_min=0.0;
  Double_t w_max=Enu-params->f_m_lep;

  Double_t Q2_min=0.0;
  Double_t Q2_max=2.0*params->f_M*w_max;
  //Double_t Q2_max=0.2;

  for (Int_t iQ2=0;iQ2<N_Q2;iQ2++) {
    XS1d[iQ2]=0.0;
    for (Int_t iW=0;iW<N_W;iW++) {
    
      Int_t index=iQ2*N_W+iW;
    
      Double_t w=iW*(w_max-w_min)/N_W+w_min;
      Double_t Q2=iQ2*(Q2_max-Q2_min)/N_Q2+Q2_min;
    
      Double_t q_bold_mag=sqrt(Q2+w*w);
      Double_t p_lep_sq=pow(Enu-w,2)-pow(params->f_m_lep,2);
      Double_t cos_theta_q=(pow(Enu,2)-p_lep_sq+pow(q_bold_mag,2))/(2.0*Enu*q_bold_mag);
      if (cos_theta_q>1.0) continue;
      Double_t sin_theta_q=sqrt(1.0-pow(cos_theta_q,2));
      //Double_t sin_theta_lep=-Q2*sin_theta_q/sqrt(p_lep_sq);
      //printf("w=%5.3g; q=%5.3g; c_q=%5.3g; s_l=%5.3g\n",w,Q2,cos_theta_q,sin_theta_q);
    
      Double_t q_lower[4];
      q_lower[0]=w;
      q_lower[1]=0.0;
      q_lower[2]=q_bold_mag*sin_theta_q;
      q_lower[3]=q_bold_mag*cos_theta_q;
    
      TV_class foo(params,Enu,q_lower);
      foo.Compute_XS();
      
      Double_t jacobian=0.5/q_bold_mag;
      XS2d    =unit_conv*jacobian*foo.f_d2sigma_dw_dq;
      XS2d_err=unit_conv*jacobian*foo.f_d2sigma_dw_dq_err;
      XS1d[iQ2] += XS2d*(w_max-w_min)/N_W;
      
      //if (XS2d>0.0) printf("enuqe=%5.3g; E=%5.3g; C=%5.3g; XS=%5.3g; err=%5.3g; \n",foo.SM_compute_enuqe(),w,Q2,XS2d,XS2d_err);
      //printf("enuqe=%5.3g; E=%5.3g; C=%5.3g; XS=%5.3g; err=%5.3g; \n",foo.SM_compute_enuqe(),w,Q2,XS2d,XS2d_err);
    }
    Q21d[iQ2]=iQ2*(Q2_max-Q2_min)/N_Q2+Q2_min;
    printf("Q2=%5.3g;xs=%5.3g\n",Q21d[iQ2],XS1d[iQ2]);
  }

  TGraph xs1d_graph(N_Q2,Q21d,XS1d);
  xs1d_graph.SetName("xs1d");
  xs1d_graph.SetTitle(";Q^{2} (GeV^{2});d#sigma/dQ^{2} (cm^{2}/GeV^{2})");
  
  //TGraph2D xs2d_graph(N_W*N_Q2,E,C,XS2d);
  //xs2d_graph.SetName("xs2d");
  //xs2d_graph.SetTitle(";E_{#mu} (GeV);cos(#theta_{#mu#nu});d#sigma/dE_{#mu}dcos(#theta) (cm^{2}/GeV)");
  //TGraph2D xs2d_err_graph(N_W*N_Q2,E,C,XS2d_err);
  //xs2d_err_graph.SetName("xs2d_errs");
  //xs2d_err_graph.SetTitle(";x-title;y-title;z-title");
  TFile f("~elaird/ccqe/xs_ted/xs.root","RECREATE");
  xs1d_graph.Write();
  //xs2d_graph.Write();
  //xs2d_err_graph.Write();
  f.Close();

  delete [] w       ;
  delete [] Q2      ;
  delete [] E       ;
  delete [] C       ;
  delete [] XS2d    ;
  delete [] XS2d_err;

  return 0;
}

////////////////////////////////////////////////////////////////////////
Int_t xs_Q2(TT_param_class *params) {

  Double_t Enu=0.8;

  const Double_t unit_conv=0.389379e-27;  

  const Int_t N_Q2=200;
  
  Double_t Q2      [N_Q2];
  Double_t XS1d    [N_Q2];
  Double_t XS1d_err[N_Q2];
  
  Double_t Q2_min=pow(params->f_M,2)-pow(params->f_M_p,2);
  Double_t Q2_max=Q2_min+2.0*params->f_M*(Enu-params->f_m_lep);

  for (Int_t iQ2=0;iQ2<N_Q2;iQ2++) {
  
    Q2      [iQ2]=iQ2*(Q2_max-Q2_min)/N_Q2+Q2_min;
    XS1d    [iQ2]=0.0;
    XS1d_err[iQ2]=0.0;

    Double_t w=(Q2[iQ2]-pow(params->f_M,2)+pow(params->f_M_p,2))/(2.0*params->f_M);
    Double_t q_bold_mag=sqrt(Q2[iQ2]+w*w);

    Double_t p_lep_sq=pow(Enu-w,2)-pow(params->f_m_lep,2);
    Double_t cos_theta_q=(pow(Enu,2)-p_lep_sq+pow(q_bold_mag,2))/(2.0*Enu*q_bold_mag);
    if (cos_theta_q>1.0) continue;
    printf("w=%5.3g; Q2=%5.3g; q_bold_mag=%5.3g; c_q=%5.3g\n",w,Q2[iQ2],q_bold_mag,cos_theta_q);
    Double_t sin_theta_q=sqrt(1.0-pow(cos_theta_q,2));
    //Double_t sin_theta_lep=-Q2[iQ2]*sin_theta_q/sqrt(p_lep_sq);

    Double_t q_lower[4];
    q_lower[0]=w;
    q_lower[1]=0.0;
    q_lower[2]=q_bold_mag*sin_theta_q;
    q_lower[3]=q_bold_mag*cos_theta_q;

    params->f_single_nucleon_mode=1;
    params->f_use_M_instead_of_MplusW=0;
    TV_class foo(params,Enu,q_lower);
    foo.Compute_XS();

    Double_t jacobian=0.5/q_bold_mag;
    XS1d    [iQ2]=unit_conv*jacobian*foo.f_d2sigma_dw_dq;
    XS1d_err[iQ2]=unit_conv*jacobian*foo.f_d2sigma_dw_dq_err;
    
    //      //if (XS2d[iQ2]>0.0) printf("enuqe=%5.3g; E=%5.3g; C=%5.3g; XS=%5.3g; err=%5.3g; \n",foo.SM_compute_enuqe(),w[iQ2],Q2[iQ2],XS2d[iQ2],XS2d_err[iQ2]);
    //      //printf("enuqe=%5.3g; E=%5.3g; C=%5.3g; XS=%5.3g; err=%5.3g; \n",foo.SM_compute_enuqe(),w[iQ2],Q2[iQ2],XS2d[iQ2],XS2d_err[iQ2]);
  }

  TGraph xs1d_graph(N_Q2,Q2,XS1d);
  xs1d_graph.SetName("xs1d");
  xs1d_graph.SetTitle(";Q^{2} (GeV^{2});d#sigma/dQ^{2} (cm^{2}/GeV^{2})");
  
  TFile f("~elaird/ccqe/xs_ted/xs.root","RECREATE");
  xs1d_graph.Write();
  f.Close();

  return 0;
}

////////////////////////////////////////////////////////////////////////
//Int_t xs_e(TT_param_class *params) {
//
//  //tree for debuging
//  TTree *tree = new TTree("tree","");
//  Int_t tree_problem;
//  Double_t tree_E_lepprime;
//  Double_t tree_p_lower[4];
//  tree->Branch("problem_branch",&tree_problem,"problem/I");
//  tree->Branch("E_lepprime_branch",&tree_E_lepprime,"E_lepprime/D");
//  tree->Branch("p_lower_branch",tree_p_lower,"p_lower[4]/D");
//
//  Double_t Enu=0.8;
//  const Double_t unit_conv=0.389379e-30;
//  Int_t N_E=1;
//  
//  Double_t E_lepprime[N_E];
//  Double_t         XS[N_E];
//  Double_t     XS_err[N_E];
//  Double_t    XS_xerr[N_E];
//  
//  Double_t E_min=Enu;//1.01*params->f_m_lep;
//  Double_t E_max=Enu;///2.0;
//
//  for (Int_t iE=0;iE<N_E;iE++) {
//    E_lepprime[iE]=iE*(E_max-E_min)/N_E+E_min;
//  }
//
//  for (Int_t iE=0;iE<N_E;iE++) {
//    TU_class foo(params,Enu,E_lepprime[iE],tree,&tree_E_lepprime,tree_p_lower,&tree_problem);
//    foo.Compute_XS();
//
//    XS    [iE]=unit_conv*foo.f_d3sigma_d3kprime;
//    XS_err[iE]=unit_conv*foo.f_d3sigma_d3kprime_err;
//    XS_xerr[iE]=0.0;
//    cout << "forbidden, problems, E, xs, err: " << foo.f_d3sigma_d3kprime_forbidden << ","
//	 << foo.f_d3sigma_d3kprime_problems << "," << E_lepprime[iE] << "," << XS[iE] << "," << XS_err[iE] << endl;
//  }
//  
//  TGraphErrors graph(N_E,E_lepprime,XS,XS_xerr,XS_err);
//  graph.SetName("xs_e");
//  graph.SetTitle(";E_{#mu} (GeV);d#sigma/dE_{#mu} (cm^{2}/MeV);z-title");
//  
//  TFile f("~elaird/ccqe/xs_ted/xs.root","RECREATE");
//  graph.Write();
//  tree->Write();
//  f.Close();
//
//  return 0;
//}

