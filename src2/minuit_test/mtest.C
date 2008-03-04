#include "TMinuit.h"

////////////////////////////////////////////////////////////////////////
int main()
{

}

//////////////////////////////////////////////////////////////////////////
//void Minuit_FCN(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
//{
//  //TV_class foo;
//  //foo.Init(f_params,f_nucleus,f_Enu,par[0],par[1],par[2],par[3]);
//  //foo.Evaluate_d4sigma_dw_dq_bold_domega_p();
//  //f=foo.f_d4sigma_dw_dq_bold_domega_p;
//}

////////////////////////////////////////////////////////////////////////
Double_t f1(Int_t factor)
{
  //static?
  Double_t vstart[4];
  Double_t step[4];

  //vstart[0]=f_rand.Uniform(w_min,w_max);
  //vstart[1]=f_rand.Uniform(qbold_min,qbold_max);
  //vstart[2]=f_rand.Uniform(cth_p_min,cth_p_max);
  //vstart[3]=f_rand.Uniform(phi_p_min,phi_p_max);
	    
  TMinuit minuit(4);

  //minuit.SetFCN(Minuit_FCN);

  //Double_t arglist[10];
  //Int_t ierflg = 0;
  //
  //arglist[0] = 1;
  //minuit.mnexcm("SET ERR", arglist ,1,ierflg);
  //
  //minuit.mnparm(0, "w",      vstart[0], step[0], 0,0,ierflg);
  //minuit.mnparm(1, "q_bold", vstart[1], step[1], 0,0,ierflg);
  //minuit.mnparm(2, "cth_p",  vstart[2], step[2], 0,0,ierflg);
  //minuit.mnparm(3, "phi_p",  vstart[3], step[3], 0,0,ierflg);
  //
  //// Now ready for minimization step
  //arglist[0] = 500;
  //arglist[1] = 1.;
  //minuit.mnexcm("MIGRAD", arglist ,2,ierflg);
  //
  //// Print results
  //Double_t amin,edm,errdef;
  //Int_t nvpar,nparx,icstat;
  //minuit.mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  //minuit.mnprin(3,amin);


  return 1.0;
}

