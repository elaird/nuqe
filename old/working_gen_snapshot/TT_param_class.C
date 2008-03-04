#define TT_param_class_cxx
#include "TT_param_class.h"

////////////////////////////////////////////////////////////////////////
void TT_param_class::Init_g()
{
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
void TT_param_class::Init_epsilon_lower()
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
Int_t TT_param_class::Power(Int_t base,Int_t exp) {
  Int_t ans=1;
  for (Int_t i=0;i<exp;i++) ans*=base;
  return ans;
}

////////////////////////////////////////////////////////////////////////
Int_t TT_param_class::Count_transpositions(const Int_t a[]) {
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
void TT_param_class::Upper(Double_t g[4][4],Double_t a_lower[4],Double_t a_upper[4]){
  //g is diagonal
  for (Int_t mu=0;mu<4;mu++) a_upper[mu]=a_lower[mu]*g[mu][mu];
}

////////////////////////////////////////////////////////////////////////
Double_t TT_param_class::G_D(Double_t q2)
{
  return 1.0/pow(1.0-q2/(f_MV*f_MV),2);
}

//////////////////////////////////////////////////////////////////////////
//Double_t TT_param_class::G_pE(Double_t q2)
//{
//  return G_D(q2);
//}
//
//////////////////////////////////////////////////////////////////////////
//Double_t TT_param_class::G_nE(Double_t q2)
//{
//  return 0.0;
//}
//
//////////////////////////////////////////////////////////////////////////
//Double_t TT_param_class::G_pM(Double_t q2)
//{
//  return f_mu_p_factor*f_mu_N*G_D(q2);
//}
//
//////////////////////////////////////////////////////////////////////////
//Double_t TT_param_class::G_nM(Double_t q2)
//{
//  return f_mu_n_factor*f_mu_N*G_D(q2);
//}
//
//////////////////////////////////////////////////////////////////////////
//Double_t TT_param_class::G_VE(Double_t q2)
//{
//  if (f_do_Paschos_FF || f_do_Benhar_FF) return 1.0/pow(1.0-q2/(f_MV*f_MV),2);
//  return G_pE(q2)-G_nE(q2);
//}
//
//////////////////////////////////////////////////////////////////////////
//Double_t TT_param_class::G_VM(Double_t q2)
//{
//  if (f_do_Paschos_FF || f_do_Benhar_FF) return 4.706*G_VE(q2);
//  return G_pM(q2)-G_nM(q2);
//}
//
//////////////////////////////////////////////////////////////////////////
//Double_t TT_param_class::F1_old(Double_t q2)
//{
//  return (G_VE(q2)-G_VM(q2)*q2/(4.0*f_M*f_M)) / (1.0-q2/(4.0*f_M*f_M));
//}
//
//////////////////////////////////////////////////////////////////////////
//Double_t TT_param_class::F2_old(Double_t q2)
//{
//  if (f_do_Benhar_FF)  return (G_VM(q2)-G_VE(q2)) / (1.0-q2/(4.0*f_M*f_M));
//  if (f_do_Paschos_FF) return (G_VM(q2)-G_VE(q2)) / (1.0-q2/(4.0*f_M*f_M));
//  //or perhaps this?
//  //if (f_do_Paschos_FF) return (G_VM(q2)-G_VE(q2)) / (1.0-q2/(4.0*f_M*f_M)) / f_xi;
//  return (G_VM(q2)-G_VE(q2)) / (1.0-q2/(4.0*f_M*f_M));
//}

////////////////////////////////////////////////////////////////////////
Double_t TT_param_class::F1(Double_t q2)
{
  return (1.0-q2*(1.0+f_xi)/(4.0*f_M*f_M))*G_D(q2) / (1.0-q2/(4.0*f_M*f_M));
}

////////////////////////////////////////////////////////////////////////
Double_t TT_param_class::F2(Double_t q2)
{
  return f_xi*G_D(q2) / (1.0-q2/(4.0*f_M*f_M));
}

////////////////////////////////////////////////////////////////////////
Double_t TT_param_class::FA(Double_t q2)
{
  return f_gA/pow(1.0-q2/(f_axial_mass*f_axial_mass),2);
}

////////////////////////////////////////////////////////////////////////
Double_t TT_param_class::FP(Double_t q2)
{
  return 2*f_M*f_M*FA(q2)/(f_Mpi*f_Mpi-q2);
}

////////////////////////////////////////////////////////////////////////
Double_t TT_param_class::H1(Double_t q2)
{
  Double_t tau=-q2/(4.0*f_M*f_M);
  return pow(FA(q2),2)*(1+tau) + tau*pow(F1(q2)+F2(q2),2);
}

////////////////////////////////////////////////////////////////////////
Double_t TT_param_class::H2(Double_t q2)
{
  Double_t tau=-q2/(4.0*f_M*f_M);
  return pow(FA(q2),2) + pow(F1(q2),2) + tau*pow(F2(q2),2);
}

////////////////////////////////////////////////////////////////////////
Double_t TT_param_class::H3(Double_t q2)
{
  return 2.0*FA(q2)*(F1(q2)+F2(q2));
}

////////////////////////////////////////////////////////////////////////
Double_t TT_param_class::H4(Double_t q2)
{
  Double_t tau=-q2/(4.0*f_M*f_M);
  return pow(F2(q2),2)*(1-tau)/4.0 + F1(q2)*F2(q2)/2.0 + FA(q2)*FP(q2) - tau*pow(FP(q2),2);
}

////////////////////////////////////////////////////////////////////////
Double_t TT_param_class::H5(Double_t q2)
{
  return H2(q2);
}

