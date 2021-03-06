#include "Rtypes.h"

Int_t power(Int_t base,Int_t exp) {
  Int_t ans=1;
  for (Int_t i=0;i<exp;i++) ans*=base;
  return ans;
}

//bubble sort!
Int_t count_transpositions(const Int_t a[]) {
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

void upper(Double_t g[4][4],Double_t a_lower[4],Double_t a_upper[4]){

  for (Int_t i=0;i<4;i++) a_upper[i]=0.0;
  for (Int_t mu=0;mu<4;mu++) {
    for (Int_t nu=0;nu<4;nu++) {
      a_upper[nu]+=a_lower[mu]*g[mu][nu];
    }
  }
}
