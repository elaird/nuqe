#include <iostream>
#include "Rtypes.h"

void g_init(Double_t a[][4]) {
  for (Int_t mu=0;mu<4;mu++) {
    for (Int_t nu=0;nu<4;nu++) {
      a[mu][nu]=0.0;
    }
  }
  a[0][0]= 1.0;
  a[1][1]=-1.0;
  a[2][2]=-1.0;
  a[3][3]=-1.0;
}

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

void epsilon_lower_init(Double_t e[][4][4][4]) {
  Int_t i[4];
  for (i[0]=0;i[0]<4;i[0]++) {
    for (i[1]=0;i[1]<4;i[1]++) {
      for (i[2]=0;i[2]<4;i[2]++) {
	for (i[3]=0;i[3]<4;i[3]++) {
	  Int_t c=0;
	  for (Int_t j=0;j<4;j++) c+=1<<i[j];

	  if (c==15) {
	    e[i[0]] [i[1]] [i[2]] [i[3]]=power(-1,count_transpositions(i));
	    //printf("%d,%d,%d,%d: %d,%g",i[0],i[1],i[2],i[3],count_transpositions(i),e[i[0]] [i[1]] [i[2]] [i[3]]);
	  }
	  else e[i[0]] [i[1]] [i[2]] [i[3]]=0.0;
	}
      }
    }
  }
}

//void epsilon_lower_init0(Double_t e[]) {
//  for (Int_t i=0;i<power(4,4);i++) {
//    e[i]=0;
//
//    Int_t index[4];
//    Int_t stuff=0;
//    for (Int_t j=0;j<4;j++) {
//      index[j]=(i%power(4,j+1)-stuff)/power(4,j);
//      stuff=index[j];
//    }
//    cout << Form("%d,%d,%d,%d",index[0],index[1],index[2],index[3]) << endl;
//    Int_t c=0;
//  }
//}

Double_t dot(Double_t g[4][4],Double_t a[4],Double_t b[4]){
  Double_t ans=0.0;
  for (Int_t mu=0;mu<4;mu++) {
    for (Int_t nu=0;nu<4;nu++) {
      ans+=a[mu]*g[mu][nu]*b[nu];
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
