#include <iostream>
#include "Rtypes.h"

using namespace std;

////////////////////////////////////////////////////////////////////////
Int_t main(Int_t argc,char **argv) {
  Bool_t help=0;

  Long64_t l6_1=1;
  Long64_t l6_2=1;

  for (Int_t i=1;i<63;i++) {
    l6_1*=2;
    l6_2*=2;
    Double_t n1=l6_1;
    Double_t n2 =1.0;
    Double_t d=l6_2;
    Double_t ratio=(n1+n2)/d;
    cout << "i=" << i << "; l6_1=" << l6_1 << "; ratio=" << ratio << endl;

    Double_t s1=1.0/l6_1;
    Double_t s2 =1.0;
    Double_t t=s1+s2;
    cout << "i=" << i << "; l6_1=" << l6_1 << "; sum=" << t << endl;
  }

  return 0;
}
