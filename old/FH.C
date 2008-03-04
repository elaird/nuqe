#include <iostream>
#include <math.h>
#include "Rtypes.h"

//form factors are from hep-ex/0308005 and Phys. Rev. C 51 409 (1995)
//(references from hep-ex/0603034)
//(BBA 2003)


const Double_t gA=-1.267;
const Double_t axial_mass=1.0;

Double_t G_pE(Double_t q2) {
  return 1.0;
}

Double_t G_nE(Double_t q2) {
  return 1.0;
}

Double_t G_pM(Double_t q2) {
  return 1.0;
}

Double_t G_nM(Double_t q2) {
  return 1.0;
}



Double_t G_VE(Double_t q2) {
  return G_pE(q2)-G_nE(q2);
}

Double_t G_VM(Double_t q2) {
  return G_pM(q2)-G_nM(q2);
}



Double_t F1(Double_t q2,Double_t M) {
  return (G_VE(q2)-G_VM(q2)*q2/(4.0*M*M)) / (1.0-q2/(4.0*M*M));
}

Double_t F2(Double_t q2,Double_t M) {
  return (G_VM(q2)-G_VE(q2)) / (1.0-q2/(4.0*M*M));
}

Double_t FA(Double_t q2) {
  return gA/pow(1.0-q2/axial_mass/axial_mass,2);
}

Double_t FP(Double_t q2,Double_t M) {
  const Double_t Mpi=0.13957;
  return 2*M*M*FA(q2)/(Mpi*Mpi-q2);
}



Double_t H1(Double_t q2,Double_t M) {
  Double_t tau=-q2/(4.0*M*M);
  return pow(FA(q2),2)*(1+tau) + tau*pow(F1(q2,M)+F2(q2,M),2);
}

Double_t H2(Double_t q2,Double_t M) {
  Double_t tau=-q2/(4.0*M*M);
  return pow(FA(q2),2) + pow(F1(q2,M),2) + tau*pow(F2(q2,M),2);
}

Double_t H3(Double_t q2,Double_t M) {
  return 2*FA(q2)*(F1(q2,M)+F2(q2,M));
}

Double_t H4(Double_t q2,Double_t M) {
  Double_t tau=-q2/(4.0*M*M);
  return pow(F2(q2,M),2)*(1-tau)/4.0 + F1(q2,M)*F2(q2,M)/2.0 + FA(q2)*FP(q2,M) - tau*pow(FP(q2,M),2);
}

Double_t H5(Double_t q2,Double_t M) {
  return H2(q2,M);
}

