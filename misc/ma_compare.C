#include "TFile.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TTree.h"
#include "TChain.h"
#include <iostream>
#include "TROOT.h"
#include "TEnv.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TF1.h"
#include "TLine.h"

////////////////////////////////////////////////////////////////////////
void Init() {
  gROOT->SetStyle("Plain");
  gEnv->SetValue("Root.Stacktrace","No");
}

////////////////////////////////////////////////////////////////////////
void XS_look() {
  Init();

  const Int_t NPoints=3;
  Int_t EnergyPoints[NPoints];

  EnergyPoints[0]=800;
  EnergyPoints[1]=1200;
  EnergyPoints[2]=0;

  TGraph g1(NPoints);
  TGraph g2(NPoints);
  TGraph g3(NPoints);
  TGraph g4(NPoints);

  for (Int_t iPoint=0;iPoint<NPoints;iPoint++) {
    Int_t energy_point=EnergyPoints[iPoint];

    //NUANCE mode
    TFile f4(Form("../../events/events_4_%d.root",energy_point));
    gROOT->cd();
    TGraph *gr4=(TGraph*)f4.Get("ccqe_rate")->Clone();
    Double_t rate4,dummy;
    gr4->GetPoint(4,dummy,rate4);
    delete gr4;

    //A-S SF mode
    TFile f23(Form("../../events/events_23_%d.root",energy_point));
    gROOT->cd();
    TGraph *gr23=(TGraph*)f23.Get("ccqe_rate")->Clone();
    Double_t rate2=0.0;
    Double_t rate3=0.0;
    gr23->GetPoint(2,dummy,rate2);
    gr23->GetPoint(3,dummy,rate3);
    delete gr23;
  
    //A-S SM mode
    TFile f1(Form("../../events/events_1_%d.root",energy_point));
    gROOT->cd();
    TGraph *gr1=(TGraph*)f1.Get("ccqe_rate")->Clone();
    Double_t rate1=0.0;
    gr1->GetPoint(1,dummy,rate1);
    delete gr1;

    g1.SetPoint(iPoint,energy_point,rate1);
    g2.SetPoint(iPoint,energy_point,rate2);
    g3.SetPoint(iPoint,energy_point,rate3);
    g4.SetPoint(iPoint,energy_point,rate4);
  }

  TCanvas can;

  g1.SetMarkerStyle(20);g1.SetMarkerColor(kBlack);
  g2.SetMarkerStyle(20);g2.SetMarkerColor(kRed);
  g3.SetMarkerStyle(20);g3.SetMarkerColor(kBlue);
  g4.SetMarkerStyle(20);g4.SetMarkerColor(kGreen);

  g1.Draw("ap");
  g2.Draw("psame");
  g3.Draw("psame");
  g4.Draw("psame");
  can.DrawClone();
}

////////////////////////////////////////////////////////////////////////
void MA_compare(Int_t mode,Int_t energy_point) {
  Init();

  //NUANCE mode
  TFile f4(Form("../../events/events_4_%d.root",energy_point));
  gROOT->cd();
  TTree *tree4=(TTree*)f4.Get("tree");
  tree4->SetName("tree4");
  TGraph *gr4=(TGraph*)f4.Get("ccqe_rate")->Clone();
  Double_t rate4,dummy;
  gr4->GetPoint(4,dummy,rate4);
  delete gr4;

  //A-S SF mode
  TFile f23(Form("../../events/events_23_%d.root",energy_point));
  gROOT->cd();
  TTree *tree23=(TTree*)f23.Get("tree");
  tree23->SetName("tree23");
  TGraph *gr23=(TGraph*)f23.Get("ccqe_rate")->Clone();
  Double_t rate23=0.0;
  Double_t dummy_x,dummy_y;
  gr23->GetPoint(2,dummy_x,dummy_y);
  rate23+=dummy_y;
  gr23->GetPoint(3,dummy_x,dummy_y);
  rate23+=dummy_y;
  delete gr23;
  
  //A-S SM mode
  TFile f1(Form("../../events/events_1_%d.root",energy_point));
  gROOT->cd();
  TTree *tree1=(TTree*)f1.Get("tree");
  tree1->SetName("tree1");
  TGraph *gr1=(TGraph*)f1.Get("ccqe_rate")->Clone();
  Double_t rate1=0.0;
  gr1->GetPoint(1,dummy_x,rate1);
  delete gr1;
  
  TString *tree_var=0;
  TString *label=0;

  Int_t N_bins=0;
  Double_t lower=0.0;
  Double_t upper=0.0;
  Double_t x_leg=0.0;

  switch (mode) {
  case 0:
    tree_var=new TString("k[0]");
    if (energy_point==0) {
      label=new TString(";E_{#nu} (GeV);event rate (arb. units)");
      x_leg=0.55;
      N_bins=70;
      lower=0.12;
      upper=3.0;
    }
    else {
      label=new TString(";E_{#nu} (GeV);d#sigma/dE_{#nu} (cm^{2}/GeV)");
      x_leg=0.55;
      N_bins=70;
      lower=0.12;
      upper=3.0;
    }
    break;
  case 1:
    tree_var=new TString("kprime[0]");
    if (energy_point==0) {
      label=new TString(";E_{#mu} (GeV);event rate (arb. units)");
      x_leg=0.55;
      N_bins=70;
      lower=0.12;
      upper=3.0;
    }
    else {
      label=new TString(";E_{#mu} (GeV);d#sigma/dE_{#mu} (cm^{2}/GeV)");
      x_leg=0.15;
      N_bins=70;
      lower=0.12;
      upper=energy_point*1.0e-3;
    }
    
    break;
  case 2:
    tree_var=new TString("-q[0]**2+q[1]**2+q[2]**2+q[3]**2");
    if (energy_point==0) label=new TString(";Q^{2} (GeV^{2});event rate (arb. units)");
    else                 label=new TString(";Q^{2} (GeV^{2});d#sigma/dQ^{2} (cm^{2}/GeV^{2})");

    x_leg=0.55;
    N_bins=70;
    lower=0.0;
    upper=1.4;
    break;
  case 3:
    tree_var=new TString("Q2qe");
    if (energy_point==0) label=new TString(";Q^{2}_{QE} (GeV^{2});event rate (arb. units)");
    else                 label=new TString(";Q^{2}_{QE} (GeV^{2});d#sigma/dQ^{2} (cm^{2}/GeV^{2})");

    x_leg=0.55;
    N_bins=70;
    lower=0.0;
    upper=1.4;
    break;
  case 4:
    tree_var=new TString("sqrt(pprime[1]**2+pprime[2]**2+pprime[3]**2)");
    label=new TString(";p_{f} (GeV/c);d#sigma/dp_{f} (cm^{2}/(GeV/c))");

    x_leg=0.15;
    N_bins=150;
    lower=0.0;
    upper=0.8;
    break;
  case 5:
    tree_var=new TString("sqrt(pprime[1]**2+pprime[2]**2+pprime[3]**2)");
    label=new TString(";p_{f} (GeV/c);d#sigma/dp_{f} (cm^{2}/(GeV/c))");

    x_leg=0.15;
    N_bins=150;
    lower=0.2;
    upper=0.6;
    //lower=0.2;
    //upper=1.3;
    break;
  case 6:
    tree_var=new TString("kprime[3]/kprime[0]");
    if (energy_point==0) label=new TString(";cos(#theta_{#mu});event rate (arb. units)");
    else                 label=new TString(";E_{#mu} (GeV);d#sigma/dE_{#mu} (cm^{2}/GeV)");
    
    x_leg=0.15;
    N_bins=70;
    lower=-1.0;
    upper=1.0;
    break;
  }

  Int_t NHistos=1;
  Int_t MA_offset=0;
  Double_t MA_14=1.03;
  Double_t MA_base=1.03;

  if (energy_point==0) {
    NHistos=4;
    MA_offset=1;
    MA_14=1.23;
    MA_base=1.03;
  }

  TH1D *h[NHistos];
  TH1D *r[NHistos];
  for (Int_t iHisto=0;iHisto<NHistos;iHisto++) {
    h[iHisto]=new TH1D(Form("h%d",iHisto),*label,N_bins,lower,upper);
    h[iHisto]->Sumw2();
  }
  TH1D p4("p4",*label,N_bins,lower,upper);
  p4.Sumw2();
  p4.SetStats(kFALSE);
  tree4->Draw(*tree_var+">>p4","process==4","goff");
  p4.Scale(rate4/p4.Integral("width"));

  TH1D p1("p1",*label,N_bins,lower,upper);
  p1.Sumw2();
  p1.SetStats(kFALSE);
  tree1->Draw(*tree_var+">>p1","process==1","goff");
  p1.Scale(rate1/p1.Integral("width"));
  p1.SetLineColor(kRed);
  p1.SetMarkerColor(kRed);

  Int_t color=3;
  for (Int_t iHisto=0;iHisto<NHistos;iHisto++) {
    tree23->Draw(*tree_var+Form(">>h%d",iHisto),Form("(MA_weights[%d]) * (process==2 || process==3)",iHisto),"goff");
    h[iHisto]->SetStats(kFALSE);
    h[iHisto]->GetXaxis()->CenterTitle();
    h[iHisto]->GetYaxis()->CenterTitle();
    h[iHisto]->GetYaxis()->SetTitleOffset(1.2);
    h[iHisto]->SetLineColor(color);
    h[iHisto]->SetMarkerColor(color);
    color++;
    h[iHisto]->Print();
  }

  TLegend legend(x_leg,0.6,x_leg+0.35,0.9);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetHeader(Form("dipole form factors, O^{16}, %d MeV #nu_{#mu}",energy_point));
  legend.AddEntry(&p4,Form("NUANCE mode, MA=%g GeV",MA_14));
  legend.AddEntry(&p1,Form("A.-S. SM mode, MA=%g GeV",MA_14));

  Double_t MA_values[NHistos];
  for (Int_t iMA=0;iMA<NHistos;iMA++) {
    MA_values[iMA]=MA_base+(iMA-MA_offset)*0.1;
  }

  TCanvas canvas("canvas","",700,1000);
  canvas.Divide(1,2);
  canvas.cd(1);
  
  p4.Draw("e");
  p1.Draw("esame");

  Double_t factor=rate23/h[0]->Integral("width");
  for (Int_t iHisto=0;iHisto<NHistos;iHisto++) {
    h[iHisto]->Scale(factor);
    h[iHisto]->Draw("esame"); 
    r[iHisto]=(TH1D*)h[iHisto]->Clone(Form("r%d",iHisto));
    legend.AddEntry(h[iHisto],Form("A.-S. SF mode, MA=%#4.3g GeV",MA_values[iHisto]));
  }
  legend.Draw("same");
  
  canvas.cd(2);
  TH2D null("null",Form(";;ratio to (NUANCE mode at %g GeV)",MA_14),1,lower,upper,1,0.6,1.2);
  null.SetStats(kFALSE);
  null.Draw();

  for (Int_t iHisto=0;iHisto<NHistos;iHisto++) {
    //r[iHisto]->Scale(1.0/r[iHisto]->Integral());
    r[iHisto]->Divide(&p4);
    r[iHisto]->Draw("same");
  }
  canvas.cd(0);
  canvas.DrawClone();

  delete tree_var;
  delete label;
  
  for (Int_t iHisto=0;iHisto<NHistos;iHisto++) {
    delete h[iHisto];
    delete r[iHisto];
  }

}
