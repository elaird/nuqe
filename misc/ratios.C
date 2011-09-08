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
void adjust_histo(TH1D*h,Int_t color) {
  //h->SetStats(kFALSE);
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleOffset(1.2);
  h->SetLineColor(color);
  h->SetMarkerColor(color);
}

////////////////////////////////////////////////////////////////////////
void ratios_look(Int_t mode,Int_t energy_point) {
  Init();

  //NUANCE mode
  TFile f4mu(Form("../../events/events_4_mu_%d.root",energy_point));
  gROOT->cd();
  TTree *tree4mu=(TTree*)f4mu.Get("tree");
  tree4mu->SetName("tree4mu");
  TGraph *gr4mu=(TGraph*)f4mu.Get("ccqe_rate")->Clone();
  Double_t rate4mu,dummy;
  gr4mu->GetPoint(4,dummy,rate4mu);
  delete gr4mu;

  TFile f4e(Form("../../events/events_4_e_%d.root",energy_point));
  gROOT->cd();
  TTree *tree4e=(TTree*)f4e.Get("tree");
  tree4e->SetName("tree4e");
  TGraph *gr4e=(TGraph*)f4e.Get("ccqe_rate")->Clone();
  Double_t rate4e;
  gr4e->GetPoint(4,dummy,rate4e);
  delete gr4e;

  //A-S SF mode
  TFile f23mu(Form("../../events/events_23_mu_%d.root",energy_point));
  gROOT->cd();
  TTree *tree23mu=(TTree*)f23mu.Get("tree");
  tree23mu->SetName("tree23mu");
  TGraph *gr23mu=(TGraph*)f23mu.Get("ccqe_rate")->Clone();
  Double_t rate23mu=0.0;
  Double_t dummy_x,dummy_y;
  gr23mu->GetPoint(2,dummy_x,dummy_y);
  rate23mu+=dummy_y;
  gr23mu->GetPoint(3,dummy_x,dummy_y);
  rate23mu+=dummy_y;
  delete gr23mu;
  
  TFile f23e(Form("../../events/events_23_e_%d.root",energy_point));
  gROOT->cd();
  TTree *tree23e=(TTree*)f23e.Get("tree");
  tree23e->SetName("tree23e");
  TGraph *gr23e=(TGraph*)f23e.Get("ccqe_rate")->Clone();
  Double_t rate23e=0.0;
  gr23e->GetPoint(2,dummy_x,dummy_y);
  rate23e+=dummy_y;
  gr23e->GetPoint(3,dummy_x,dummy_y);
  rate23e+=dummy_y;
  delete gr23e;
  
  //A-S SM mode
  TFile f1mu(Form("../../events/events_1_mu_%d.root",energy_point));
  gROOT->cd();
  TTree *tree1mu=(TTree*)f1mu.Get("tree");
  tree1mu->SetName("tree1mu");
  TGraph *gr1mu=(TGraph*)f1mu.Get("ccqe_rate")->Clone();
  Double_t rate1mu=0.0;
  gr1mu->GetPoint(1,dummy_x,rate1mu);
  delete gr1mu;
  
  TFile f1e(Form("../../events/events_1_e_%d.root",energy_point));
  gROOT->cd();
  TTree *tree1e=(TTree*)f1e.Get("tree");
  tree1e->SetName("tree1e");
  TGraph *gr1e=(TGraph*)f1e.Get("ccqe_rate")->Clone();
  Double_t rate1e=0.0;
  gr1e->GetPoint(1,dummy_x,rate1e);
  delete gr1e;
  
  TString *tree_var=0;
  TString *label=0;

  Int_t N_bins=0;
  Double_t lower=0.0;
  Double_t upper=0.0;
  Double_t x_leg=0.0;

  switch (mode) {
  case 1:
    tree_var=new TString("k[0]");
    label=new TString(";E_{#nu} (GeV);(<events>/POT/target) / GeV");

    x_leg=0.15;
    N_bins=150;
    lower=0.0;
    upper=2.5;
    break;
  case 2:
    tree_var=new TString("enuqe");
    label=new TString(";E_{#nu}^{QE} (GeV);(<events>/POT/target) / GeV");

    x_leg=0.15;
    N_bins=150;
    lower=0.0;
    upper=2.5;
    break;
  }

  TH1D h4mu ("h4mu" ,*label,N_bins,lower,upper); h4mu.Sumw2();
  TH1D h4e  ("h4e"  ,*label,N_bins,lower,upper); h4e.Sumw2();
  TH1D h23mu("h23mu",*label,N_bins,lower,upper); h23mu.Sumw2();
  TH1D h23e ("h23e" ,*label,N_bins,lower,upper); h23e.Sumw2();
  TH1D h1mu ("h1mu" ,*label,N_bins,lower,upper); h1mu.Sumw2();
  TH1D h1e  ("h1e"  ,*label,N_bins,lower,upper); h1e.Sumw2();

  TH1D r4 ("r4" ,*label,N_bins,lower,upper); r4.Sumw2();
  TH1D r23("r23",*label,N_bins,lower,upper); r23.Sumw2();
  TH1D r1 ("r1" ,*label,N_bins,lower,upper); r1.Sumw2();
  TH1D r_r23_r4("r_r23_r4",*label,N_bins,lower,upper); r_r23_r4.Sumw2();
  r4.GetYaxis()->SetTitle("ratio");
  r23.GetYaxis()->SetTitle("ratio");
  r1.GetYaxis()->SetTitle("ratio");
  r_r23_r4.GetYaxis()->SetTitle("ratio");

  tree4mu ->Draw(*tree_var+">>h4mu" ,"process==4","goff");
  tree4e  ->Draw(*tree_var+">>h4e"  ,"process==4","goff");
  tree23mu->Draw(*tree_var+">>h23mu","process==2 || process==3","goff");
  tree23e ->Draw(*tree_var+">>h23e" ,"process==2 || process==3","goff");
  tree1mu ->Draw(*tree_var+">>h1mu" ,"process==1","goff");
  tree1e  ->Draw(*tree_var+">>h1e"  ,"process==1","goff");

  h4mu.Scale(rate4mu/h4mu.Integral("width"));
  h4e.Scale(rate4e/h4e.Integral("width"));

  h23mu.Scale(rate23mu/h23mu.Integral("width"));
  h23e.Scale(rate23e/h23e.Integral("width"));

  h1mu.Scale(rate1mu/h1mu.Integral("width"));
  h1e.Scale(rate1e/h1e.Integral("width"));

  r23.Divide(&h23mu,&h23e);
  r4.Divide(&h4mu,&h4e);
  r1.Divide(&h1mu,&h1e);
  r_r23_r4.Divide(&r23,&r4);

  r23.SetMinimum(0.0);
  r23.SetMaximum(2.0);
  r4.SetMinimum(0.0);
  r4.SetMaximum(2.0);
  r1.SetMinimum(0.0);
  r1.SetMaximum(2.0);
  r_r23_r4.SetMinimum(0.0);
  r_r23_r4.SetMaximum(2.0);

  adjust_histo(&h4mu,kBlack);
  adjust_histo(&h4e,kRed);

  adjust_histo(&h23mu,kBlue);
  adjust_histo(&h23e,kGreen);

  adjust_histo(&h1mu,kMagenta);
  adjust_histo(&h1e,kCyan);

  adjust_histo(&r4,kBlack);
  adjust_histo(&r23,kBlue);
  adjust_histo(&r1,kMagenta);
  adjust_histo(&r_r23_r4,kBlack);


//TLegend legend(x_leg,0.6,x_leg+0.35,0.9);
//legend.SetBorderSize(0);
//legend.SetFillStyle(0);
//legend.SetHeader(Form("dipole form factors, O^{16}, %d MeV #nu_{#mu}",energy_point));
//legend.AddEntry(&p4,"NUANCE mode, MA=1.23 GeV");
//legend.AddEntry(&p1,"A.-S. SM mode, MA=1.23 GeV");
//

  TLine line;
  TString psfile=Form("new_ratios_%d.ps",mode);
  TCanvas canvas("canvas","",700,1000);
  canvas.Print(psfile+"[","Landscape");
  canvas.Divide(1,2);

  canvas.cd(1);
  h4mu.Draw();
  h4e.Draw("same");
  canvas.cd(2);
  r4.Draw();
  line.DrawLine(lower,1.0,upper,1.0);
  canvas.cd(0);
  canvas.Print(psfile,"Landscape");
  
  canvas.cd(1);
  h23mu.Draw();
  h23e.Draw("same");
  canvas.cd(2);
  r23.Draw();
  line.DrawLine(lower,1.0,upper,1.0);
  canvas.cd(0);
  canvas.Print(psfile,"Landscape");
  
  canvas.cd(1);
  h1mu.Draw();
  h1e.Draw("same");
  canvas.cd(2);
  r1.Draw();
  line.DrawLine(lower,1.0,upper,1.0);
  canvas.cd(0);
  canvas.Print(psfile,"Landscape");
  
  canvas.Divide(1,1);
  canvas.cd(0);
  r_r23_r4.Draw();
  line.DrawLine(lower,1.0,upper,1.0);
  canvas.Print(psfile,"Landscape");
  
  canvas.Print(psfile+"]","Landscape");

  gROOT->ProcessLine(".! gv "+psfile+" &");
  delete tree_var;
  delete label;

}

////////////////////////////////////////////////////////////////////////
void nu_antinu_ratio_look(Int_t mode,Int_t energy_point) {
  Init();

  //NUANCE mode
  TFile f4mu(Form("../../events/events_4_mu_%d.root",energy_point));
  gROOT->cd();
  TTree *tree4mu=(TTree*)f4mu.Get("tree");
  tree4mu->SetName("tree4mu");
  TGraph *gr4mu=(TGraph*)f4mu.Get("ccqe_rate")->Clone();
  Double_t rate4mu,dummy;
  gr4mu->GetPoint(4,dummy,rate4mu);
  delete gr4mu;

  TFile f4mubar(Form("../../events/events_4_mubar_%d.root",energy_point));
  gROOT->cd();
  TTree *tree4mubar=(TTree*)f4mubar.Get("tree");
  tree4mubar->SetName("tree4mubar");
  TGraph *gr4mubar=(TGraph*)f4mubar.Get("ccqe_rate")->Clone();
  Double_t rate4mubar;
  gr4mubar->GetPoint(4,dummy,rate4mubar);
  delete gr4mubar;

  //A-S SF mode
  TFile f23mu(Form("../../events/events_23_mu_%d.root",energy_point));
  gROOT->cd();
  TTree *tree23mu=(TTree*)f23mu.Get("tree");
  tree23mu->SetName("tree23mu");
  TGraph *gr23mu=(TGraph*)f23mu.Get("ccqe_rate")->Clone();
  Double_t rate23mu=0.0;
  Double_t dummy_x,dummy_y;
  gr23mu->GetPoint(2,dummy_x,dummy_y);
  rate23mu+=dummy_y;
  gr23mu->GetPoint(3,dummy_x,dummy_y);
  rate23mu+=dummy_y;
  delete gr23mu;
  
  TFile f23mubar(Form("../../events/events_23_mubar_%d.root",energy_point));
  gROOT->cd();
  TTree *tree23mubar=(TTree*)f23mubar.Get("tree");
  tree23mubar->SetName("tree23mubar");
  TGraph *gr23mubar=(TGraph*)f23mubar.Get("ccqe_rate")->Clone();
  Double_t rate23mubar=0.0;
  gr23mubar->GetPoint(2,dummy_x,dummy_y);
  rate23mubar+=dummy_y;
  gr23mubar->GetPoint(3,dummy_x,dummy_y);
  rate23mubar+=dummy_y;
  delete gr23mubar;
  
  //A-S SM mode
  TFile f1mu(Form("../../events/events_1_mu_%d.root",energy_point));
  gROOT->cd();
  TTree *tree1mu=(TTree*)f1mu.Get("tree");
  tree1mu->SetName("tree1mu");
  TGraph *gr1mu=(TGraph*)f1mu.Get("ccqe_rate")->Clone();
  Double_t rate1mu=0.0;
  gr1mu->GetPoint(1,dummy_x,rate1mu);
  delete gr1mu;
  
  TFile f1mubar(Form("../../events/events_1_mubar_%d.root",energy_point));
  gROOT->cd();
  TTree *tree1mubar=(TTree*)f1mubar.Get("tree");
  tree1mubar->SetName("tree1mubar");
  TGraph *gr1mubar=(TGraph*)f1mubar.Get("ccqe_rate")->Clone();
  Double_t rate1mubar=0.0;
  gr1mubar->GetPoint(1,dummy_x,rate1mubar);
  delete gr1mubar;
  
  TString *tree_var=0;
  TString *label=0;

  Int_t N_bins=0;
  Double_t lower=0.0;
  Double_t upper=0.0;
  Double_t x_leg=0.0;

  switch (mode) {
  case 1:
    tree_var=new TString("k[0]");
    label=new TString(";E_{#nu} (GeV);(<events>/POT/target) / GeV");

    x_leg=0.15;
    N_bins=150;
    lower=0.0;
    upper=2.5;
    break;
  case 2:
    tree_var=new TString("enuqe");
    label=new TString(";E_{#nu}^{QE} (GeV);(<events>/POT/target) / GeV");

    x_leg=0.15;
    N_bins=150;
    lower=0.0;
    upper=2.5;
    break;
  }

  TH1D h4mu ("h4mu" ,*label,N_bins,lower,upper); h4mu.Sumw2();
  TH1D h4mubar  ("h4mubar"  ,*label,N_bins,lower,upper); h4mubar.Sumw2();
  TH1D h23mu("h23mu",*label,N_bins,lower,upper); h23mu.Sumw2();
  TH1D h23mubar ("h23mubar" ,*label,N_bins,lower,upper); h23mubar.Sumw2();
  TH1D h1mu ("h1mu" ,*label,N_bins,lower,upper); h1mu.Sumw2();
  TH1D h1mubar  ("h1mubar"  ,*label,N_bins,lower,upper); h1mubar.Sumw2();

  TH1D r4 ("r4" ,*label,N_bins,lower,upper); r4.Sumw2();
  TH1D r23("r23",*label,N_bins,lower,upper); r23.Sumw2();
  TH1D r1 ("r1" ,*label,N_bins,lower,upper); r1.Sumw2();
  TH1D r_r23_r4("r_r23_r4",*label,N_bins,lower,upper); r_r23_r4.Sumw2();
  r4.GetYaxis()->SetTitle("ratio");
  r23.GetYaxis()->SetTitle("ratio");
  r1.GetYaxis()->SetTitle("ratio");
  r_r23_r4.GetYaxis()->SetTitle("ratio");

  tree4mu ->Draw(*tree_var+">>h4mu" ,"process==4","goff");
  tree4mubar  ->Draw(*tree_var+">>h4mubar"  ,"process==4","goff");
  tree23mu->Draw(*tree_var+">>h23mu","process==2 || process==3","goff");
  tree23mubar ->Draw(*tree_var+">>h23mubar" ,"process==2 || process==3","goff");
  tree1mu ->Draw(*tree_var+">>h1mu" ,"process==1","goff");
  tree1mubar  ->Draw(*tree_var+">>h1mubar"  ,"process==1","goff");

  h4mu.Scale(rate4mu/h4mu.Integral("width"));
  h4mubar.Scale(rate4mubar/h4mubar.Integral("width"));

  h23mu.Scale(rate23mu/h23mu.Integral("width"));
  h23mubar.Scale(rate23mubar/h23mubar.Integral("width"));

  h1mu.Scale(rate1mu/h1mu.Integral("width"));
  h1mubar.Scale(rate1mubar/h1mubar.Integral("width"));

  r23.Divide(&h23mu,&h23mubar);
  r4.Divide(&h4mu,&h4mubar);
  r1.Divide(&h1mu,&h1mubar);
  r_r23_r4.Divide(&r23,&r4);

  r23.SetMinimum(1.0);
  r23.SetMaximum(6.0);
  r4.SetMinimum(1.0);
  r4.SetMaximum(6.0);
  r1.SetMinimum(1.0);
  r1.SetMaximum(6.0);
  r_r23_r4.SetMinimum(0.0);
  r_r23_r4.SetMaximum(2.0);

  adjust_histo(&h4mu,kBlack);
  adjust_histo(&h4mubar,kRed);

  adjust_histo(&h23mu,kBlue);
  adjust_histo(&h23mubar,kGreen);

  adjust_histo(&h1mu,kMagenta);
  adjust_histo(&h1mubar,kCyan);

  adjust_histo(&r4,kBlack);
  adjust_histo(&r23,kBlue);
  adjust_histo(&r1,kMagenta);
  adjust_histo(&r_r23_r4,kBlack);


//TLegend legend(x_leg,0.6,x_leg+0.35,0.9);
//legend.SetBorderSize(0);
//legend.SetFillStyle(0);
//legend.SetHeader(Form("dipole form factors, O^{16}, %d MeV #nu_{#mu}",energy_point));
//legend.AddEntry(&p4,"NUANCE mode, MA=1.23 GeV");
//legend.AddEntry(&p1,"A.-S. SM mode, MA=1.23 GeV");
//

  TLine line;
  TString psfile=Form("new_ratios_%d.ps",mode);
  TCanvas canvas("canvas","",700,1000);
  canvas.Print(psfile+"[","Landscape");
  canvas.Divide(1,2);

  canvas.cd(1);
  h4mu.Draw();
  h4mubar.Draw("same");
  canvas.cd(2);
  r4.Draw();
  line.DrawLine(lower,1.0,upper,1.0);
  canvas.cd(0);
  canvas.Print(psfile,"Landscape");
  
  canvas.cd(1);
  h23mu.Draw();
  h23mubar.Draw("same");
  canvas.cd(2);
  r23.Draw();
  line.DrawLine(lower,1.0,upper,1.0);
  canvas.cd(0);
  canvas.Print(psfile,"Landscape");
  
  canvas.cd(1);
  h1mu.Draw();
  h1mubar.Draw("same");
  canvas.cd(2);
  r1.Draw();
  line.DrawLine(lower,1.0,upper,1.0);
  canvas.cd(0);
  canvas.Print(psfile,"Landscape");
  
  canvas.Divide(1,1);
  canvas.cd(0);
  r_r23_r4.Draw();
  line.DrawLine(lower,1.0,upper,1.0);
  canvas.Print(psfile,"Landscape");
  
  canvas.Print(psfile+"]","Landscape");

  gROOT->ProcessLine(".! gv "+psfile+" &");
  delete tree_var;
  delete label;

}

////////////////////////////////////////////////////////////////////////
void mag_p_look(Int_t mode,Int_t energy_point) {
  Init();

  //A-S SF mode
  TFile f23mu(Form("../../events/events_23_mu_%d.root",energy_point));
  gROOT->cd();
  TTree *tree23mu=(TTree*)f23mu.Get("tree");
  tree23mu->SetName("tree23mu");
  TGraph *gr23mu=(TGraph*)f23mu.Get("ccqe_rate")->Clone();
  Double_t rate23mu=0.0;
  Double_t dummy_x,dummy_y;
  gr23mu->GetPoint(2,dummy_x,dummy_y);
  rate23mu+=dummy_y;
  gr23mu->GetPoint(3,dummy_x,dummy_y);
  rate23mu+=dummy_y;
  delete gr23mu;
  
  //A-S SM mode
  TFile f1mu(Form("../../events/events_1_mu_%d.root",energy_point));
  gROOT->cd();
  TTree *tree1mu=(TTree*)f1mu.Get("tree");
  tree1mu->SetName("tree1mu");
  TGraph *gr1mu=(TGraph*)f1mu.Get("ccqe_rate")->Clone();
  Double_t rate1mu=0.0;
  gr1mu->GetPoint(1,dummy_x,rate1mu);
  delete gr1mu;
  
  TString *tree_var=0;
  TString *label=0;

  Int_t N_bins=0;
  Double_t lower=0.0;
  Double_t upper=0.0;
  Double_t x_leg=0.0;

  switch (mode) {
  case 1:
    tree_var=new TString("k[0]");
    label=new TString(";E_{#nu} (GeV);(<events>/POT/target) / GeV");

    x_leg=0.15;
    N_bins=150;
    lower=0.0;
    upper=2.5;
    break;
  case 2:
    tree_var=new TString("enuqe");
    label=new TString(";E_{#nu}^{QE} (GeV);(<events>/POT/target) / GeV");

    x_leg=0.15;
    N_bins=150;
    lower=0.0;
    upper=2.5;
    break;
  case 3:
    tree_var=new TString("mag_p");
    label=new TString(";p_{i} (GeV/c);(<events>/POT/target) / (GeV/c)");

    x_leg=0.15;
    N_bins=100;
    lower=0.0;
    upper=0.5;
    break;
  }

  TH1D h23mu("h23mu",*label,N_bins,lower,upper); h23mu.Sumw2();
  TH1D h1mu ("h1mu" ,*label,N_bins,lower,upper); h1mu.Sumw2();

  TH1D r23_1("r23_1",*label,N_bins,lower,upper); r23_1.Sumw2();
  r23_1.GetYaxis()->SetTitle("ratio");

  tree23mu->Draw(*tree_var+">>h23mu","(process==2 || process==3)","goff");
  tree1mu ->Draw(*tree_var+">>h1mu" ,"(process==1)","goff");

  h23mu.Scale(rate23mu/h23mu.Integral("width"));
  h1mu.Scale(rate1mu/h1mu.Integral("width"));

  r23_1.Divide(&h23mu,&h1mu);
  r23_1.SetMinimum(0.0);
  r23_1.SetMaximum(2.0);

  adjust_histo(&h23mu,kBlue);
  adjust_histo(&h1mu,kMagenta);
  adjust_histo(&r23_1,kBlue);

//TLegend legend(x_leg,0.6,x_leg+0.35,0.9);
//legend.SetBorderSize(0);
//legend.SetFillStyle(0);
//legend.SetHeader(Form("dipole form factors, O^{16}, %d MeV #nu_{#mu}",energy_point));
//legend.AddEntry(&p4,"NUANCE mode, MA=1.23 GeV");
//legend.AddEntry(&p1,"A.-S. SM mode, MA=1.23 GeV");
//

  TLine line;
  TString psfile=Form("new_ratios_%d.ps",mode);
  TCanvas canvas("canvas","",700,1000);
  canvas.Print(psfile+"[","Landscape");

  h1mu.Draw();
  h23mu.Draw("same");
  canvas.Print(psfile,"Landscape");
  
  r23_1.Draw();
  canvas.Print(psfile,"Landscape");
  
  canvas.Print(psfile+"]","Landscape");

  gROOT->ProcessLine(".! gv "+psfile+" &");
  delete tree_var;
  delete label;

}
