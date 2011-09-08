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
#include "benhar_rip/BenharDrawer.h"
#include "butkevich_numbers/ButkevichDrawer.h"

////////////////////////////////////////////////////////////////////////
class var_stuff {
public :
  
  TString tree_var;
  TString label;
  double  xLeg;
  int     nBins;
  double  lower;
  double  upper;

//  var_stuff();
//  virtual ~var_stuff();
};
////////////////////////////////////////////////////////////////////////
class MA_looker {
public :

  MA_looker(int energyPoint);
  virtual ~MA_looker();
  virtual void adjust_histo(TH1& histo);
  virtual void XS_look();
  virtual void compare_simple();
  virtual void make_benhar_and_nuance_comparisons();
  virtual void make_ratio(TH1D *ratio,TGraph& graph);
  virtual void makePsFile();

private :
  TCanvas* f_canvas;
  TString f_rootFile123_ma103;
  TString f_rootFile123_ma123;
  TString f_rootFile4_ma103;
  TString f_rootFile4_ma123;
  TString f_psFile;
  TString f_psOptions;
  int     f_energyPoint;
  double  f_nEvents;

  double  f_xTitleSize;
  double  f_xTitleOffset;
  double  f_yTitleSize;
  double  f_yTitleOffset;
  
  double f_xLegRightSide;
  double f_xLegLeftSide;

  std::vector<var_stuff> f_varStuffs;
  

};
////////////////////////////////////////////////////////////////////////
MA_looker::MA_looker(int energyPoint) {
  f_energyPoint=energyPoint;
  f_nEvents=1.0e10;
  //f_nEvents=1.0e4;

  gROOT->SetStyle("Plain");
  f_xTitleSize=0.06;
  f_xTitleOffset=0.75;
  f_yTitleSize=0.06;
  f_yTitleOffset=0.75;

  f_xLegRightSide=0.53;
  f_xLegLeftSide=0.15;

  f_canvas = new TCanvas("f_canvas","",700,1000);
  f_rootFile123_ma103="~ted/work/boone/ccqe/events/may09_events_2M_proc123_mu";
  f_rootFile123_ma123="~ted/work/boone/ccqe/events/may09_events_2M_proc123_ma123_mu";
  f_rootFile4_ma103="~ted/work/boone/ccqe/events/may09_events_2M_proc4_mu";
  f_rootFile4_ma123="~ted/work/boone/ccqe/events/may09_events_2M_proc4_ma123_mu";
  f_psFile="may09.ps";
  f_psOptions="";

  bool drawQ2    = true;
  bool drawQ2QE  = true;
  bool drawEnu   = true;
  bool drawEnuQE = true;
  bool drawEmu   = true;
  bool drawUz    = true;
  bool drawPi    = true;
  bool drawPf    = true;

  {
    if (drawQ2) {
      var_stuff foo0;                                            var_stuff foo1;						      
      foo0.tree_var="-q[0]**2+q[1]**2+q[2]**2+q[3]**2";          foo1.tree_var="-q[0]**2+q[1]**2+q[2]**2+q[3]**2";
      foo0.label=";Q^{2} (GeV^{2});event rate (arb. units)";     foo1.label=";Q^{2} (GeV^{2});d#sigma/dQ^{2} (cm^{2}/GeV^{2})";
      foo0.xLeg=f_xLegRightSide;				 foo1.xLeg=f_xLegRightSide;						      
      foo0.nBins=70;					         foo1.nBins=70;						      
      foo0.lower=0.0;					         foo1.lower=0.0;						      
      foo0.upper=1.4;                                            foo1.upper=1.4;
      
      if (f_energyPoint==0) f_varStuffs.push_back(foo0);
      else                  f_varStuffs.push_back(foo1);
    }
  }
  {
    if (drawQ2QE) {
      var_stuff foo0;                                            var_stuff foo1;						      
      foo0.tree_var="Q2qe";                                      foo1.tree_var="Q2qe";
      foo0.label=";Q^{2}_{QE} (GeV^{2});event rate (arb. units)";foo1.label=";Q^{2}_{QE} (GeV^{2});d#sigma/dQ^{2}_{QE} (cm^{2}/GeV^{2})";
      foo0.xLeg=f_xLegRightSide;				 foo1.xLeg=f_xLegRightSide;						      
      foo0.nBins=70;					         foo1.nBins=70;						      
      foo0.lower=0.0;					         foo1.lower=0.0;						      
      foo0.upper=1.4;                                            foo1.upper=1.4;
      
      if (f_energyPoint==0) f_varStuffs.push_back(foo0);
      else                  f_varStuffs.push_back(foo1);
    }
  }
  {
    if (drawEnu) {
      var_stuff foo0;                                          var_stuff foo1;                                                 
      foo0.tree_var="k[0]";				       foo1.tree_var="k[0]";					   
      foo0.label=";E_{#nu} (GeV);event rate (arb. units)";     foo1.label=";E_{#nu} (GeV);d#sigma/dE_{#nu} (cm^{2}/GeV)";	   
      foo0.xLeg=f_xLegRightSide;			       foo1.xLeg=f_xLegRightSide;						   
      foo0.nBins=70;					       foo1.nBins=70;						   
      foo0.lower=0.0;					       foo1.lower=0.0;						   
      foo0.upper=2.5;					       foo1.upper=2.5;						   
      
      if (f_energyPoint==0) f_varStuffs.push_back(foo0);
      else                  f_varStuffs.push_back(foo1);
    }
  }							       
  {
    if (drawEnuQE) {
      var_stuff foo0;                                            var_stuff foo1;                                                 
      foo0.tree_var="enuqe";				         foo1.tree_var="enuqe";					   
      foo0.label=";E_{#nu}^{QE} (GeV);event rate (arb. units)";  foo1.label=";E_{#nu}^{QE} (GeV);d#sigma/dE_{#nu}^{QE} (cm^{2}/GeV)";	   
      foo0.xLeg=f_xLegRightSide;				 foo1.xLeg=f_xLegRightSide;						   
      foo0.nBins=70;					         foo1.nBins=70;						   
      foo0.lower=0.0;					         foo1.lower=0.0;						   
      foo0.upper=2.5;					         foo1.upper=2.5;						   
      
      if (f_energyPoint==0) f_varStuffs.push_back(foo0);
      else                  f_varStuffs.push_back(foo1);
    }
  }							       
  {
    if (drawEmu) {
      var_stuff foo0;                                          var_stuff foo1;						      
      foo0.tree_var="kprime[0]";			       foo1.tree_var="kprime[0]";					      
      foo0.label=";E_{#mu} (GeV);event rate (arb. units)";     foo1.label=";E_{#mu} (GeV);d#sigma/dE_{#mu} (cm^{2}/GeV)";   
      foo0.xLeg=f_xLegRightSide;			       foo1.xLeg=f_xLegLeftSide;						      
      foo0.nBins=70;					       foo1.nBins=70;						      
      foo0.lower=0.0;					       foo1.lower=0.0;						      
      foo0.upper=2.5;                                          foo1.upper=f_energyPoint*1.0e-3;                                   
      
      if (f_energyPoint==0) f_varStuffs.push_back(foo0);
      else                  f_varStuffs.push_back(foo1);
    }
  }
  {
    if (drawUz) {
      var_stuff foo0;                                                         var_stuff foo1;						      
      foo0.tree_var="kprime[3]/sqrt(kprime[1]**2+kprime[2]**2+kprime[3]**2)"; foo1.tree_var="kprime[3]/sqrt(kprime[1]**2+kprime[2]**2+kprime[3]**2)";
      foo0.label=";cos(#theta_{#mu});event rate (arb. units)";                foo1.label=";cos(#theta_{#mu};d#sigma/dcos(#theta_{#mu}) (cm^{2})";
      foo0.xLeg=f_xLegLeftSide;					              foo1.xLeg=f_xLegLeftSide;						      
      foo0.nBins=70;					                      foo1.nBins=150;						      
      foo0.lower=-1.0;					                      foo1.lower=-1.0;						      
      foo0.upper= 1.0;                                                        foo1.upper= 1.0;
      
      if (f_energyPoint==0) f_varStuffs.push_back(foo0);
      else                  f_varStuffs.push_back(foo1);
    }
  }
  {
    if (drawPi) {
      var_stuff foo0;                                             var_stuff foo1;						      
      foo0.tree_var="sqrt(p[1]**2+p[2]**2+p[3]**2)";              foo1.tree_var="sqrt(p[1]**2+p[2]**2+p[3]**2)";
      foo0.label=";p_{initial} (GeV/c);event rate (arb. units)";  foo1.label=";p_{initial} (GeV/c);d#sigma/dp_{initial} (cm^{2}/(GeV/c))";
      foo0.xLeg=f_xLegRightSide;				  foo1.xLeg=f_xLegRightSide;						      
      foo0.nBins=100;					          foo1.nBins=100;						      
      foo0.lower=0.0;					          foo1.lower=0.0;						      
      foo0.upper=0.5;                                             foo1.upper=0.5;
      
      if (f_energyPoint==0) f_varStuffs.push_back(foo0);
      else                  f_varStuffs.push_back(foo1);
    }
  }
  {
    if (drawPf) {
      var_stuff foo0;                                               var_stuff foo1;						      
      foo0.tree_var="sqrt(pprime[1]**2+pprime[2]**2+pprime[3]**2)"; foo1.tree_var="sqrt(pprime[1]**2+pprime[2]**2+pprime[3]**2)";
      foo0.label=";p_{final} (GeV/c);event rate (arb. units)";      foo1.label=";p_{final} (GeV/c);d#sigma/dp_{final} (cm^{2}/(GeV/c))";
      foo0.xLeg=f_xLegRightSide;				    foo1.xLeg=f_xLegRightSide;						      
      foo0.nBins=100;					            foo1.nBins=100;						      
      foo0.lower=0.1;					            foo1.lower=0.1;						      
      foo0.upper=1.6;                                               foo1.upper=1.6;
      
      if (f_energyPoint==0) f_varStuffs.push_back(foo0);
      else                  f_varStuffs.push_back(foo1);
    }
  }

}
////////////////////////////////////////////////////////////////////////
MA_looker::~MA_looker() {
  delete f_canvas;
}
////////////////////////////////////////////////////////////////////////
void MA_looker::adjust_histo(TH1& histo) {
  histo.SetStats(kFALSE);
  histo.GetXaxis()->CenterTitle();
  histo.GetXaxis()->SetTitleSize(f_xTitleSize);
  histo.GetXaxis()->SetTitleOffset(f_xTitleOffset);
  histo.GetYaxis()->CenterTitle();
  histo.GetYaxis()->SetTitleSize(f_yTitleSize);
  histo.GetYaxis()->SetTitleOffset(f_yTitleOffset);
}
////////////////////////////////////////////////////////////////////////
void MA_looker::compare_simple() {

  //A-S SF mode
  TFile f123_ma103(f_rootFile123_ma103+Form("_%d.root",f_energyPoint));
  gROOT->cd();
  TTree *tree123_ma103=(TTree*)f123_ma103.Get("tree");
  tree123_ma103->SetName("tree123_ma103");
  Double_t rate1_ma103=0.0;
  Double_t rate23_ma103=0.0;
  {
    TGraph *gr=(TGraph*)f123_ma103.Get("ccqe_rate")->Clone();
    Double_t dummy_x,dummy_y;
    gr->GetPoint(1,dummy_x,rate1_ma103);
    gr->GetPoint(2,dummy_x,dummy_y);
    rate23_ma103+=dummy_y;
    gr->GetPoint(3,dummy_x,dummy_y);
    rate23_ma103+=dummy_y;
    delete gr;
  }

  TFile f123_ma123(f_rootFile123_ma123+Form("_%d.root",f_energyPoint));
  gROOT->cd();
  TTree *tree123_ma123=(TTree*)f123_ma123.Get("tree");
  tree123_ma123->SetName("tree123_ma123");
  Double_t rate1_ma123=0.0;
  Double_t rate23_ma123=0.0;
  {
    TGraph *gr=(TGraph*)f123_ma123.Get("ccqe_rate")->Clone();
    Double_t dummy_x,dummy_y;
    gr->GetPoint(1,dummy_x,rate1_ma123);
    gr->GetPoint(2,dummy_x,dummy_y);
    rate23_ma123+=dummy_y;
    gr->GetPoint(3,dummy_x,dummy_y);
    rate23_ma123+=dummy_y;
    delete gr;
  }

  //NUANCE
  TFile f4_ma103(f_rootFile4_ma103+Form("_%d.root",f_energyPoint));
  gROOT->cd();
  TTree *tree4_ma103=(TTree*)f4_ma103.Get("tree");
  tree4_ma103->SetName("tree4_ma103");
  Double_t rate4_ma103=0.0;
  {
    TGraph *gr=(TGraph*)f4_ma103.Get("ccqe_rate")->Clone();
    Double_t dummy_x,dummy_y;
    gr->GetPoint(4,dummy_x,rate4_ma103);
    delete gr;
  }
  TFile f4_ma123(f_rootFile4_ma123+Form("_%d.root",f_energyPoint));
  gROOT->cd();
  TTree *tree4_ma123=(TTree*)f4_ma123.Get("tree");
  tree4_ma123->SetName("tree4_ma123");
  Double_t rate4_ma123=0.0;
  {
    TGraph *gr=(TGraph*)f4_ma123.Get("ccqe_rate")->Clone();
    Double_t dummy_x,dummy_y;
    gr->GetPoint(4,dummy_x,rate4_ma123);
    delete gr;
  }

  f_canvas->cd(0);
  f_canvas->Divide(1,3);

  for (unsigned int iVar=0;iVar<f_varStuffs.size();iVar++) {
    TString tree_var=f_varStuffs.at(iVar).tree_var;
    TString label=f_varStuffs.at(iVar).label;
    Int_t N_bins=f_varStuffs.at(iVar).nBins;
    Double_t lower=f_varStuffs.at(iVar).lower;
    Double_t upper=f_varStuffs.at(iVar).upper;
    Double_t xleg=f_varStuffs.at(iVar).xLeg;


    TH1D p23_ma103("p23_ma103",label,N_bins,lower,upper);
    p23_ma103.Sumw2();
    p23_ma103.SetStats(kFALSE);
    tree123_ma103->Draw(tree_var+">>p23_ma103","(process==2 || process==3)","goff",f_nEvents,1);
    p23_ma103.Scale(rate23_ma103/p23_ma103.Integral("width"));
    p23_ma103.SetLineColor(kRed);
    p23_ma103.SetMarkerColor(kRed);

    TH1D p1("p1",label,N_bins,lower,upper);
    p1.Sumw2();
    p1.SetStats(kFALSE);
    tree123_ma103->Draw(tree_var+">>p1","process==1","goff",f_nEvents,1);
    p1.Scale(rate1_ma103/p1.Integral("width"));
    p1.SetLineColor(kBlue);
    p1.SetMarkerColor(kBlue);

    TH1D p4_ma103("p4_ma103",label,N_bins,lower,upper);
    p4_ma103.Sumw2();
    p4_ma103.SetStats(kFALSE);
    tree4_ma103->Draw(tree_var+">>p4_ma103","process==4","goff",f_nEvents,1);
    p4_ma103.Scale(rate4_ma103/p4_ma103.Integral("width"));
    p4_ma103.SetLineColor(kBlack);
    p4_ma103.SetMarkerColor(kBlack);

    TH1D p4_ma123("p4_ma123",label,N_bins,lower,upper);
    p4_ma123.Sumw2();
    p4_ma123.SetStats(kFALSE);
    tree4_ma123->Draw(tree_var+">>p4_ma123","process==4","goff",f_nEvents,1);
    p4_ma123.Scale(rate4_ma123/p4_ma123.Integral("width"));
    p4_ma123.SetLineColor(kMagenta);
    p4_ma123.SetMarkerColor(kMagenta);

    TLegend legend(xleg,0.53,xleg+0.48,0.88);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    TString headerString="dipole form factors; O^{16}; ";
    if (f_energyPoint==0) headerString+=Form("#nu_{#mu}; MB flux",f_energyPoint);
    else headerString+=Form("%d MeV #nu_{#mu}",f_energyPoint);
    legend.SetHeader(headerString);
    legend.AddEntry(&p23_ma103,"A.-S. SF mode, MA=1.03 GeV");
    legend.AddEntry(&p1       ,"A.-S. SM mode, MA=1.03 GeV");
    legend.AddEntry(&p4_ma103,"NUANCE mode, MA=1.03 GeV");
    legend.AddEntry(&p4_ma123,"NUANCE mode, MA=1.23 GeV");
    
    f_canvas->cd(1);
    adjust_histo(p4_ma123);
    p4_ma123.Draw("e");
    p4_ma103.Draw("esame");
    p1.Draw("esame");
    p23_ma103.Draw("esame");
    
    legend.Draw("same");
  
    f_canvas->cd(2);
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetTickx();
    gPad->SetTicky();

    TString xtitle=p4_ma123.GetXaxis()->GetTitle();

    TH2D null_ma103("null_ma103",";"+xtitle+Form(";ratio to (NUANCE mode @ 1.03)"),1,lower,upper,1,0.5,1.5);
    adjust_histo(null_ma103);
    null_ma103.Draw();
    
    TH1D* ratio_23_4_ma103=(TH1D*)p23_ma103.Clone("ratio_23_4_ma103");
    TH1D* ratio_1_4_ma103=(TH1D*)p1.Clone("ratio_1_4_ma103");
    
    ratio_1_4_ma103->Divide(&p4_ma103);
    ratio_1_4_ma103->Draw("same");
    
    ratio_23_4_ma103->Divide(&p4_ma103);
    ratio_23_4_ma103->Draw("same");

    f_canvas->cd(3);
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetTickx();
    gPad->SetTicky();
    
    TH2D null_ma123("null_ma123",";"+xtitle+Form(";ratio to (NUANCE mode @ 1.23)"),1,lower,upper,1,0.5,1.5);
    null_ma123.SetStats(kFALSE);
    adjust_histo(null_ma123);
    null_ma123.Draw();
    
    TH1D* ratio_23_4_ma123=(TH1D*)p23_ma103.Clone("ratio_23_4_ma123");
    TH1D* ratio_1_4_ma123=(TH1D*)p1.Clone("ratio_1_4_ma123");
    
    ratio_1_4_ma123->Divide(&p4_ma123);
    ratio_1_4_ma123->Draw("same");
    
    ratio_23_4_ma123->Divide(&p4_ma123);
    ratio_23_4_ma123->Draw("same");
    
    f_canvas->Print(f_psFile);
  } //end loop over variables

}
////////////////////////////////////////////////////////////////////////
void MA_looker::make_benhar_and_nuance_comparisons() {
  bool onePerPage=true;

  BenharDrawer bed;
  ButkevichDrawer bud;

  TCanvas canvas("canvas","",500,700);
  if (onePerPage) {
    canvas.Divide(1,3);
    f_psOptions="Portrait";
  }
  
  TString psFile="benhar_compare.ps";
  canvas.Print(psFile+"[",f_psOptions);

  TString tree_var="-q[0]**2+q[1]**2+q[2]**2+q[3]**2";
  TString label=";Q^{2} (GeV^{2});d#sigma/dQ^{2} (cm^{2}/GeV^{2})";
  Int_t N_bins=70;
  Double_t lower=0.0;
  Double_t upper=2.0;

  //NUANCE 1.03 FP=0
  TH1D n_nfp_ma103("n_nfp_ma103",label,N_bins,lower,upper);
  n_nfp_ma103.Sumw2();
  n_nfp_ma103.SetStats(kFALSE);
  n_nfp_ma103.SetLineColor(kCyan);
  n_nfp_ma103.SetMarkerColor(kCyan);
  //NUANCE 1.20 FP=0
  TH1D n_nfp_ma120("n_nfp_ma120",label,N_bins,lower,upper);
  n_nfp_ma120.Sumw2();
  n_nfp_ma120.SetStats(kFALSE);
  n_nfp_ma120.SetLineColor(kCyan+3);
  n_nfp_ma120.SetMarkerColor(kCyan+3);
  //NUANCE 1.03 dipole FP
  TH1D n_fp_ma103("n_fp_ma103",label,N_bins,lower,upper);
  n_fp_ma103.Sumw2();
  n_fp_ma103.SetStats(kFALSE);
  n_fp_ma103.SetLineColor(kGreen);
  n_fp_ma103.SetMarkerColor(kGreen);

  {
    int N_files=40;
    for (Int_t iFile=0;iFile<N_files;iFile++) {
      TString filename_nfp=Form("../../nuance_from_fnal/numu_oxygen_dipole_1200/nuance%d.root",iFile+1);
      TString filename_fp=Form("../../nuance_from_fnal/numu_oxygen_dipole_1200_fp/nuance%d.root",iFile+1);

      TFile file_nfp(filename_nfp);
      gROOT->cd();
      TFile file_fp(filename_fp);
      gROOT->cd();
      TTree *tree_nfp=(TTree*)file_nfp.Get("h3");
      TTree *tree_fp=(TTree*)file_fp.Get("h3");
      tree_fp->Draw("-qsq/1.0e6>>+n_fp_ma103","cc && bound && channel==1","goff");
      tree_nfp->Draw("-qsq/1.0e6>>+n_nfp_ma103","cc && bound && channel==1","goff");
      file_fp.Close();
      file_nfp.Close();
    }
    n_nfp_ma103.Scale(0.0100135*8*1.0e-36/n_nfp_ma103.Integral("width"));
    n_fp_ma103.Scale(0.0100417*8*1.0e-36/n_fp_ma103.Integral("width"));
  }
  {
    int N_files=40;
    for (Int_t iFile=0;iFile<N_files;iFile++) {
      TString filename_nfp=Form("../../nuance_from_fnal/numu_oxygen_dipole_1200_nofp_ma120/nuance%d.root",iFile+1);

      TFile file_nfp(filename_nfp);
      gROOT->cd();
      TTree *tree_nfp=(TTree*)file_nfp.Get("h3");
      tree_nfp->Draw("-qsq/1.0e6>>+n_nfp_ma120","cc && bound && channel==1","goff");
      file_nfp.Close();
    }
    n_nfp_ma120.Scale(0.0117485*8*1.0e-36/n_nfp_ma120.Integral("width"));
  }

  //NUANCE mode (FP!=0)
  TFile f4_fp_ma103(Form("~ted/work/boone/ccqe/events/may09_events_2M_proc4_ma103_fp_mu_%d.root",f_energyPoint));
  gROOT->cd();
  TTree *tree4_fp_ma103=(TTree*)f4_fp_ma103.Get("tree");
  tree4_fp_ma103->SetName("tree4_fp_ma103");
  Double_t rate4_fp_ma103=0.0;
  {
    TGraph *gr=(TGraph*)f4_fp_ma103.Get("ccqe_rate")->Clone();
    Double_t dummy_x,dummy_y;
    gr->GetPoint(4,dummy_x,rate4_fp_ma103);
    delete gr;
  }

  //NUANCE mode (FP==0)
  TFile f4_nfp_ma103(Form("~ted/work/boone/ccqe/events/may09_events_2M_proc4_ma103_nfp_mu_%d.root",f_energyPoint));
  gROOT->cd();
  TTree *tree4_nfp_ma103=(TTree*)f4_nfp_ma103.Get("tree");
  tree4_nfp_ma103->SetName("tree4_nfp_ma103");
  Double_t rate4_nfp_ma103=0.0;
  {
    TGraph *gr=(TGraph*)f4_nfp_ma103.Get("ccqe_rate")->Clone();
    Double_t dummy_x,dummy_y;
    gr->GetPoint(4,dummy_x,rate4_nfp_ma103);
    delete gr;
  }

  TFile f4_nfp_ma120(Form("~ted/work/boone/ccqe/events/may09_events_2M_proc4_ma120_nfp_mu_%d.root",f_energyPoint));
  gROOT->cd();
  TTree *tree4_nfp_ma120=(TTree*)f4_nfp_ma120.Get("tree");
  tree4_nfp_ma120->SetName("tree4_nfp_ma120");
  Double_t rate4_nfp_ma120=0.0;
  {
    TGraph *gr=(TGraph*)f4_nfp_ma120.Get("ccqe_rate")->Clone();
    Double_t dummy_x,dummy_y;
    gr->GetPoint(4,dummy_x,rate4_nfp_ma120);
    delete gr;
  }

  //A-S SF mode
  TFile f123_ma100(Form("~ted/work/boone/ccqe/events/may09_events_2M_proc123_ma100_fp_mu_%d.root",f_energyPoint));
  gROOT->cd();
  TTree *tree123_ma100=(TTree*)f123_ma100.Get("tree");
  tree123_ma100->SetName("tree123_ma100");
  Double_t rate1_ma100=0.0;
  Double_t rate23_ma100=0.0;
  {
    TGraph *gr=(TGraph*)f123_ma100.Get("ccqe_rate")->Clone();
    Double_t dummy_x,dummy_y;
    gr->GetPoint(1,dummy_x,rate1_ma100);
    gr->GetPoint(2,dummy_x,dummy_y);
    rate23_ma100+=dummy_y;
    gr->GetPoint(3,dummy_x,dummy_y);
    rate23_ma100+=dummy_y;
    delete gr;
  }

  TFile f123_ma120(Form("~ted/work/boone/ccqe/events/may09_events_2M_proc123_ma120_fp_mu_%d.root",f_energyPoint));
  gROOT->cd();
  TTree *tree123_ma120=(TTree*)f123_ma120.Get("tree");
  tree123_ma120->SetName("tree123_ma120");
  Double_t rate1_ma120=0.0;
  Double_t rate23_ma120=0.0;
  {
    TGraph *gr=(TGraph*)f123_ma120.Get("ccqe_rate")->Clone();
    Double_t dummy_x,dummy_y;
    gr->GetPoint(1,dummy_x,rate1_ma120);
    gr->GetPoint(2,dummy_x,dummy_y);
    rate23_ma120+=dummy_y;
    gr->GetPoint(3,dummy_x,dummy_y);
    rate23_ma120+=dummy_y;
    delete gr;
  }

  TH1D p4_fp_ma103("p4_fp_ma103",label,N_bins,lower,upper);
  p4_fp_ma103.Sumw2();
  p4_fp_ma103.SetStats(kFALSE);
  tree4_fp_ma103->Draw(tree_var+">>p4_fp_ma103","(process==4)","goff",f_nEvents,1);
  p4_fp_ma103.Scale(rate4_fp_ma103/p4_fp_ma103.Integral("width"));
  p4_fp_ma103.SetLineColor(kYellow);
  p4_fp_ma103.SetMarkerColor(kYellow);

  TH1D p4_nfp_ma103("p4_nfp_ma103",label,N_bins,lower,upper);
  p4_nfp_ma103.Sumw2();
  p4_nfp_ma103.SetStats(kFALSE);
  tree4_nfp_ma103->Draw(tree_var+">>p4_nfp_ma103","(process==4)","goff",f_nEvents,1);
  p4_nfp_ma103.Scale(rate4_nfp_ma103/p4_nfp_ma103.Integral("width"));
  p4_nfp_ma103.SetLineColor(kGray);
  p4_nfp_ma103.SetMarkerColor(kGray);

  TH1D p4_nfp_ma120("p4_nfp_ma120",label,N_bins,lower,upper);
  p4_nfp_ma120.Sumw2();
  p4_nfp_ma120.SetStats(kFALSE);
  tree4_nfp_ma120->Draw(tree_var+">>p4_nfp_ma120","(process==4)","goff",f_nEvents,1);
  p4_nfp_ma120.Scale(rate4_nfp_ma120/p4_nfp_ma120.Integral("width"));
  p4_nfp_ma120.SetLineColor(kBlue-8);
  p4_nfp_ma120.SetMarkerColor(kBlue-8);

  TH1D p23_ma100("p23_ma100",label,N_bins,lower,upper);
  p23_ma100.Sumw2();
  p23_ma100.SetStats(kFALSE);
  tree123_ma100->Draw(tree_var+">>p23_ma100","(process==2 || process==3)","goff",f_nEvents,1);
  p23_ma100.Scale(rate23_ma100/p23_ma100.Integral("width"));
  p23_ma100.SetLineColor(kBlue);
  p23_ma100.SetMarkerColor(kBlue);

  TH1D p23_ma120("p23_ma120",label,N_bins,lower,upper);
  p23_ma120.Sumw2();
  p23_ma120.SetStats(kFALSE);
  tree123_ma120->Draw(tree_var+">>p23_ma120","(process==2 || process==3)","goff",f_nEvents,1);
  p23_ma120.Scale(rate23_ma120/p23_ma120.Integral("width"));
  p23_ma120.SetLineColor(kMagenta);
  p23_ma120.SetMarkerColor(kMagenta);

  TH1D p1_ma100("p1_ma100",label,N_bins,lower,upper);
  p1_ma100.Sumw2();
  p1_ma100.SetStats(kFALSE);
  tree123_ma100->Draw(tree_var+">>p1_ma100","process==1","goff",f_nEvents,1);
  p1_ma100.Scale(rate1_ma100/p1_ma100.Integral("width"));
  p1_ma100.SetLineColor(kRed);
  p1_ma100.SetMarkerColor(kRed);

  TH1D p1_ma120("p1_ma120",label,N_bins,lower,upper);
  p1_ma120.Sumw2();
  p1_ma120.SetStats(kFALSE);
  tree123_ma120->Draw(tree_var+">>p1_ma120","process==1","goff",f_nEvents,1);
  p1_ma120.Scale(rate1_ma120/p1_ma120.Integral("width"));
  p1_ma120.SetLineColor(kBlack);
  p1_ma120.SetMarkerColor(kBlack);

  TLegend legend(0.4,0.53,0.9,0.88);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  TString headerString="dipole form factors; O^{16}; ";
  if (f_energyPoint==0) headerString+=Form("#nu_{#mu}; MB flux",f_energyPoint);
  else headerString+=Form("%d MeV #nu_{#mu}",f_energyPoint);
  legend.SetHeader(headerString);
  legend.AddEntry(&p1_ma100  ,"A.-S. SM mode, M_{A}=1.00 GeV");
  legend.AddEntry(&p1_ma120  ,"A.-S. SM mode, M_{A}=1.20 GeV");
  legend.AddEntry(&p23_ma100 ,"A.-S. SF mode, M_{A}=1.00 GeV");
  legend.AddEntry(&p23_ma120 ,"A.-S. SF mode, M_{A}=1.20 GeV");
  //legend.AddEntry(&n_fp_ma103,"NUANCE, M_{A}=1.03 GeV");
  legend.AddEntry(&n_nfp_ma103,"NUANCE (FP=0), M_{A}=1.03 GeV");
  //legend.AddEntry(&p4_fp_ma103,"NUANCE mode, M_{A}=1.03 GeV");
  legend.AddEntry(&p4_nfp_ma103,"NUANCE mode (FP=0), M_{A}=1.03 GeV");

  TH1D *null=(TH1D*)p1_ma100.Clone("null");
  null->SetStats(false);
  null->GetXaxis()->CenterTitle();
  null->GetYaxis()->CenterTitle();
  null->GetYaxis()->SetTitleOffset(1.3);
  null->Reset();
  null->SetMaximum(0.14e-36);
  
  canvas.cd(1);
  null->Draw();

  std::vector<TH1D> histograms;
  std::vector<TGraph> be_grs;
  std::vector<TGraph> bu_grs;
  std::vector<TString> models;
  std::vector<double> mas;

  histograms.push_back(p1_ma100);  be_grs.push_back(bed.graph_fg_ma10_); bu_grs.push_back(bud.graph_fg_ma10_); models.push_back("SM"); mas.push_back(1.00);
  histograms.push_back(p1_ma120);  be_grs.push_back(bed.graph_fg_ma12_); bu_grs.push_back(bud.graph_fg_ma12_); models.push_back("SM"); mas.push_back(1.20);
  histograms.push_back(p23_ma100); be_grs.push_back(bed.graph_sf_ma10_); bu_grs.push_back(bud.graph_sf_ma10_); models.push_back("SF"); mas.push_back(1.00);
  histograms.push_back(p23_ma120); be_grs.push_back(bed.graph_sf_ma12_); bu_grs.push_back(bud.graph_sf_ma12_); models.push_back("SF"); mas.push_back(1.20);

  unsigned int nStop=4;

  for (unsigned int iElement=0;iElement<nStop;iElement++) {
    bool yes_butkevich=bu_grs.at(iElement).GetN()>0;

    histograms.at(iElement).Draw("esame");
    if (onePerPage) {
      be_grs.at(iElement).Draw("lsame");
      if (yes_butkevich) {
	bu_grs.at(iElement).SetMarkerStyle(20);	
	bu_grs.at(iElement).Draw("psame");
      }
      TLegend legend_local(0.4,0.53,0.9,0.88);
      legend_local.SetBorderSize(0);
      legend_local.SetFillStyle(0);
      legend_local.SetHeader(headerString);
      legend_local.AddEntry( &(histograms.at(iElement)) ,"A.-S. " +models.at(iElement)+ " mode, M_{A}="+Form("%4.3g GeV",mas.at(iElement)),"lpe");
      legend_local.AddEntry( &(be_grs.at(iElement))     ,"Benhar "+models.at(iElement)+" curve, M_{A}="+Form("%4.3g GeV",mas.at(iElement)),"l");
      if (yes_butkevich) {
	legend_local.AddEntry( &(bu_grs.at(iElement))     ,"Butkevich "+models.at(iElement)+" curve, M_{A}="+Form("%4.3g GeV",mas.at(iElement)),"p");
      }
      legend_local.Draw("same");

      canvas.cd(2);
      TH1D *ratio=(TH1D*)histograms.at(iElement).Clone("ratio");
      ratio->GetYaxis()->SetTitle("A.-S. / Benhar");
      make_ratio(ratio,be_grs.at(iElement));
      ratio->Draw();
      TLine line;
      line.DrawLine(0.0,1.0,2.0,1.0);

      TH1D *ratio2=0;
      if (yes_butkevich) {
	canvas.cd(3);
	ratio2=(TH1D*)histograms.at(iElement).Clone("ratio2");
	ratio2->GetYaxis()->SetTitle("A.-S. / Butkevich");
	make_ratio(ratio2,bu_grs.at(iElement));
	ratio2->Draw();
	line.DrawLine(0.0,1.0,2.0,1.0);
      }

      canvas.cd(0);
      canvas.Print(psFile);

      delete ratio;
      if (ratio2) delete ratio2;

      canvas.cd(1);
      null->Draw();
    }
    else if (iElement==nStop-1) {
      //n_fp_ma103.Draw("esame");
      n_nfp_ma103.Draw("esame");
      //p4_fp_ma103.Draw("esame");
      p4_nfp_ma103.Draw("esame");
      legend.Draw("same");
      canvas.Print(psFile);
    }
  }

  if (onePerPage) {
    n_nfp_ma103.Draw("esame");
    p4_nfp_ma103.Draw("esame");

    TLegend legend_local(0.4,0.53,0.9,0.88);
    legend_local.SetBorderSize(0);
    legend_local.SetFillStyle(0);
    legend_local.SetHeader(headerString);
    legend_local.AddEntry(&n_nfp_ma103,"NUANCE (FP=0), M_{A}=1.03 GeV","lpe");
    legend_local.AddEntry(&p4_nfp_ma103,"NUANCE mode (FP=0), M_{A}=1.03 GeV","lpe");
    legend_local.Draw("same");

    canvas.cd(2);
    TH1D *ratio=(TH1D*)n_nfp_ma103.Clone("ratio");
    ratio->Divide(&p4_nfp_ma103);
    ratio->GetXaxis()->CenterTitle();
    ratio->GetYaxis()->CenterTitle();
    ratio->GetYaxis()->SetTitle("NUANCE / NUANCE-mode");
    ratio->SetMinimum(0.60);
    ratio->SetMaximum(1.40);
    ratio->Draw();
    TLine line;
    line.DrawLine(0.0,1.0,2.0,1.0);
    canvas.cd(0);
    canvas.Print(psFile);
    delete ratio;
  }
  if (onePerPage) {
    canvas.cd(1);
    null->Draw();
    n_nfp_ma120.Draw("esame");
    p4_nfp_ma120.Draw("esame");
    bud.graph_fg_ma12_.SetMarkerStyle(20);
    bud.graph_fg_ma12_.Draw("psame");

    TLegend legend_local(0.4,0.53,0.9,0.88);
    legend_local.SetBorderSize(0);
    legend_local.SetFillStyle(0);
    legend_local.SetHeader(headerString);
    legend_local.AddEntry(&n_nfp_ma120,"NUANCE (FP=0), M_{A}=1.20 GeV","lpe");
    legend_local.AddEntry(&p4_nfp_ma120,"NUANCE mode (FP=0), M_{A}=1.20 GeV","lpe");
    legend_local.AddEntry(&(bud.graph_fg_ma12_),"Butkevich SM curve, M_{A}=1.20 GeV","p");
    legend_local.Draw("same");

    canvas.cd(2);
    TH1D *ratio=(TH1D*)n_nfp_ma120.Clone("ratio");
    ratio->Divide(&p4_nfp_ma120);
    ratio->GetXaxis()->CenterTitle();
    ratio->GetYaxis()->CenterTitle();
    ratio->GetYaxis()->SetTitle("NUANCE / NUANCE-mode");
    ratio->SetMinimum(0.60);
    ratio->SetMaximum(1.40);
    ratio->Draw();
    TLine line;
    line.DrawLine(0.0,1.0,2.0,1.0);

    canvas.cd(3);
    TH1D *ratio2=(TH1D*)n_nfp_ma120.Clone("ratio2");
    make_ratio(ratio2,bud.graph_fg_ma12_);
    ratio2->GetXaxis()->CenterTitle();
    ratio2->GetYaxis()->CenterTitle();
    ratio2->GetYaxis()->SetTitle("NUANCE / Butkevich");
    ratio2->SetMinimum(0.60);
    ratio2->SetMaximum(1.40);
    ratio2->Draw();
    line.DrawLine(0.0,1.0,2.0,1.0);

    canvas.cd(0);
    canvas.Print(psFile);
    delete ratio;
    delete ratio2;
  }

  canvas.Print(psFile+"]",f_psOptions);
  delete null;
}
////////////////////////////////////////////////////////////////////////
void MA_looker::make_ratio(TH1D *ratio,TGraph& graph) {
  TH1D *denominator=(TH1D*)ratio->Clone("denominator");
  double xmin,xmax,ymin,ymax;
  graph.ComputeRange(xmin,ymin,xmax,ymax);
  for (int iBin=1;iBin<=denominator->GetNbinsX();iBin++) {
    double x=denominator->GetBinCenter(iBin);
    double y=graph.Eval(x);
    if (x>=xmin && x<=xmax) {
      denominator->SetBinContent(iBin,y);
    }
    else {
      ratio->SetBinContent(iBin,0.0);
      ratio->SetBinError(iBin,0.0);
    }
  }
  ratio->Divide(denominator);
  ratio->GetXaxis()->CenterTitle();
  ratio->GetYaxis()->CenterTitle();
  ratio->SetMinimum(0.60);
  ratio->SetMaximum(1.40);

  delete denominator;
}
////////////////////////////////////////////////////////////////////////
void MA_looker::makePsFile() {
  f_canvas->Print(f_psFile+"[",f_psOptions);
  compare_simple();
  f_canvas->Print(f_psFile+"]",f_psOptions);
}
////////////////////////////////////////////////////////////////////////
int main() {
  //MA_looker ml(0);
  //ml.makePsFile();

  MA_looker ml(1200);
  ml.make_benhar_and_nuance_comparisons();
  return 0;
}
////////////////////////////////////////////////////////////////////////
void MA_looker::XS_look() {

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


  g1.SetMarkerStyle(20);g1.SetMarkerColor(kBlack);
  g2.SetMarkerStyle(20);g2.SetMarkerColor(kRed);
  g3.SetMarkerStyle(20);g3.SetMarkerColor(kBlue);
  g4.SetMarkerStyle(20);g4.SetMarkerColor(kGreen);

  g1.Draw("ap");
  g2.Draw("psame");
  g3.Draw("psame");
  g4.Draw("psame");
  f_canvas->Print(f_psFile);
}

////////////////////////////////////////////////////////////////////////
void MA_compare(Int_t mode,Int_t energy_point) {

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
