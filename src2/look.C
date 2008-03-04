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
void chi2_compare(TH1D *h1,TH1D *h2,Int_t bin_i,Int_t bin_f,TH1D* contribs) {
  Double_t chi2=0.0;
  Int_t ndf=0;
  for (Int_t iBin=bin_i;iBin<=bin_f;iBin++) {
    Double_t err1=h1->GetBinError(iBin);
    Double_t err2=h2->GetBinError(iBin);
    Double_t err=sqrt(err1*err1+err2*err2);
    if (err1>0.0 && err2>0.0) {
      Double_t ans=(h1->GetBinContent(iBin)-h2->GetBinContent(iBin))/err;
      if (contribs) contribs->SetBinContent(iBin,ans*ans);
      chi2+=ans*ans;
      ndf++;
    }
  }
  TString str=Form("#chi^{2}=%5.3g, ndf=%d, prob.=%5.3g",chi2,ndf,TMath::Prob(chi2,ndf));
  contribs->SetStats(kFALSE);
  contribs->GetYaxis()->SetTitle("#chi^{2} contribution per bin");
  contribs->SetTitle(str);
}

////////////////////////////////////////////////////////////////////////
void histo_graph(TH1D *histo,TGraph *graph,Int_t N_to_fill) {
  Int_t N_Points=graph->GetN();
  Double_t x,y;

  Double_t y_min=0.0;
  Double_t y_max=0.0;
  Double_t x_min=0.0;
  Double_t x_max=0.0;

  for (Int_t i=0;i<N_Points;i++) {
    graph->GetPoint(i,x,y);
    if (y>y_max) y_max=y;
    if (x>x_max) x_max=x;
    if (y<y_min) y_min=y;
    if (x<x_min) x_min=x;
  }
  Double_t histo_x_min=histo->GetXaxis()->GetXmin();
  Double_t histo_x_max=histo->GetXaxis()->GetXmax();
  if (histo_x_min>x_min) x_min=histo_x_min;
  if (histo_x_max<x_max) x_max=histo_x_max;

  //reverse graph order if necessary
  TGraph graph2(N_Points);
  Double_t x0,xNm1;
  graph->GetPoint(0,x0,y);
  graph->GetPoint(N_Points-1,xNm1,y);
  if (xNm1<x0) {
    for (Int_t iPoint=0;iPoint<N_Points;iPoint++) {
      Double_t x,y;
      graph->GetPoint(N_Points-iPoint,x,y);
      graph2.SetPoint(iPoint,x,y);
    }
    graph=&graph2;
  }

  TRandom3 r;
  r.SetSeed(0);

  Int_t N_tried=0;
  Int_t N_filled=0;
  while (N_filled<N_to_fill) {
    N_tried++;
    //rejection method
    x=r.Uniform(x_min,x_max);
    y=r.Uniform(0.0,y_max);
    
    if (y<graph->Eval(x)) {
      histo->Fill(x);
      N_filled++;
    }
  }
  printf("x_max=%8.6g; x_min=%8.6g; y_max=%8.6g; N_filled=%d; N_tried=%d\n",x_max,x_min,y_max,N_filled,N_tried);
  Double_t norm=(x_max-x_min)*y_max*N_filled/N_tried;
  histo->Scale(norm/histo->Integral("width"));

}

////////////////////////////////////////////////////////////////////////
void single_compare(Int_t mode) {
  Init();
  
  //LS
  TFile f_LS("../ref/katori/single_nofp_mm.root");
  gROOT->cd();
  TGraph *gr_LS=(TGraph*)f_LS.Get("xs_Q2")->Clone();
  f_LS.Close();

  gr_LS->SetTitle(";Q^{2} (GeV^{2});d#sigma/dQ^{2} (cm^{2}/GeV^{2})");
  gr_LS->GetXaxis()->CenterTitle();
  gr_LS->GetYaxis()->CenterTitle();
  gr_LS->GetYaxis()->SetTitleOffset(1.2);

  gr_LS->SetMarkerStyle(20);
  gr_LS->SetMarkerSize(0.3);


  //NUANCE
  Int_t N_files=20;
  TFile *files_d2_nofp[N_files];
  TTree *trees_d2_nofp[N_files];

  for (Int_t iFile=0;iFile<N_files;iFile++) {
    files_d2_nofp[iFile]=new TFile(Form("$CONDOR_TMP/nuance/d2/nuance%d.root",iFile+1));
    //files_d2_nofp[iFile]=new TFile(Form("../nuance_from_fnal/d2/nuance%d.root",iFile+1));
    gROOT->cd();
    trees_d2_nofp[iFile]=(TTree*)(files_d2_nofp[iFile]->Get("h3"));
  }


  //mine
  TFile f("../../../events/single_nofp_mm.root");
  gROOT->cd();
  TTree *tree=(TTree*)f.Get("tree");
  TGraph *gr=(TGraph*)f.Get("ccqe_rate")->Clone();
  Int_t gr_size=gr->GetN();
  Double_t *graph_proc = new Double_t[gr_size];
  Double_t *graph_xs   = new Double_t[gr_size];
  for (Int_t iProc=0;iProc<gr_size;iProc++) {
    gr->GetPoint(iProc,graph_proc[iProc],graph_xs[iProc]);
  }
  delete gr;


  //real stuff
  Int_t N_bins=100;
  Double_t lower=0.0;
  Double_t upper=1.2;
  TH1D h1("h1",";Q^{2} (GeV^{2});d#sigma/dQ^{2} (cm^{2}/GeV^{2})",N_bins,lower,upper);
  h1.Sumw2();
  TH1D h2("h2",";Q^{2} (GeV^{2});d#sigma/dQ^{2} (cm^{2}/GeV^{2})",N_bins,lower,upper);
  h2.Sumw2();

  h1.SetStats(kFALSE);
  h1.GetXaxis()->CenterTitle();
  h1.GetYaxis()->CenterTitle();
  h1.GetYaxis()->SetTitleOffset(1.2);
  h1.SetLineColor(kRed);
  h1.SetMarkerColor(kRed);

  if (mode==90) {
    tree->Draw("-q[0]**2+q[1]**2+q[2]**2+q[3]**2>>h1","process==0","goff");
    h1.Scale(graph_xs[0]/h1.Integral("width"));
    h1.DrawCopy("e");

    histo_graph(&h2,gr_LS,10000000);
    h2.DrawCopy("same");
    chi2_compare(&h1,&h2,1,2*N_bins/3,0);
    
    TLegend legend(0.4,0.5,0.9,0.9,"#nu_{#mu} + n #rightarrow #mu + p (E_{#nu}=0.8 GeV)");
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.AddEntry(&h1,"my generated events (reduction)","p");
    legend.AddEntry(&h2,"Llewellyn-Smith (BooNE memo)","p");
    legend.DrawClone("same");
  }
  if (mode==91) {
    gr_LS->DrawClone("ap");
    tree->SetMarkerColor(kRed);
    tree->Draw("xs/(0.93827*2.0):-q[0]**2+q[1]**2+q[2]**2+q[3]**2","process==0","same");
  }

  if (mode==92) {
    for (Int_t iTree=0;iTree<N_files;iTree++) {
      trees_d2_nofp[iTree]->Draw("-qsq/1.0e6>>+h1","cc && channel==1","goff");
    }
    h1.Scale(0.0104954*1.0e-36/h1.Integral("width"));
    h1.SetMarkerColor(kBlue);
    h1.DrawCopy("e");

    histo_graph(&h2,gr_LS,10000000);
    h2.DrawCopy("same");
    
    chi2_compare(&h1,&h2,1,2*N_bins/3,0);
  }

  delete gr_LS;
  delete graph_proc;
  delete graph_xs;
  for (Int_t iFile=0;iFile<N_files;iFile++) {
    delete files_d2_nofp[iFile];
  }

}


////////////////////////////////////////////////////////////////////////
void polish_compare(Int_t proc,Int_t mode,TString file) {
  Init();

  //mine
  TFile f(file);
  gROOT->cd();
  TTree *tree=(TTree*)f.Get("tree");
  TGraph *gr=(TGraph*)f.Get("ccqe_rate")->Clone();
  Int_t gr_size=gr->GetN();
  Double_t *graph_proc = new Double_t[gr_size];
  Double_t *graph_xs   = new Double_t[gr_size];
  for (Int_t iProc=0;iProc<gr_size;iProc++) {
    gr->GetPoint(iProc,graph_proc[iProc],graph_xs[iProc]);
  }
  delete gr;

  
  //polish
  TFile f_polish("../ref/polish/O16/graph.root");
  gROOT->cd();
  
  TGraph *gr_pol=0;
  TString *var=0;
  TString *title=0;

  TString *model=0;
  TString legend_string="";
  switch(proc) {
  case 1:
    model=new TString("sm");
    legend_string="Ankowski-Sobczyk Fermi Gas Model";
    break;
  case 2:
    model=new TString("as_mf");
    break;
  case 3:
    model=new TString("as_corr");
    break;
  case 23:
    model=new TString("as");
    legend_string="Ankowski-Sobczyk Spectral Function Model";
    break;
  }
  
  Int_t N_bins=0;
  Double_t lower=0.0;
  Double_t upper=0.0;
  Double_t x_leg=0.1;

  if (mode==50) {
    gr_pol=(TGraph*)f_polish.Get("polish_Q2_"+*model)->Clone();
    N_bins=100;
    lower=0.0;
    upper=1.4;
    var=new TString("q[1]**2+q[2]**2+q[3]**2-q[0]**2");
    title=new TString(";Q^{2} (GeV^{2});d#sigma/dQ^{2} (cm^{2}/GeV^{2})");
  }
  
  if (mode==60) {
    gr_pol=(TGraph*)f_polish.Get("polish_e_"+*model)->Clone();
    N_bins=100;
    lower=0.1;
    upper=0.8;
    x_leg=0.1;
    var=new TString("kprime[0]");
    title=new TString(";E_{#mu} (GeV);d#sigma/dE_{#mu} (cm^{2}/GeV)");
  }
  
  if (mode==70) {
    gr_pol=(TGraph*)f_polish.Get("polish_nopb_e_"+*model)->Clone();
    N_bins=100;
    lower=0.1;
    upper=0.8;
    var=new TString("kprime[0]");
    title=new TString(";E_{#mu} (GeV);d#sigma/dE_{#mu} (cm^{2}/GeV)");
  }
  
  f_polish.Close();

  TH1D h1("h1",*title,N_bins,lower,upper);
  h1.Sumw2();
  TH1D h2("h2",*title,N_bins,lower,upper);
  h2.Sumw2();
  TH1D *contribs=(TH1D*)h1.Clone("contribs");
  
  TString cut=Form("process==%d",proc);
  Double_t xs=graph_xs[proc];
  if (proc==23) {
    cut="process==2 || process==3";
    xs=graph_xs[2]+graph_xs[3];
  }
  tree->Draw(*var+">>h1",cut,"goff");
  h1.Scale(xs/h1.Integral("width"));
  
  printf("proc=%d; norm=%8.6g\n",proc,xs);
  
  TCanvas canvas("canvas","",700,1000);
  canvas.Divide(1,2);
  
  canvas.cd(1);
  h1.SetStats(kFALSE);
  h1.GetXaxis()->CenterTitle();
  h1.GetYaxis()->CenterTitle();
  h1.GetYaxis()->SetTitleOffset(1.2);
  h1.SetLineColor(kRed);
  h1.SetMarkerColor(kRed);
  h1.Draw("e");

  //Double_t previous_x=0.0;
  //for (Int_t iPoint=0;iPoint<gr_pol->GetN();iPoint++) {
  //  Double_t this_x,this_y;
  //  gr_pol->GetPoint(iPoint,this_x,this_y);
  //  cout << iPoint << "," << this_x-previous_x << endl;
  //  previous_x=this_x;
  //}
  gr_pol->Draw("psame");
  histo_graph(&h2,gr_pol,4000000);
  h2.Draw("same");
  
  TLegend legend(x_leg,0.60,x_leg+0.65,0.9,"#nu_{#mu} + ^{16}O #rightarrow #mu + p + X (E_{#nu}=0.8 GeV)");
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.AddEntry(&h1,"my generated events ("+legend_string+")","l");
  legend.AddEntry(&h2,"Ankowski-Sobczyk paper ("+legend_string+")","l");
  legend.Draw("same");
  
  chi2_compare(&h1,&h2,1,N_bins,contribs);
  
  canvas.cd(2);
  contribs->GetXaxis()->CenterTitle();
  contribs->GetYaxis()->CenterTitle();
  contribs->Draw();
  
  canvas.cd(0);
  canvas.DrawClone();
  
  if (gr_pol) delete gr_pol;
  if (model)  delete model;
  if (title)  delete title;
  if (var)    delete var;
  delete graph_proc;
  delete graph_xs;
}

////////////////////////////////////////////////////////////////////////
void NUANCE_compare(Int_t mode,Int_t energy_point,TString file) {
  Init();

  //mine
  TFile f(file);
  TTree *tree=(TTree*)f.Get("tree");
  TGraph *gr=(TGraph*)f.Get("ccqe_rate")->Clone();
  Int_t process=4;
  Double_t xs,dummy;
  gr->GetPoint(process,dummy,xs);
  delete gr;

  //NUANCE
  /*argh--chains are unusable in this ROOT version*/
  Int_t N_files=0;
  if (energy_point==300) N_files=40;
  if (energy_point==800) N_files=20;
  TFile *nu_files[N_files];
  TTree *nu_trees[N_files];
  for (Int_t iFile=0;iFile<N_files;iFile++) {
    TString file;
    if (energy_point==300) file=Form("../../nuance_from_fnal/oxnofp3/nuance%d.root",iFile+1);
    if (energy_point==800) file=Form("../../nuance_from_fnal/oxyg_no_fp2/nuance%d.root",iFile+1);
    nu_files[iFile]=new TFile(file);
    //nu_files[iFile]=new TFile(Form("$CONDOR_TMP/nuance/oxyg_no_fp2_300/nuance%d.root",iFile+1));
    //nu_files[iFile]=new TFile(Form("$CONDOR_TMP/nuance/oxyg_no_fp2/nuance%d.root",iFile+1));
    //nu_files[iFile]=new TFile(Form("$CONDOR_TMP/nuance/oxyg_no_fp_nor_eb/nuance%d.root",iFile+1));
    gROOT->cd();
    nu_trees[iFile]=(TTree*)(nu_files[iFile]->Get("h3"));
  }

  TString *tree_var=0;
  TString *nuance_var=0;
  TCut *nuance_cut=0;
  TString *label=0;

  Int_t N_bins=0;
  Double_t lower=0.0;
  Double_t upper=0.0;
  Double_t x_leg=0.0;

  switch (mode) {
  case 1:
    tree_var=new TString("kprime[0]");
    nuance_var=new TString("p_lepton[0][3]/1.0e3");
    nuance_cut=new TCut("cc && bound && channel==1");
    label=new TString(";E_{#mu} (GeV);d#sigma/dE_{#mu} (cm^{2}/GeV)");
    
    x_leg=0.15;
    N_bins=250;
    lower=0.12;
    upper=energy_point/1.0e3;
    break;
  case 2:
    tree_var=new TString("-q[0]**2+q[1]**2+q[2]**2+q[3]**2");
    nuance_var=new TString("-qsq/1.0e6");
    nuance_cut=new TCut("cc && bound && channel==1");
    label=new TString(";Q^{2} (GeV^{2});d#sigma/dQ^{2} (cm^{2}/GeV^{2})");

    x_leg=0.55;
    N_bins=250;
    lower=0.0;
    if (energy_point==300) upper=0.8;
    if (energy_point==800) upper=1.4;
    break;
  case 3:
    tree_var=new TString("mag_p");
    nuance_var=new TString("p_targ[4]/1.0e3");
    nuance_cut=new TCut("cc && bound && channel==1");
    label=new TString(";p_{i} (GeV/c);d#sigma/dp_{i} (cm^{2}/(GeV/c))");

    x_leg=0.15;
    N_bins=150;
    lower=0.0;
    upper=0.25;
    break;
  case 4:
    tree_var=new TString("sqrt(pprime[1]**2+pprime[2]**2+pprime[3]**2)");
    nuance_var=new TString("p_hadron[][4]/1.0e3");
    nuance_cut=new TCut("cc && bound && channel==1 && n_hadrons==1 && n_leptons==1");
    label=new TString(";p_{f} (GeV/c);d#sigma/dp_{f} (cm^{2}/(GeV/c))");

    x_leg=0.15;
    N_bins=150;
    lower=0.0;
    upper=0.8;
    break;
  case 5:
    tree_var=new TString("sqrt(pprime[1]**2+pprime[2]**2+pprime[3]**2)");
    nuance_var=new TString("sqrt((p_neutrino[0]-p_lepton[0][0]+p_targ[0])**2+(p_neutrino[1]-p_lepton[0][1]+p_targ[1])**2+(p_neutrino[2]-p_lepton[0][2]+p_targ[2])**2)/1.0e3");
    nuance_cut=new TCut("cc && bound && channel==1 && n_leptons==1");
    label=new TString(";p_{f} (GeV/c);d#sigma/dp_{f} (cm^{2}/(GeV/c))");

    x_leg=0.15;
    N_bins=150;
    lower=0.2;
    upper=0.6;
    //lower=0.2;
    //upper=1.3;
    break;
  }

  TCanvas canvas("canvas","",700,1000);
  canvas.Divide(1,2);
  canvas.cd(1);

  TH1D nh("nh",*label,N_bins,lower,upper);
  nh.Sumw2();
  nh.SetStats(kFALSE);
  for (Int_t iFile=0;iFile<N_files;iFile++) {
    nu_trees[iFile]->Draw(*nuance_var+">>+nh",*nuance_cut,"goff");
  }
  if (energy_point==300) nh.Scale(0.0033647*8*1.0e-36/nh.Integral("width"));
  if (energy_point==800) nh.Scale(0.0093045*8*1.0e-36/nh.Integral("width"));
  nh.Draw("e");

  TH1D *contribs=new TH1D("contribs",*label,N_bins,lower,upper);
  
  Int_t NHistos=1;
  Int_t base_index=0;
  TH1D *h[NHistos];
  for (Int_t iHisto=0;iHisto<NHistos;iHisto++) {
    h[iHisto]=new TH1D(Form("h%d",iHisto),*label,N_bins,lower,upper);
    h[iHisto]->Sumw2();
  }
  
  Int_t color=2;
  Double_t factor=1.0;
  for (Int_t iHisto=0;iHisto<NHistos;iHisto++) {
    tree->Draw(*tree_var+Form(">>h%d",iHisto),Form("(MA_weights[%d]) * (process==%d)",iHisto,process),"goff");
    if (iHisto==0) factor=xs/h[0]->Integral("width");
    h[iHisto]->Scale(factor);
    h[iHisto]->SetStats(kFALSE);
    h[iHisto]->GetXaxis()->CenterTitle();
    h[iHisto]->GetYaxis()->CenterTitle();
    h[iHisto]->GetYaxis()->SetTitleOffset(1.2);
    h[iHisto]->SetLineColor(color);
    h[iHisto]->SetMarkerColor(color);
    if (++color==5) color++;
    h[iHisto]->Draw("esame"); 
    h[iHisto]->Print();
  }

  TLegend legend(x_leg,0.7,x_leg+0.45,0.9,Form("#nu_{#mu} + ^{16}O #rightarrow #mu + p + X (E_{#nu}=%d MeV)",energy_point));
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  //legend.SetHeader("S-M, dipole, FP=0");
  legend.AddEntry(h[base_index],"my generated events (NUANCE mode)");
  legend.AddEntry(&nh,"NUANCE");
  legend.Draw("same");
  
  chi2_compare(h[base_index],&nh,1,N_bins,contribs);
  
  //TLatex text(0.20,0.92,comp_str);
  //text.SetNDC();
  //text.DrawClone();
  
  canvas.cd(2);
  contribs->Draw();

  canvas.cd(0);
  canvas.DrawClone();

  f.Close();

  delete tree_var;
  delete nuance_var;
  delete nuance_cut;
  delete label;
  delete contribs;
  
  for (Int_t iHisto=0;iHisto<NHistos;iHisto++) {
    delete h[iHisto];
  }

  for (Int_t iFile=0;iFile<N_files;iFile++) {
    delete nu_files[iFile];
  }


}

////////////////////////////////////////////////////////////////////////
void flux_compare(Int_t mode,TString linlog) {
  Init();

  //mine
  TFile f("~/work/boone/ccqe/xs_ted/../events/events_mb_4.root");
  TTree *tree=(TTree*)f.Get("tree");
  TGraph *gr=(TGraph*)f.Get("ccqe_rate")->Clone();
  Int_t process=4;
  Double_t rate,dummy;
  gr->GetPoint(process,dummy,rate);
  delete gr;
  Int_t NFiles=40;
  TFile *nu_files[NFiles];
  TTree *nu_trees[NFiles];
  for (Int_t iFile=0;iFile<NFiles;iFile++) {
    TString file;
    file="../../nuance_from_fnal/oxygen_dipole_mbflux_e"+linlog+Form("/nuance%d.root",iFile+1);
    nu_files[iFile]=new TFile(file);
    gROOT->cd();
    nu_trees[iFile]=(TTree*)(nu_files[iFile]->Get("h3"));
  }

  const Int_t NLines=29;
  Double_t linex[NLines];
  if (linlog=="log") {
    Int_t iLine;
    iLine= 0; linex[iLine]=0.1371614;
    iLine= 1; linex[iLine]=0.1480541;
    iLine= 2; linex[iLine]=0.1725035;
    iLine= 3; linex[iLine]=0.2009903;
    iLine= 4; linex[iLine]=0.2341815;
    iLine= 5; linex[iLine]=0.2728536;
    iLine= 6; linex[iLine]=0.3179121;
    iLine= 7; linex[iLine]=0.3704115;
    iLine= 8; linex[iLine]=0.4315804;
    iLine= 9; linex[iLine]=0.5028507;
    iLine=10; linex[iLine]=0.5858902;
    iLine=11; linex[iLine]=0.6826431;
    iLine=12; linex[iLine]=0.7953730;
    iLine=13; linex[iLine]=0.9267194;
    iLine=14; linex[iLine]=1.0797561;
    iLine=15; linex[iLine]=1.2580642;
    iLine=16; linex[iLine]=1.4658186;
    iLine=17; linex[iLine]=1.7078801;
    iLine=18; linex[iLine]=1.9899162;
    iLine=19; linex[iLine]=2.3185259;
    iLine=20; linex[iLine]=2.7014029;
    iLine=21; linex[iLine]=3.1475057;
    iLine=22; linex[iLine]=3.6672788;
    iLine=23; linex[iLine]=4.2728839;
    iLine=24; linex[iLine]=4.9784997;
    iLine=25; linex[iLine]=5.8006397;
    iLine=26; linex[iLine]=6.7585426;
    iLine=27; linex[iLine]=7.8746354;
    iLine=28; linex[iLine]=8.5000022;
  }
  if (linlog=="lin") {
    Int_t iLine;
    iLine= 0; linex[iLine]=0.1371614;
    iLine= 1; linex[iLine]=0.2920288;
    iLine= 2; linex[iLine]=0.6017636;
    iLine= 3; linex[iLine]=0.9114983;
    iLine= 4; linex[iLine]=1.2212331;
    iLine= 5; linex[iLine]=1.5309679;
    iLine= 6; linex[iLine]=1.8407026;
    iLine= 7; linex[iLine]=2.1504375;
    iLine= 8; linex[iLine]=2.4601721;
    iLine= 9; linex[iLine]=2.7699069;
    iLine=10; linex[iLine]=3.0796416;
    iLine=11; linex[iLine]=3.3893764;
    iLine=12; linex[iLine]=3.6991110;
    iLine=13; linex[iLine]=4.0088459;
    iLine=14; linex[iLine]=4.3185805;
    iLine=15; linex[iLine]=4.6283154;
    iLine=16; linex[iLine]=4.9380502;
    iLine=17; linex[iLine]=5.2477851;
    iLine=18; linex[iLine]=5.5575195;
    iLine=19; linex[iLine]=5.8672543;
    iLine=20; linex[iLine]=6.1769892;
    iLine=21; linex[iLine]=6.4867241;
    iLine=22; linex[iLine]=6.7964589;
    iLine=23; linex[iLine]=7.1061933;
    iLine=24; linex[iLine]=7.4159282;
    iLine=25; linex[iLine]=7.7256630;
    iLine=26; linex[iLine]=8.0353979;
    iLine=27; linex[iLine]=8.3451328;
    iLine=28; linex[iLine]=      8.5;
  }

  TString *tree_var=0;
  TString *nuance_var=0;
  TCut *nuance_cut=0;
  TString *label=0;

  Int_t N_bins=0;
  Double_t lower=0.0;
  Double_t upper=0.0;
  Double_t x_leg=0.0;

  switch (mode) {
  case 1:
    tree_var=new TString("kprime[0]");
    nuance_var=new TString("p_lepton[0][3]/1.0e3");
    nuance_cut=new TCut("cc && bound && channel==1");
    label=new TString(";E_{#mu} (GeV);d#sigma/dE_{#mu} (cm^{2}/GeV)");
    
    x_leg=0.15;
    N_bins=250;
    lower=0.12;
    upper=2.0;
    break;
  case 2:
    tree_var=new TString("-q[0]**2+q[1]**2+q[2]**2+q[3]**2");
    nuance_var=new TString("-qsq/1.0e6");
    nuance_cut=new TCut("cc && bound && channel==1");
    label=new TString(";Q^{2} (GeV^{2});d#sigma/dQ^{2} (cm^{2}/GeV^{2})");

    x_leg=0.55;
    N_bins=250;
    lower=0.0;
    upper=2.0;
    break;
  case 3:
    tree_var=new TString("mag_p");
    nuance_var=new TString("p_targ[4]/1.0e3");
    nuance_cut=new TCut("cc && bound && channel==1");
    label=new TString(";p_{i} (GeV/c);d#sigma/dp_{i} (cm^{2}/(GeV/c))");

    x_leg=0.15;
    N_bins=150;
    lower=0.0;
    upper=0.25;
    break;
  case 4:
    tree_var=new TString("sqrt(pprime[1]**2+pprime[2]**2+pprime[3]**2)");
    nuance_var=new TString("p_hadron[][4]/1.0e3");
    nuance_cut=new TCut("cc && bound && channel==1 && n_hadrons==1 && n_leptons==1");
    label=new TString(";p_{f} (GeV/c);d#sigma/dp_{f} (cm^{2}/(GeV/c))");

    x_leg=0.15;
    N_bins=150;
    lower=0.0;
    upper=0.8;
    break;
  case 5:
    tree_var=new TString("sqrt(pprime[1]**2+pprime[2]**2+pprime[3]**2)");
    nuance_var=new TString("sqrt((p_neutrino[0]-p_lepton[0][0]+p_targ[0])**2+(p_neutrino[1]-p_lepton[0][1]+p_targ[1])**2+(p_neutrino[2]-p_lepton[0][2]+p_targ[2])**2)/1.0e3");
    nuance_cut=new TCut("cc && bound && channel==1 && n_leptons==1");
    label=new TString(";p_{f} (GeV/c);d#sigma/dp_{f} (cm^{2}/(GeV/c))");

    x_leg=0.15;
    N_bins=150;
    lower=0.2;
    upper=1.3;
    break;
  case 6:
    tree_var=new TString("k[0]");
    nuance_var=new TString("p_neutrino[3]/1.0e3");
    nuance_cut=new TCut("cc && bound && channel==1 && n_leptons==1");
    label=new TString(";E_{#nu} (GeV);d#sigma/dE_{#nu} (cm^{2}/GeV)");

    x_leg=0.15;
    N_bins=150;
    lower=0.0;
    upper=2.0;
    break;
  }

  TCanvas canvas("canvas","",700,1000);
  canvas.Divide(1,2);
  canvas.cd(1);

  TH1D nh("nh",*label,N_bins,lower,upper);
  nh.Sumw2();
  nh.SetStats(kFALSE);
  for (Int_t iFile=0;iFile<NFiles;iFile++) {
    nu_trees[iFile]->Draw(*nuance_var+">>+nh",*nuance_cut,"goff");
  }
  //for elin
  //nh.Scale(1.100e-09*8*1.0e-36/nh.Integral("width"));
  //for elog
  //nh.Scale(1.124e-09*8*1.0e-36/nh.Integral("width"));
  nh.Scale(1.0/nh.Integral("width"));
  nh.Draw("e");
  nh.Print();

  TH1D *contribs=new TH1D("contribs",*label,N_bins,lower,upper);
  
  Int_t NHistos=1;
  Int_t base_index=0;
  TH1D *h[NHistos];
  for (Int_t iHisto=0;iHisto<NHistos;iHisto++) {
    h[iHisto]=new TH1D(Form("h%d",iHisto),*label,N_bins,lower,upper);
    h[iHisto]->Sumw2();
  }
  
  Int_t color=2;
  Double_t factor=1.0;
  for (Int_t iHisto=0;iHisto<NHistos;iHisto++) {
    tree->Draw(*tree_var+Form(">>h%d",iHisto),Form("(MA_weights[%d]) * (process==%d)",iHisto,process),"goff");
    //if (iHisto==0) factor=rate/h[0]->Integral("width");
    if (iHisto==0) factor=1.0/h[0]->Integral("width");
    h[iHisto]->Scale(factor);
    h[iHisto]->SetStats(kFALSE);
    h[iHisto]->GetXaxis()->CenterTitle();
    h[iHisto]->GetYaxis()->CenterTitle();
    h[iHisto]->GetYaxis()->SetTitleOffset(1.2);
    h[iHisto]->SetLineColor(color);
    h[iHisto]->SetMarkerColor(color);
    if (++color==5) color++;
    h[iHisto]->Draw("esame"); 
    h[iHisto]->Print();
  }

  Double_t y1=nh.GetYaxis()->GetXmin();
  Double_t y2=nh.GetYaxis()->GetXmax();
  Double_t xmin=nh.GetXaxis()->GetXmin();
  Double_t xmax=nh.GetXaxis()->GetXmax();
  for (Int_t iLine=0;iLine<NLines;iLine++) {
    TLine line;
    line.SetLineColor(kGreen);
    Double_t x=linex[iLine];
    if (x<xmin || x>xmax) continue;
    line.DrawLine(x,y1,x,y2);
  }

  TLegend legend(x_leg,0.7,x_leg+0.45,0.9,"#nu_{#mu} + ^{16}O #rightarrow #mu + p + X (E_{#nu} from MB flux)");
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  //legend.SetHeader("S-M, dipole, FP=0");
  legend.AddEntry(h[base_index],"my generated events (NUANCE mode)");
  legend.AddEntry(&nh,"NUANCE");
  legend.Draw("same");
  
  chi2_compare(h[base_index],&nh,1,N_bins,contribs);
  
  //TLatex text(0.20,0.92,comp_str);
  //text.SetNDC();
  //text.DrawClone();
  
  canvas.cd(2);
  contribs->Draw();

  canvas.cd(0);
  canvas.DrawClone();

  f.Close();

  delete tree_var;
  delete nuance_var;
  delete nuance_cut;
  delete label;
  delete contribs;
  
  for (Int_t iHisto=0;iHisto<NHistos;iHisto++) {
    delete h[iHisto];
  }

}

////////////////////////////////////////////////////////////////////////
void MA_compare1(Int_t mode,Int_t energy_point) {
  Init();

  TFile g(Form("../../events_6_nofp_%3d.root",energy_point));
  gROOT->cd();
  TTree *gtree=(TTree*)g.Get("tree");
  gtree->SetName("gtree");
  TGraph *ggr=(TGraph*)g.Get("ccqe_rate")->Clone();
  Int_t gprocess=1;
  Double_t gxs,gdummy;
  ggr->GetPoint(gprocess,gdummy,gxs);
  delete ggr;

  //mine
  TFile f(Form("../../events_5_nofp_%3d.root",energy_point));
  gROOT->cd();
  TTree *tree=(TTree*)f.Get("tree");
  TGraph *gr=(TGraph*)f.Get("ccqe_rate")->Clone();
  Int_t process=4;
  Double_t xs,dummy;
  gr->GetPoint(process,dummy,xs);
  delete gr;

  //NUANCE
  /*argh--chains are unusable in this ROOT version*/
  Int_t N_files=0;
  if (energy_point==300) N_files=40;
  if (energy_point==800) N_files=20;
  TFile *nu_files[N_files];
  TTree *nu_trees[N_files];
  for (Int_t iFile=0;iFile<N_files;iFile++) {
    TString file;
    if (energy_point==300) file=Form("../../nuance_from_fnal/oxnofp3/nuance%d.root",iFile+1);
    if (energy_point==800) file=Form("../../nuance_from_fnal/oxyg_no_fp2/nuance%d.root",iFile+1);
    nu_files[iFile]=new TFile(file);
    //nu_files[iFile]=new TFile(Form("$CONDOR_TMP/nuance/oxyg_no_fp2_300/nuance%d.root",iFile+1));
    //nu_files[iFile]=new TFile(Form("$CONDOR_TMP/nuance/oxyg_no_fp2/nuance%d.root",iFile+1));
    //nu_files[iFile]=new TFile(Form("$CONDOR_TMP/nuance/oxyg_no_fp_nor_eb/nuance%d.root",iFile+1));
    gROOT->cd();
    nu_trees[iFile]=(TTree*)(nu_files[iFile]->Get("h3"));
  }

  TString *tree_var=0;
  TString *nuance_var=0;
  TCut *nuance_cut=0;
  TString *label=0;

  Int_t N_bins=0;
  Double_t lower=0.0;
  Double_t upper=0.0;
  Double_t x_leg=0.0;

  switch (mode) {
  case 1:
    tree_var=new TString("kprime[0]");
    nuance_var=new TString("p_lepton[0][3]/1.0e3");
    nuance_cut=new TCut("cc && bound && channel==1");
    label=new TString(";E_{#mu} (GeV);d#sigma/dE_{#mu} (cm^{2}/GeV)");
    
    x_leg=0.15;
    N_bins=70;
    lower=0.12;
    upper=energy_point/1.0e3;
    break;
  case 2:
    tree_var=new TString("-q[0]**2+q[1]**2+q[2]**2+q[3]**2");
    nuance_var=new TString("-qsq/1.0e6");
    nuance_cut=new TCut("cc && bound && channel==1");
    label=new TString(";Q^{2} (GeV^{2});d#sigma/dQ^{2} (cm^{2}/GeV^{2})");

    x_leg=0.55;
    N_bins=70;
    lower=0.0;
    if (energy_point==300) upper=0.8;
    if (energy_point==800) upper=1.4;
    break;
  case 3:
    tree_var=new TString("mag_p");
    nuance_var=new TString("p_targ[4]/1.0e3");
    nuance_cut=new TCut("cc && bound && channel==1");
    label=new TString(";p_{i} (GeV/c);d#sigma/dp_{i} (cm^{2}/(GeV/c))");

    x_leg=0.15;
    N_bins=150;
    lower=0.0;
    upper=0.25;
    break;
  case 4:
    tree_var=new TString("sqrt(pprime[1]**2+pprime[2]**2+pprime[3]**2)");
    nuance_var=new TString("p_hadron[][4]/1.0e3");
    nuance_cut=new TCut("cc && bound && channel==1 && n_hadrons==1 && n_leptons==1");
    label=new TString(";p_{f} (GeV/c);d#sigma/dp_{f} (cm^{2}/(GeV/c))");

    x_leg=0.15;
    N_bins=150;
    lower=0.0;
    upper=0.8;
    break;
  case 5:
    tree_var=new TString("sqrt(pprime[1]**2+pprime[2]**2+pprime[3]**2)");
    nuance_var=new TString("sqrt((p_neutrino[0]-p_lepton[0][0]+p_targ[0])**2+(p_neutrino[1]-p_lepton[0][1]+p_targ[1])**2+(p_neutrino[2]-p_lepton[0][2]+p_targ[2])**2)/1.0e3");
    nuance_cut=new TCut("cc && bound && channel==1 && n_leptons==1");
    label=new TString(";p_{f} (GeV/c);d#sigma/dp_{f} (cm^{2}/(GeV/c))");

    x_leg=0.15;
    N_bins=150;
    lower=0.2;
    upper=0.6;
    //lower=0.2;
    //upper=1.3;
    break;
  }

  TCanvas canvas("canvas","",700,1000);
  canvas.Divide(1,2);
  canvas.cd(1);

  TH1D g0("g0",*label,N_bins,lower,upper);
  g0.Sumw2();
  gtree->Draw(*tree_var+">>g0","","goff");
  g0.Scale(gxs/g0.Integral("width"));
  g0.SetStats(kFALSE);
  g0.Draw("e");

  TH1D nh("nh",*label,N_bins,lower,upper);
  nh.Sumw2();
  nh.SetStats(kFALSE);
  for (Int_t iFile=0;iFile<N_files;iFile++) {
    nu_trees[iFile]->Draw(*nuance_var+">>+nh",*nuance_cut,"goff");
  }
  if (energy_point==300) nh.Scale(0.0033647*8*1.0e-36/nh.Integral("width"));
  if (energy_point==800) nh.Scale(0.0093045*8*1.0e-36/nh.Integral("width"));
  nh.SetLineColor(kRed);
  nh.SetMarkerColor(kRed);
  nh.Draw("esame");

  Int_t NHistos=5;
  TH1D *h[NHistos];
  TH1D *r[NHistos];
  for (Int_t iHisto=0;iHisto<NHistos;iHisto++) {
    h[iHisto]=new TH1D(Form("h%d",iHisto),*label,N_bins,lower,upper);
    h[iHisto]->Sumw2();
  }

  Int_t color=3;
  for (Int_t iHisto=0;iHisto<NHistos;iHisto++) {
    tree->Draw(*tree_var+Form(">>h%d",iHisto),Form("(MA_weights[%d]) * (process==%d)",iHisto,process),"goff");
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
  legend.SetHeader("S-M, dipole, FP=0");
  legend.AddEntry(&g0,"my events (A-S mode,MA=1.03)");
  legend.AddEntry(&nh,"NUANCE (MA=1.03)");

  Double_t MA_values[NHistos];
  for (Int_t iMA=0;iMA<NHistos;iMA++) {
    Double_t base=1.03;
    MA_values[iMA]=base+(iMA-2)*0.03;
  }

  Double_t factor=xs/h[1]->Integral("width");
  for (Int_t iHisto=0;iHisto<NHistos;iHisto++) {
    h[iHisto]->Scale(factor);
    h[iHisto]->Draw("esame"); 
    r[iHisto]=(TH1D*)h[iHisto]->Clone(Form("r%d",iHisto));
    legend.AddEntry(h[iHisto],Form("my events (NUANCE mode,MA=%#4.3g)",MA_values[iHisto]));
  }

  legend.Draw("same");
  
  canvas.cd(2);
  TH2D null("null",";;ratio to black",1,0.0,1.4,1,0.8,1.2);
  null.SetStats(kFALSE);
  null.Draw();

  for (Int_t iHisto=0;iHisto<NHistos;iHisto++) {
    //r[iHisto]->Scale(1.0/r[iHisto]->Integral());
    r[iHisto]->Divide(&g0);
    r[iHisto]->Draw("same");
  }
  canvas.cd(0);
  canvas.DrawClone();

  f.Close();

  delete tree_var;
  delete nuance_var;
  delete nuance_cut;
  delete label;
  
  for (Int_t iHisto=0;iHisto<NHistos;iHisto++) {
    delete h[iHisto];
    delete r[iHisto];
  }

  for (Int_t iFile=0;iFile<N_files;iFile++) {
    delete nu_files[iFile];
  }


}

////////////////////////////////////////////////////////////////////////
void MA_compare_800() {
  Init();

  //polish
  TFile f_polish("../ref/polish/O16/graph.root");
  gROOT->cd();
  
  TGraph *gr_pol_benhar=(TGraph*)f_polish.Get("benhar_q2")->Clone();
  TGraph *gr_pol_fg=(TGraph*)f_polish.Get("FG_q2")->Clone();
  f_polish.Close();


  //mine
  Int_t mode=2;

  //TFile f("../../events_5_nofp_800.root");
  TFile f("../../events.root");
  gROOT->cd();
  TTree *tree=(TTree*)f.Get("tree");
  TGraph *gr=(TGraph*)f.Get("ccqe_rate")->Clone();
  Int_t process=4;
  Double_t xs,dummy;
  gr->GetPoint(process,dummy,xs);
  delete gr;

  TString *tree_var=0;
  TString *label=0;

  Int_t N_bins=0;
  Double_t lower=0.0;
  Double_t upper=0.0;
  Double_t x_leg=0.0;

  switch (mode) {
  case 1:
    tree_var=new TString("kprime[0]");
    label=new TString(";E_{#mu} (GeV);d#sigma/dE_{#mu} (cm^{2}/GeV)");
    
    x_leg=0.15;
    N_bins=70;
    lower=0.12;
    upper=0.8;
    break;
  case 2:
    tree_var=new TString("-q[0]**2+q[1]**2+q[2]**2+q[3]**2");
    label=new TString(";Q^{2} (GeV^{2});d#sigma/dQ^{2} (cm^{2}/GeV^{2})");

    x_leg=0.55;
    N_bins=70;
    lower=0.0;
    upper=1.4;
    break;
  case 3:
    tree_var=new TString("mag_p");
    label=new TString(";p_{i} (GeV/c);d#sigma/dp_{i} (cm^{2}/(GeV/c))");

    x_leg=0.15;
    N_bins=150;
    lower=0.0;
    upper=0.25;
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
  }

  TCanvas canvas("canvas","",700,1000);
  canvas.Divide(1,2);
  canvas.cd(1);

  TH1D g_benhar("g_benhar",*label,N_bins,lower,upper);
  TH1D g_fg    ("g_fg",    *label,N_bins,lower,upper);
  g_benhar.Sumw2();
  g_fg.Sumw2();

  histo_graph(&g_benhar,gr_pol_benhar,4000000);
  histo_graph(&g_fg,gr_pol_fg,4000000);
  if (gr_pol_benhar) delete gr_pol_benhar;
  if (gr_pol_fg)     delete gr_pol_fg;

  g_fg.SetStats(kFALSE);
  g_benhar.SetStats(kFALSE);
  g_fg.SetLineColor(kRed);
  g_fg.SetMarkerColor(kRed);
  g_fg.Draw("e");
  g_benhar.Draw("esame");

  Int_t NHistos=5;
  TH1D *h[NHistos];
  TH1D *r[NHistos];
  for (Int_t iHisto=0;iHisto<NHistos;iHisto++) {
    h[iHisto]=new TH1D(Form("h%d",iHisto),*label,N_bins,lower,upper);
    h[iHisto]->Sumw2();
  }

  Int_t color=3;
  for (Int_t iHisto=0;iHisto<NHistos;iHisto++) {
    tree->Draw(*tree_var+Form(">>h%d",iHisto),Form("(MA_weights[%d]) * (process==%d)",iHisto,process),"goff");
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
  legend.SetHeader("dipole form factors, O^{16}, 800 MeV #nu_{#mu}");
  legend.AddEntry(&g_benhar,"A-S (Benhar, MA=1.03)");
  legend.AddEntry(&g_fg    ,"A-S (FG,     MA=1.03)");

  Double_t MA_values[NHistos];
  for (Int_t iMA=0;iMA<NHistos;iMA++) {
    Double_t base=1.03;
    MA_values[iMA]=base+(iMA-2)*0.03;
  }

  Double_t factor=xs/h[0]->Integral("width");
  for (Int_t iHisto=0;iHisto<NHistos;iHisto++) {
    h[iHisto]->Scale(factor);
    h[iHisto]->Draw("esame"); 
    r[iHisto]=(TH1D*)h[iHisto]->Clone(Form("r%d",iHisto));
    legend.AddEntry(h[iHisto],Form("my events (NUANCE mode,MA=%#4.3g)",MA_values[iHisto]));
  }
  TH1D *r_fg=(TH1D*)g_fg.Clone("r_fg");
  legend.Draw("same");
  
  canvas.cd(2);
  TH2D null("null",";;ratio to A-S Benhar",1,lower,upper,1,0.8,2.0);
  null.SetStats(kFALSE);
  null.Draw();

  r_fg->Divide(&g_benhar);
  r_fg->Draw("same");
  for (Int_t iHisto=0;iHisto<NHistos;iHisto++) {
    //r[iHisto]->Scale(1.0/r[iHisto]->Integral());
    r[iHisto]->Divide(&g_benhar);
    r[iHisto]->Draw("same");
  }
  canvas.cd(0);
  canvas.DrawClone();

  f.Close();

  delete tree_var;
  delete label;
  
  for (Int_t iHisto=0;iHisto<NHistos;iHisto++) {
    delete h[iHisto];
    delete r[iHisto];
  }

}

////////////////////////////////////////////////////////////////////////
void MA_compare(Int_t mode,Int_t energy_point) {
  Init();

  //NUANCE mode
  TFile f4(Form("../../events_4_%d.root",energy_point));
  gROOT->cd();
  TTree *tree4=(TTree*)f4.Get("tree");
  tree4->SetName("tree4");
  TGraph *gr4=(TGraph*)f4.Get("ccqe_rate")->Clone();
  Double_t xs4,dummy;
  gr4->GetPoint(4,dummy,xs4);
  delete gr4;

  //A-S SF mode
  TFile f23(Form("../../events_23_%d.root",energy_point));
  gROOT->cd();
  TTree *tree23=(TTree*)f23.Get("tree");
  tree23->SetName("tree23");
  TGraph *gr23=(TGraph*)f23.Get("ccqe_rate")->Clone();
  Double_t xs23=0.0;
  {
    Double_t dummy_x,dummy_y;
    gr23->GetPoint(2,dummy_x,dummy_y);
    xs23+=dummy_y;
    gr23->GetPoint(3,dummy_x,dummy_y);
    xs23+=dummy_y;
  }
  delete gr23;
  
  TString *tree_var=0;
  TString *label=0;

  Int_t N_bins=0;
  Double_t lower=0.0;
  Double_t upper=0.0;
  Double_t x_leg=0.0;

  switch (mode) {
  case 1:
    tree_var=new TString("kprime[0]");
    label=new TString(";E_{#mu} (GeV);d#sigma/dE_{#mu} (cm^{2}/GeV)");
    
    x_leg=0.15;
    N_bins=70;
    lower=0.12;
    upper=0.8;
    break;
  case 2:
    tree_var=new TString("-q[0]**2+q[1]**2+q[2]**2+q[3]**2");
    label=new TString(";Q^{2} (GeV^{2});d#sigma/dQ^{2} (cm^{2}/GeV^{2})");

    x_leg=0.55;
    N_bins=70;
    lower=0.0;
    upper=1.4;
    break;
  case 3:
    tree_var=new TString("Q2qe");
    label=new TString(";Q^{2}_{QE} (GeV^{2});d#sigma/dQ^{2} (cm^{2}/GeV^{2})");

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
  }

  Int_t NHistos=4;
  TH1D *h[NHistos];
  TH1D *r[NHistos];
  for (Int_t iHisto=0;iHisto<NHistos;iHisto++) {
    h[iHisto]=new TH1D(Form("h%d",iHisto),*label,N_bins,lower,upper);
    h[iHisto]->Sumw2();
  }
  TH1D as23("as23",*label,N_bins,lower,upper);
  as23.Sumw2();
  as23.SetStats(kFALSE);
  tree23->Draw(*tree_var+">>as23","process==2 || process==3","goff");
  as23.Scale(xs23/as23.Integral("width"));

  Int_t color=3;
  for (Int_t iHisto=0;iHisto<NHistos;iHisto++) {
    tree4->Draw(*tree_var+Form(">>h%d",iHisto),Form("(MA_weights[%d]) * (process==4)",iHisto),"goff");
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
  legend.AddEntry(&as23,"Ankowski-Sobczyk S.F., MA=1.03");

  Double_t MA_values[NHistos];
  for (Int_t iMA=0;iMA<NHistos;iMA++) {
    Double_t base=1.03;
    MA_values[iMA]=base+(iMA-1)*0.1;
  }

  TCanvas canvas("canvas","",700,1000);
  canvas.Divide(1,2);
  canvas.cd(1);
  
  as23.Draw("e");

  Double_t factor=xs4/h[0]->Integral("width");
  for (Int_t iHisto=0;iHisto<NHistos;iHisto++) {
    h[iHisto]->Scale(factor);
    h[iHisto]->Draw("esame"); 
    r[iHisto]=(TH1D*)h[iHisto]->Clone(Form("r%d",iHisto));
    legend.AddEntry(h[iHisto],Form("NUANCE mode, MA=%#4.3g",MA_values[iHisto]));
  }
  legend.Draw("same");
  
  canvas.cd(2);
  TH2D null("null",";;ratio to A-S SF",1,lower,upper,1,0.8,2.0);
  null.SetStats(kFALSE);
  null.Draw();

  for (Int_t iHisto=0;iHisto<NHistos;iHisto++) {
    //r[iHisto]->Scale(1.0/r[iHisto]->Integral());
    r[iHisto]->Divide(&as23);
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

////////////////////////////////////////////////////////////////////////
void q2_look() {
  TFile f("../../events.root");
  gROOT->cd();
  TTree *tree=(TTree*)f.Get("tree");

  tree->Draw("Q2qe:-q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3]","Q2qe<2.0");

  TF1 g("g","x",0.0,1.1);
  g.SetLineWidth(1);
  g.SetLineColor(kRed);
  g.Draw("same");
}

////////////////////////////////////////////////////////////////////////
void compare1() {
  Init();

  //mine
  TFile f1("../../events_1.root");
  TFile f4("../../events_4.root");
  gROOT->cd();
  TTree *tree1=(TTree*)f1.Get("tree");
  TGraph *gr1=(TGraph*)f1.Get("ccqe_rate");

  TTree *tree4=(TTree*)f4.Get("tree");
  TGraph *gr4=(TGraph*)f4.Get("ccqe_rate");

  
  Int_t N_bins=10;
  Double_t lower=0.0;
  Double_t upper=1.2;

  TString var="q[1]**2+q[2]**2+q[3]**2-q[0]**2";
  TString title=";Q^{2} (GeV^{2});d#sigma/dQ^{2} (cm^{2}/GeV^{2})";

  TH1D h1("h1",title,N_bins,lower,upper);
  h1.Sumw2();
  TH1D h4("h4",title,N_bins,lower,upper);
  h4.Sumw2();
  TH1D hr("hr",title,N_bins,lower,upper);
  hr.Sumw2();
  hr.GetYaxis()->SetTitle("ratio");
  
  Double_t dummy,xs1,xs4;

  gr1->GetPoint(1,dummy,xs1);
  gr4->GetPoint(4,dummy,xs4);

  tree1->Draw(var+">>h1","","goff");
  h1.Scale(xs1/h1.Integral("width"));
  tree4->Draw(var+">>h4","","goff");
  h4.Scale(xs4/h4.Integral("width"));

  hr.Divide(&h1,&h4);
  hr.DrawClone("");
}

////////////////////////////////////////////////////////////////////////
void compare6() {
  Init();

  TFile f1("../../events_3.root");
  TTree *tree1=(TTree*)f1.Get("tree");
  TGraph *gr1=(TGraph*)f1.Get("ccqe_rate")->Clone();
  Int_t process1=4;
  Double_t xs1,dummy;
  gr1->GetPoint(process1,dummy,xs1);
  delete gr1;

  TFile f2("../../events_4.root");
  TTree *tree2=(TTree*)f2.Get("tree");
  TGraph *gr2=(TGraph*)f2.Get("ccqe_rate")->Clone();
  Int_t process2=1;
  Double_t xs2;
  gr2->GetPoint(process2,dummy,xs2);
  delete gr2;

  Double_t x_leg=0.15;
  Int_t N_bins=100;
  Double_t lower=0.0;
  Double_t upper=1.2;

  TH1D h1("h1","",N_bins,lower,upper);
  TH1D h2("h2","",N_bins,lower,upper);
  h1.Sumw2();
  h2.Sumw2();
  TH1D *contribs=(TH1D*)h1.Clone("contribs");

  tree1->Draw("-q[0]**2+q[1]**2+q[2]**2+q[3]**2>>h1",Form("process==%d",process1),"goff");
  tree2->Draw("-q[0]**2+q[1]**2+q[2]**2+q[3]**2>>h2",Form("process==%d",process2),"goff");
  h1.Scale(xs1/h1.Integral("width"));
  h2.Scale(xs2/h2.Integral("width"));

  TCanvas canvas("canvas","",700,1000);

  canvas.Divide(1,2);
  canvas.cd(1);

  h1.SetStats(kFALSE);
  h1.GetXaxis()->CenterTitle();
  h1.GetYaxis()->CenterTitle();
  h1.GetYaxis()->SetTitleOffset(1.2);
  h1.SetLineColor(kRed);
  h1.SetMarkerColor(kRed);
  h1.Draw("e");

  h2.SetLineColor(kBlue);
  h2.SetMarkerColor(kBlue);
  h2.Draw("esame");

  TLegend legend(x_leg,0.8,x_leg+0.35,0.9);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetHeader("S-M, dipole, FP=0");
  legend.AddEntry(&h1,"my events (no OS,Ep-Eb)");
  legend.AddEntry(&h2,"NUANCE");
  legend.Draw("same");

  chi2_compare(&h1,&h2,1,N_bins,contribs);

  //TLatex text(0.20,0.92,comp_str);
  //text.SetNDC();
  //text.DrawClone();

  canvas.cd(2);
  contribs->Draw();

  canvas.cd(0);
  canvas.DrawClone();

  delete contribs;

}




  //if (mode==200) {
  //  Int_t N_bins=10;
  //  Double_t lower=0.0;
  //  Double_t upper=1.4;
  //  TH1D h_sm("h_sm","#nu_{#mu} + ^{16}O #rightarrow #mu + p + X (E_{#nu}=0.8 GeV);Q^{2} (GeV^{2});events/bin/arb. POT",N_bins,lower,upper);
  //  h_sm.Sumw2();
  //  TH1D h_as("h_as","#nu_{#mu} + ^{16}O #rightarrow #mu + p + X (E_{#nu}=0.8 GeV);Q^{2} (GeV^{2});events/bin/arb .POT",N_bins,lower,upper);
  //  h_as.Sumw2();
  //  TH1D h_ratio("h_ratio","#nu_{#mu} + ^{16}O #rightarrow #mu + p + X (E_{#nu}=0.8 GeV);Q^{2} (GeV^{2});ratio",N_bins,lower,upper);
  //  h_ratio.Sumw2();
  //
  //  tree->Draw("-q[0]**2+q[1]**2+q[3]**2>>h_sm","process==1","goff");
  //  tree->Draw("-q[0]**2+q[1]**2+q[3]**2>>h_as","process>1","goff");
  //
  //  h_ratio.Divide(&h_sm,&h_as);
  //
  //  h_ratio.SetStats(kFALSE);
  //  h_ratio.GetXaxis()->CenterTitle();
  //  h_ratio.GetYaxis()->CenterTitle();
  //  h_ratio.GetYaxis()->SetTitleOffset(1.2);
  //  h_ratio.SetLineColor(kRed);
  //  h_ratio.SetMarkerColor(kRed);
  //  h_ratio.DrawCopy("e");
  //
  //}
