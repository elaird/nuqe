Int_t compare_plots(Int_t mode=0) {
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gEnv->SetValue("Root.Stacktrace","no");


  if (mode==1) {
    TFile f("xs.root");
    TFile f_TK("~elaird/xs_ted/single_n/katori/single.root");
    TFile f_single_mod("~elaird/xs_ted/single_mod.root");
    TFile f_single("~elaird/xs_ted/single.root");
    gROOT->cd();

    TGraph *gr           =(TGraph*)f.Get("xs1d")->Clone();
    TGraph *gr_single_mod=(TGraph*)f_single_mod.Get("xs1d")->Clone();
    TGraph *gr_single    =(TGraph*)f_single.    Get("xs1d")->Clone();
    TGraph *gr_TK        =(TGraph*)f_TK.Get("xs_Q2")->Clone();

    f.Close();
    f_TK.Close();
    f_single.Close();
    f_single_mod.Close();

    gr_TK->SetTitle("#nu_{#mu} + n #rightarrow #mu + p (E_{#nu}=0.8 GeV);Q^{2} (GeV^{2});d#sigma/dQ^{2} (cm^{2}/GeV^{2})");
    gr_TK->GetXaxis()->CenterTitle();
    gr_TK->GetYaxis()->CenterTitle();
    //gr_TK->GetYaxis()->SetTitleOffset(1.2);

    gr_TK->SetMarkerStyle(20);
    gr_TK->SetMarkerSize(0.3);
    gr_TK->DrawClone("ap");

    //gr_single_mod->SetMarkerStyle(20);
    //gr_single_mod->SetMarkerSize(0.3);
    //gr_single_mod->SetMarkerColor(kRed);
    //gr_single_mod->DrawClone("psame");
    //
    //gr_single->SetMarkerStyle(20);
    //gr_single->SetMarkerSize(0.3);
    //gr_single->SetMarkerColor(kBlue);
    //gr_single->DrawClone("psame");

    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.3);
    gr->SetMarkerColor(kGreen);
    gr->DrawClone("psame");

    TLegend legend(0.4,0.7,0.9,0.9);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.AddEntry(gr_TK,"LLewellyn-Smith (BooNE memo)","p");
    legend.AddEntry(gr,"my calculation (S-M, p_{F}=3 MeV, E_{b}=0 MeV)","p");
    legend.DrawClone("same");

    delete gr;
    delete gr_single_mod;
    delete gr_single;
    delete gr_TK;
  }

  if (mode==2) {

    TFile f("~elaird/ccqe/xs_ted/xs.root");
    TFile f_polish("~elaird/ccqe/xs_ted/polish_plots/O16_E_SM/graph.root");
    gROOT->cd();

    TGraph *gr           =(TGraph*)f.Get("xs1d")->Clone();
    TGraph *gr_pol1=(TGraph*)f_polish.Get("polish_e_sm_pb")->Clone();
    TGraph *gr_pol2=(TGraph*)f_polish.Get("polish_e_sm_no_pb")->Clone();

    f.Close();
    f_polish.Close();

    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.3);
    gr->SetMarkerColor(kMagenta);
    gr->DrawClone("ap");

    gr_pol1->SetMarkerStyle(20);
    gr_pol1->SetMarkerSize(0.3);
    gr_pol1->SetLineColor(kBlack);
    gr_pol1->SetMarkerColor(kBlack);
    gr_pol1->DrawClone("lpsame");

    gr_pol2->SetMarkerStyle(20);
    gr_pol2->SetMarkerSize(0.3);
    gr_pol2->SetLineColor(kBlue);
    gr_pol2->SetMarkerColor(kBlue);
    gr_pol2->DrawClone("lpsame");

    delete gr;
    delete gr_pol1;
    delete gr_pol2;
  }

  if (mode==3) {
    TFile f("~elaird/ccqe/xs_ted/xs.root");//ref_plots/Q2_SM.root");
    TFile f_polish("~elaird/ccqe/xs_ted/polish_plots/O16_E_SM/graph.root");
    gROOT->cd();

    TGraph *gr    =(TGraph*)f.Get("xs1d")->Clone();
    TGraph *gr_pol=(TGraph*)f_polish.Get("polish_Q2_sm")->Clone();

    f.Close();
    f_polish.Close();

    gr_pol->SetTitle("#nu_{#mu} + O^{16} #rightarrow #mu + p + O^{15} (E_{#nu}=0.8 GeV);Q^{2} (GeV^{2});d#sigma/dQ^{2} (cm^{2}/GeV^{2})");
    gr_pol->GetXaxis()->CenterTitle();
    gr_pol->GetYaxis()->CenterTitle();
    gr_pol->GetYaxis()->SetTitleOffset(1.2);

    gr_pol->SetMarkerStyle(20);
    gr_pol->SetMarkerSize(0.3);
    gr_pol->SetLineColor(kBlack);
    gr_pol->SetMarkerColor(kBlack);
    gr_pol->DrawClone("alp");

    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.3);
    gr->SetMarkerColor(kMagenta);
    gr->DrawClone("psame");

    TLegend legend(0.55,0.8,0.9,0.9);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.AddEntry(gr_pol,"Ankowski-Sobczyk paper (S-M)","p");
    legend.AddEntry(gr,"my calculation (S-M)","p");
    legend.DrawClone("same");

    delete gr_pol;
    delete gr;
  }

  if (mode==4) {
    TFile f_polish("~elaird/ccqe/xs_ted/polish_plots/O16_E_SM/graph.root");
    gROOT->cd();
    TGraph *gr_pol=(TGraph*)f_polish.Get("polish_Q2_sm")->Clone();
    f_polish.Close();

    Int_t N_bins=150;
    Double_t lower=0.0;
    Double_t upper=1.4;
    TH1D h1("h1","#nu_{#mu} + O^{16} #rightarrow #mu + p + O^{15} (E_{#nu}=0.8 GeV);Q^{2} (GeV^{2});d#sigma/dQ^{2} (cm^{2}/GeV^{2})",N_bins,lower,upper);h1.Sumw2();
    TH1D h2("h2","#nu_{#mu} + O^{16} #rightarrow #mu + p + O^{15} (E_{#nu}=0.8 GeV);Q^{2} (GeV^{2});d#sigma/dQ^{2} (cm^{2}/GeV^{2})",N_bins,lower,upper);h2.Sumw2();

    TChain chain("chain");
    Int_t N_files=5;
    for (Int_t iFile=0;iFile<N_files;iFile++) {
      chain.Add(Form("$CONDOR_TMP/nuance/oxyg/nuance%d.root/h3",iFile+1));
    }
    chain.Draw("-qsq/1.0e6>>h2","cc && bound && channel==1","goff");
    h2.Scale(0.0093004*8*1.0e-36/h2.Integral("width"));

    //for (Int_t iBin=1;iBin<=N_bins;iBin++) {
    //  Double_t content=gr_pol->Eval(h1.GetBinCenter(iBin));
    //  h1.SetBinContent(iBin,content);
    //}

    h2.SetStats(kFALSE);
    h2.GetXaxis()->CenterTitle();
    h2.GetYaxis()->CenterTitle();
    h2.GetYaxis()->SetTitleOffset(1.2);
    h2.SetLineColor(kRed);
    h2.SetMarkerColor(kRed);
    h2.DrawCopy("e");

    gr_pol->SetMarkerStyle(20);
    gr_pol->SetMarkerSize(0.3);
    gr_pol->SetLineColor(kBlack);
    gr_pol->SetMarkerColor(kBlack);
    gr_pol->DrawClone("lpsame");


    TLegend legend(0.55,0.8,0.9,0.9);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.AddEntry(&h1,"Ankowski-Sobczyk paper (S-M, dipole)");
    legend.AddEntry(&h2,"NUANCE (S-M, dipole)");
    legend.DrawClone("same");
    
    delete gr_pol;
  }

  return 0;
}
