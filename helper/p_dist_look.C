void p_dist_look(Int_t mode=0) {
  gROOT->SetStyle("Plain");

  TFile f("p_dists.root");
  gROOT->cd();
  TGraph *sm=(TGraph*)f.Get("sm")->Clone();
  TGraph *mf=(TGraph*)f.Get("mf")->Clone();
  TGraph *corr=(TGraph*)f.Get("corr")->Clone();
  TGraph *tot=(TGraph*)f.Get("tot")->Clone();
  f.Close();

  TFile g("../ref/polish/O16/graph.root");
  gROOT->cd();
  TGraph *as_mf=(TGraph*)g.Get("polish_n_mf")->Clone();
  TGraph *as_corr=(TGraph*)g.Get("polish_n_corr")->Clone();
  f.Close();

  if (mode==1) {
    tot->DrawClone("ap");
    sm->DrawClone("psame");
    mf->DrawClone("psame");
    corr->DrawClone("psame");
    gPad->SetLogy();
  }

  if (mode==2) {
    mf->DrawClone("ap");
    corr->DrawClone("psame");
    as_mf->SetMarkerColor(kRed);
    as_mf->DrawClone("psame");
    as_corr->SetMarkerColor(kBlue);
    as_corr->DrawClone("psame");
    gPad->SetLogy();
    
  }

  if (mode==3) {
    Int_t n_points=tot->GetN();
    Int_t n_bins=n_points;

    Double_t p_max,y;
    tot->GetPoint(n_points-1,p_max,y);

    TH1D h1("h1","",n_bins,0.0,p_max);
    TH1D h2("h2","",n_bins,0.0,p_max);
    TH1D h3("h3","",n_bins,0.0,p_max);
    TH1D h4("h4","",n_bins,0.0,p_max);

    for (Int_t iBin=1;iBin<=n_bins;iBin++) {
      Double_t value1=sm->Eval(h1.GetBinCenter(iBin));
      value1*=4.0*TMath::Pi()*pow(h1.GetBinCenter(iBin),2);
      h1.SetBinContent(iBin,value1);

      Double_t value2=mf->Eval(h2.GetBinCenter(iBin));
      value2*=4.0*TMath::Pi()*pow(h2.GetBinCenter(iBin),2);
      h2.SetBinContent(iBin,value2);

      Double_t value3=tot->Eval(h3.GetBinCenter(iBin));
      value3*=4.0*TMath::Pi()*pow(h3.GetBinCenter(iBin),2);
      h3.SetBinContent(iBin,value3);

      Double_t value4=corr->Eval(h4.GetBinCenter(iBin));
      value4*=4.0*TMath::Pi()*pow(h4.GetBinCenter(iBin),2);
      h4.SetBinContent(iBin,value4);
    }

    h1.SetStats(kFALSE);
    h1.SetTitle(";p (GeV);N_{target nucleons}/GeV");
    h2.SetLineColor(kRed);
    h3.SetLineColor(kBlue);
    h4.SetLineColor(kGreen);

    h1.DrawCopy();
    h2.DrawCopy("same");
    h4.DrawCopy("same");
    h3.DrawCopy("same");

    TLegend legend(0.45,0.6,0.9,0.9,"^{16}O");
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.AddEntry(&h1,"Smith-Moniz","l");
    legend.AddEntry(&h2,"Ciofi degli Atti-Simula mean field (MF)","l");
    legend.AddEntry(&h4,"Ciofi degli Atti-Simula correlated (corr)","l");
    legend.AddEntry(&h3,"Ciofi degli Atti-Simula total","l");
    legend.DrawClone("same");

    gPad->SetLogy(0);
    cout << "h1.int: " << h1.Integral("width") << endl;
    cout << "h2.int: " << h2.Integral("width") << endl;
    cout << "h3.int: " << h3.Integral("width") << endl;
    cout << "h4.int: " << h4.Integral("width") << endl;
  }
  
  delete sm;
  delete mf;
  delete corr;
  delete tot;
}
