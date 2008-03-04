Int_t draw(Int_t mode=0) {
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gEnv->SetValue("Root.Stacktrace","no");

  TFile f("single.root");
  gROOT->cd();

  if (mode==1)  {
    TGraph *gr=(TGraph*)f.Get("xs_Q2_0")->Clone();
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.3);
    gr->DrawClone("ap");
    Double_t x,y;
    gr->GetPoint(0,x,y);
    cout << x << "," << y << endl;
  }

  if (mode==2)  {
    TGraph *gr=(TGraph*)f.Get("xs_E_mu_0")->Clone();
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.3);
    gr->DrawClone("ap");
  }

  if (mode==3)  {
    TGraph *gr=(TGraph*)f.Get("XS")->Clone();
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.3);
    gr->DrawClone("ap");
  }

  f.Close();

  return 0;
}
