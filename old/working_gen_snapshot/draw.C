Int_t draw(Int_t mode=0) {
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gEnv->SetValue("Root.Stacktrace","no");

  TFile f("~elaird/ccqe/xs_ted/xs.root");
  gROOT->cd();

  if (mode==1)  {
    TGraph2D *gr=(TGraph2D*)f.Get("xs2d")->Clone();
    gr->DrawClone("p");
  }

  if (mode==2)  {
    TGraph *gr=(TGraph*)f.Get("xs1d")->Clone();
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.3);
    gr->DrawClone("ap");
  }

  if (mode==3) {
    TF1 g("g","[0]*sqrt(1.0-x*x)+[1]*x-[2]",-1.0,1.0);
    g.SetNpx(1000);
    //example of low problem
    g.SetParameters(-0.0982392,0.437679,-0.442767);
    //example of high problem
    //g.SetParameters(0.184844,1.49364,1.49455);
    g.SetLineWidth(1);
    g.DrawClone();
    
  }

  f.Close();

  return 0;
}
