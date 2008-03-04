void gen_look(Int_t mode) {
  gROOT->SetStyle("Plain");

  TFile f("xs.root");
  gROOT->cd();

  if (mode==-10) {
    TTree *tree=(TTree*)f.Get("tree")->Clone();
    f.Close();
    
    TH1D h("h","",8,-TMath::Pi(),TMath::Pi());
    h.Sumw2();
    h.SetMinimum(0.0);
    h.SetStats(kFALSE);
    
    tree->Draw("phi>>h","","goff");
    h.Fit("pol1","f");
    h.DrawCopy("e");

    delete tree;
  }

  if (mode==10) {
    TTree *tree=(TTree*)f.Get("tree")->Clone();
    f.Close();
    
    tree->Draw("inv_weight","inv_weight<1.0");

    delete tree;
  }

  if (mode==20) {
    gROOT->ProcessLine("TFile f(\"xs.root\");");
    gROOT->ProcessLine("TTree *tree=(TTree*)f.Get(\"tree\")->Clone();");
  }

  if (mode==30) {
    TH1D *hv =(TH1D*)f.Get("integral_histo")->Clone();
    TH1D *hvc=(TH1D*)f.Get("integral_check_histo")->Clone();
    f.Close();

    hvc->Sumw2();
    hv->Scale(1.0/hv->Integral());
    //hvc->Scale(1.0/hvc->Integral());

    hv->DrawClone("e");

    delete hv;
    delete hvc;
  }

  if (mode==40) {
    TH1D *hv =(TH1D*)f.Get("integral_histo")->Clone();
    TH1D *hvc=(TH1D*)f.Get("integral_check_histo")->Clone();
    f.Close();

    hvc->Sumw2();
    hv->Scale(1.0/hv->Integral());
    hvc->Scale(1.0/hvc->Integral());

    hvc->Divide(hv);
    Int_t i_max=hvc->GetNbinsX();
    Double_t chi2=0.0;
    Int_t bins_counted=0;
    for (Int_t iBin=1;iBin<=i_max;iBin++) {
      Double_t content=hvc->GetBinContent(iBin);
      Double_t error=hvc->GetBinError(iBin);
      if (content>0.0) {
	bins_counted++;
	Double_t value=(1.0-content)/error;
	chi2+=pow(value,2);
      }
    }
    hvc->DrawClone("e");
    Double_t prob=TMath::Prob(chi2,bins_counted);
    printf("counted: %d,chi2: %g,prob: %g\n",bins_counted,chi2,prob);
    delete hv;
    delete hvc;
  }

  if (mode==50) {
    TFile f_polish("polish_plots/O16_E_SM/graph.root");
    gROOT->cd();
    TTree *tree=(TTree*)f.Get("tree")->Clone();
    TGraph *gr_pol=(TGraph*)f_polish.Get("polish_Q2_sm")->Clone();

    f.Close();
    f_polish.Close();

    Int_t N_bins=50;
    Double_t lower=0.0;
    Double_t upper=1.4;
    TH1D h("h","#nu_{#mu} + O^{16} #rightarrow #mu + p + O^{15} (E_{#nu}=0.8 GeV);Q^{2} (GeV^{2});d#sigma/dQ^{2} (cm^{2}/GeV^{2})",N_bins,lower,upper);h.Sumw2();

    tree->Draw("q[1]**2+q[2]**2+q[3]**2-q[0]**2>>h","(1.0/inv_weight)","goff");
    h.Scale(7.287e-38/h.Integral("width"));

    h.SetStats(kFALSE);
    h.GetXaxis()->CenterTitle();
    h.GetYaxis()->CenterTitle();
    h.GetYaxis()->SetTitleOffset(1.2);
    h.SetLineColor(kRed);
    h.SetMarkerColor(kRed);
    h.DrawCopy("e");

    gr_pol->SetTitle("#nu_{#mu} + O^{16} #rightarrow #mu + p + O^{15} (E_{#nu}=0.8 GeV);Q^{2} (GeV^{2});d#sigma/dQ^{2} (cm^{2}/GeV^{2})");
    gr_pol->GetXaxis()->CenterTitle();
    gr_pol->GetYaxis()->CenterTitle();
    gr_pol->GetYaxis()->SetTitleOffset(1.2);

    gr_pol->SetMarkerStyle(20);
    gr_pol->SetMarkerSize(0.3);
    gr_pol->SetLineColor(kBlack);
    gr_pol->SetMarkerColor(kBlack);
    gr_pol->DrawClone("lpsame");

    TLegend legend(0.55,0.8,0.9,0.9);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.AddEntry(gr_pol,"Ankowski-Sobczyk paper (S-M)","p");
    legend.AddEntry(&h,"my generated events (S-M)","p");
    legend.DrawClone("same");

    delete gr_pol;
    delete tree;
  }

  if ((mode/10)==6) {
    TFile f_polish("polish_plots/O16_E_SM/graph.root");
    gROOT->cd();
    TTree *tree=(TTree*)f.Get("tree")->Clone();
    TGraph *gr_pol=0;
    
    TString pb_string;
    if (mode==60) {
      pb_string=" (S-M)";
      gr_pol=(TGraph*)f_polish.Get("polish_e_sm_pb")->Clone();
    }
    else {
      pb_string=" (S-M), no PB";
      gr_pol=(TGraph*)f_polish.Get("polish_e_sm_no_pb")->Clone();
    }

    f.Close();
    f_polish.Close();

    Int_t N_bins=100;
    Double_t lower=0.0;
    Double_t upper=0.8;
    TH1D h("h","#nu_{#mu} + O^{16} #rightarrow #mu + p + O^{15} (E_{#nu}=0.8 GeV);E_{#mu} (GeV);d#sigma/dE_{#mu} (cm^{2}/GeV)",N_bins,lower,upper);
    h.Sumw2();

    tree->Draw("0.8-q[0]>>h","(1.0/inv_weight)","goff");
    h.Scale(7.287e-38/h.Integral("width"));

    h.SetStats(kFALSE);
    h.GetXaxis()->CenterTitle();
    h.GetYaxis()->CenterTitle();
    h.GetYaxis()->SetTitleOffset(1.2);
    h.SetLineColor(kRed);
    h.SetMarkerColor(kRed);
    h.DrawCopy("e");

    //gr_pol->SetTitle("#nu_{#mu} + O^{16} #rightarrow #mu + p + O^{15} (E_{#nu}=0.8 GeV);Q^{2} (GeV^{2});d#sigma/dQ^{2} (cm^{2}/GeV^{2})");
    //gr_pol->GetXaxis()->CenterTitle();
    //gr_pol->GetYaxis()->CenterTitle();
    //gr_pol->GetYaxis()->SetTitleOffset(1.2);

    gr_pol->SetMarkerStyle(20);
    gr_pol->SetMarkerSize(0.3);
    gr_pol->SetLineColor(kBlack);
    gr_pol->SetMarkerColor(kBlack);
    gr_pol->DrawClone("lpsame");

    TLegend legend(0.1,0.8,0.45,0.9);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.AddEntry(gr_pol,"Ankowski-Sobczyk paper"+pb_string,"p");
    legend.AddEntry(&h,"my generated events"+pb_string,"p");
    legend.DrawClone("same");

    delete gr_pol;
    delete tree;
  }


}
