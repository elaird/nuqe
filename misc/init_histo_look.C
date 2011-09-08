#include <iostream>
#include "TROOT.h"
#include "TFile.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TCanvas.h"

void init_histo_look() {
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(111110);
  TFile f("../events/events_123_5000.root");
  gROOT->cd();
  TH1D *h=(TH1D*)f.Get("init_histo_proc3_bin0");

  TCanvas can;
  can.SetLogy();
  h->Draw();
  can.DrawClone();
}
