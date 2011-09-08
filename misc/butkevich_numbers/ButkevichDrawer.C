#include "TGraph.h"
#include "TROOT.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "ButkevichDrawer.h"
///////////////////////////////////////////////////////////////////
ButkevichDrawer::~ButkevichDrawer () {
}
///////////////////////////////////////////////////////////////////
ButkevichDrawer::ButkevichDrawer () {
  x0_ph_=0.0;
  x1_ph_=2.0;

  y0_ph_=0.0;
  y1_ph_=14.0e-38;

  int i=0;
  graph_fg_ma10_.SetPoint(i++,0.0199999996,7.27387412 );
  graph_fg_ma10_.SetPoint(i++,0.0399999991,9.55174007 );
  graph_fg_ma10_.SetPoint(i++,0.0599999987,10.9884588 );
  graph_fg_ma10_.SetPoint(i++,0.0799999982,11.7327709 );
  graph_fg_ma10_.SetPoint(i++,0.100000001 ,12.2761743 );
  graph_fg_ma10_.SetPoint(i++,0.119999997 ,12.3691472 );
  graph_fg_ma10_.SetPoint(i++,0.140000001 ,12.3633703 );
  graph_fg_ma10_.SetPoint(i++,0.159999996 ,12.1793533 );
  graph_fg_ma10_.SetPoint(i++,0.180000007 ,11.8525565 );
  graph_fg_ma10_.SetPoint(i++,0.200000003 ,11.4479998 );
  graph_fg_ma10_.SetPoint(i++,0.219999999 ,10.9937701 );
  graph_fg_ma10_.SetPoint(i++,0.25        ,10.3420188 );
  graph_fg_ma10_.SetPoint(i++,0.300000012 ,9.33688934 );
  graph_fg_ma10_.SetPoint(i++,0.400000006 ,7.61881705 );
  graph_fg_ma10_.SetPoint(i++,0.5         ,6.24334608 );
  graph_fg_ma10_.SetPoint(i++,0.600000024 ,5.14622892 );
  graph_fg_ma10_.SetPoint(i++,0.699999988 ,4.26946668 );
  graph_fg_ma10_.SetPoint(i++,0.800000012 ,3.56562251 );
  graph_fg_ma10_.SetPoint(i++,0.899999976 ,2.99736448 );
  graph_fg_ma10_.SetPoint(i++,1.          ,2.53577786 );
  graph_fg_ma10_.SetPoint(i++,1.10000002  ,2.15864362 );
  graph_fg_ma10_.SetPoint(i++,1.20000005  ,1.79060548 );
  graph_fg_ma10_.SetPoint(i++,1.29999995  ,1.40638075 );
  graph_fg_ma10_.SetPoint(i++,1.39999998  ,1.0499232  );
  graph_fg_ma10_.SetPoint(i++,1.5         ,0.742385396);
  graph_fg_ma10_.SetPoint(i++,1.60000002  ,0.491628114);
  graph_fg_ma10_.SetPoint(i++,1.70000005  ,0.298013103);
  graph_fg_ma10_.SetPoint(i++,1.79999995  ,0.157894259);

  i=0;
  graph_fg_ma12_.SetPoint(i++,0.0199999996,7.38005844 );
  graph_fg_ma12_.SetPoint(i++,0.0399999991,9.82302812 );
  graph_fg_ma12_.SetPoint(i++,0.0599999987,11.4447471 );
  graph_fg_ma12_.SetPoint(i++,0.0799999982,12.3669387 );
  graph_fg_ma12_.SetPoint(i++,0.100000001 ,13.0870021 );
  graph_fg_ma12_.SetPoint(i++,0.119999997 ,13.3286576 );
  graph_fg_ma12_.SetPoint(i++,0.140000001 ,13.4595955 );
  graph_fg_ma12_.SetPoint(i++,0.159999996 ,13.389613  );
  graph_fg_ma12_.SetPoint(i++,0.180000007 ,13.1529349 );
  graph_fg_ma12_.SetPoint(i++,0.200000003 ,12.818601  );
  graph_fg_ma12_.SetPoint(i++,0.219999999 ,12.4166647 );
  graph_fg_ma12_.SetPoint(i++,0.25        ,11.82557   );
  graph_fg_ma12_.SetPoint(i++,0.300000012 ,10.8825037 );
  graph_fg_ma12_.SetPoint(i++,0.400000006 ,9.18569804 );
  graph_fg_ma12_.SetPoint(i++,0.5         ,7.75025531 );
  graph_fg_ma12_.SetPoint(i++,0.600000024 ,6.55366078 );
  graph_fg_ma12_.SetPoint(i++,0.699999988 ,5.56154077 );
  graph_fg_ma12_.SetPoint(i++,0.800000012 ,4.73954982 );
  graph_fg_ma12_.SetPoint(i++,0.899999976 ,4.05735041 );
  graph_fg_ma12_.SetPoint(i++,1.          ,3.48951745 );
  graph_fg_ma12_.SetPoint(i++,1.10000002  ,3.01535068 );
  graph_fg_ma12_.SetPoint(i++,1.20000005  ,2.53546827 );
  graph_fg_ma12_.SetPoint(i++,1.29999995  ,2.01614944 );
  graph_fg_ma12_.SetPoint(i++,1.39999998  ,1.52225047 );
  graph_fg_ma12_.SetPoint(i++,1.5         ,1.08762294 );
  graph_fg_ma12_.SetPoint(i++,1.60000002  ,0.727223664);
  graph_fg_ma12_.SetPoint(i++,1.70000005  ,0.44478686 );
  graph_fg_ma12_.SetPoint(i++,1.79999995  ,0.237632128);

  adjust_graph(graph_fg_ma10_);
  adjust_graph(graph_fg_ma12_);
  //adjust_graph(graph_sf_ma10_);
  //adjust_graph(graph_sf_ma12_);

  graph_fg_ma10_.SetLineColor(kRed);
  graph_fg_ma10_.SetMarkerColor(kRed);

  graph_fg_ma12_.SetLineColor(kBlack);
  graph_fg_ma12_.SetMarkerColor(kBlack);

  //graph_sf_ma10_.SetLineColor(kBlue);
  //graph_sf_ma10_.SetMarkerColor(kBlue);
  //
  //graph_sf_ma12_.SetLineColor(kMagenta);
  //graph_sf_ma12_.SetMarkerColor(kMagenta);

}
///////////////////////////////////////////////////////////////////
void ButkevichDrawer::draw() {
  gROOT->SetStyle("Plain");

  TH2D null("null","",1,x0_ph_,x1_ph_,1,y0_ph_,y1_ph_);
  null.SetStats(false);
  null.SetTitle(";Q^{2} (GeV^{2});d#sigma/dQ^{2} (cm^{2}/GeV^{2})");
  null.GetXaxis()->CenterTitle();
  null.GetYaxis()->CenterTitle();
  null.GetYaxis()->SetTitleOffset(1.3);

  TCanvas can("can","canvas",600,600);
  null.Draw();

  graph_fg_ma10_.Draw("lpsame");
  graph_fg_ma12_.Draw("lpsame");
  //graph_sf_ma10_.Draw("lpsame");
  //graph_sf_ma12_.Draw("lpsame");

  can.DrawClone();
}
///////////////////////////////////////////////////////////////////
void ButkevichDrawer::adjust_graph(TGraph& graph) {
  for (int iPoint=0;iPoint<graph.GetN();iPoint++) {
    double x,y;
    graph.GetPoint(iPoint,x,y);
    graph.SetPoint(iPoint,x,y*1.0e-38);
  }
}
