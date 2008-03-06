{
  gROOT->ProcessLine(".L look.C+");

  TString ps_file="../val.ps";
  TString options="";

  //single_compare(50,"../../events/events_val0.root");
  //gPad->Print(ps_file+"(",options);
  //
  polish_compare(1,50,"../../events/events_val1.root");
  gPad->Print(ps_file+"(",options);

  polish_compare(1,60,"../../events/events_val1.root");
  gPad->Print(ps_file,options);

  polish_compare(2,60,"../../events/events_val2.root");
  gPad->Print(ps_file,options);

  polish_compare(3,60,"../../events/events_val3.root");
  gPad->Print(ps_file,options);

  polish_compare(23,60,"../../events/events_val4.root");
  gPad->Print(ps_file,options);

  NUANCE_compare(1,800,"../../events/events_val5.root");
  gPad->Print(ps_file,options);

  NUANCE_compare(2,800,"../../events/events_val5.root");
  gPad->Print(ps_file,options);

  NUANCE_compare(1,300,"../../events/events_val6.root");
  gPad->Print(ps_file,options);

  NUANCE_compare(2,300,"../../events/events_val6.root");
  gPad->Print(ps_file,options);

  gPad->Print(ps_file+"]",options);
  gROOT->ProcessLine(".! gv "+ps_file+" &");
}
