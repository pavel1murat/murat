//

// 1. plot

namespace {

  const char   *Fn1 = "hist/egun_stnmaker.track_ana.hist";
  const char   *Fn2 = "hist/../v4_2_6/e00s0160.track_ana.hist";

  int Initialized;

};



//-----------------------------------------------------------------------------
int plot_trk_1_p() {

  TH1F* h1  = gh1(Fn1,"TrackAna","trk_1/p");
  TH1F* h2  = gh1(Fn2,"TrackAna","trk_1/p");

  double h1_bin = h1->GetBinWidth(1);
  double h2_bin = h2->GetBinWidth(1);

  printf("%10.4f %10.4f\n",h1->GetBinWidth(1),h2->GetBinWidth(1));

  h2->SetFillColor(2);
  h2->SetFillStyle(3001);

  h2->GetXaxis()->SetRangeUser(90.,110.);

  h2->DrawNormalized("",h1->GetEntries()*h1_bin/h2_bin);

  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1);
  h1->SetMarkerColor(kBlue+3);
  h1->Draw("ep,sames");


  c1->Print("trk_1_p.eps");
}


//-----------------------------------------------------------------------------
int plot_trk_1_ep() {

  TH1F* h1  = gh1(Fn1,"TrackAna","trk_1/ep");
  TH1F* h2  = gh1(Fn2,"TrackAna","trk_1/ep");

  double h1_bin = h1->GetBinWidth(1);
  double h2_bin = h2->GetBinWidth(1);

  printf("%10.4f %10.4f\n",h1->GetBinWidth(1),h2->GetBinWidth(1));

  h2->SetFillColor(2);
  h2->SetFillStyle(3001);

  //  h2->GetXaxis()->SetRangeUser(90.,110.);

  double y1_max = h1->GetMaximum()/h1->GetEntries();
  double y2_max = h2->GetMaximum()/h2->GetEntries()*(h2_bin/h1_bin);

  if (y1_max > y2_max) h2->SetMaximum(h2->GetMaximum()*y1_max/y2_max*1.1);

  h2->DrawNormalized("",h1->GetEntries()*h1_bin/h2_bin);

  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1);
  h1->SetMarkerColor(kBlue+3);
  h1->Draw("ep,sames");


  c1->Print("trk_1_ep.eps");
}

//-----------------------------------------------------------------------------
// cls_2/energy: cluster energy for events with "Set C" track
//-----------------------------------------------------------------------------
int plot_cls_2_energy() {

  const char* hist_name      = "cls_2/energy";
  const char* plot_file_name = "cls_2_energy.eps";
  const char* xtitle         = "E(cluster)";

  TH1F* h1  = gh1(Fn1,"TrackAna",hist_name);
  TH1F* h2  = gh1(Fn2,"TrackAna",hist_name);

  double h1_bin = h1->GetBinWidth(1);
  double h2_bin = h2->GetBinWidth(1);

  printf("%10.4f %10.4f\n",h1->GetBinWidth(1),h2->GetBinWidth(1));

  h2->SetFillColor(2);
  h2->SetFillStyle(3001);

  double y1_max = h1->GetMaximum()/h1->Integral();
  double y2_max = h2->GetMaximum()/h2->Integral()*(h2_bin/h1_bin);

  if (y1_max > y2_max) h2->SetMaximum(h2->GetMaximum()*y1_max/y2_max*1.1);

  h2->GetXaxis()->SetTitle(xtitle);
  h2->GetXaxis()->SetRangeUser(0.,150.);
  h2->DrawNormalized("",h1->GetEntries()*h1_bin/h2_bin);

  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1);
  h1->SetMarkerColor(kBlue+3);
  h1->Draw("ep,sames");


  c1->Print("cls_2_energy.eps");
}


//-----------------------------------------------------------------------------
// cls_2/energy: cluster energy for events with "Set C" track
//-----------------------------------------------------------------------------
int plot_trk_1_dx() {

  const char* hist_name      = "trk_1/dx";
  const char* plot_file_name = "trk_1_dx.eps";

  TH1F* h1  = gh1(Fn1,"TrackAna",hist_name);
  TH1F* h2  = gh1(Fn2,"TrackAna",hist_name);

  double h1_bin = h1->GetBinWidth(1);
  double h2_bin = h2->GetBinWidth(1);

  printf("%10.4f %10.4f\n",h1->GetBinWidth(1),h2->GetBinWidth(1));

  h2->SetFillColor(2);
  h2->SetFillStyle(3001);

  //  h2->GetXaxis()->SetRangeUser(10.,160.);

  double y1_max = h1->GetMaximum()/h1->Integral();
  double y2_max = h2->GetMaximum()/h2->Integral()*(h2_bin/h1_bin);

  if (y1_max > y2_max) h2->SetMaximum(h2->GetMaximum()*y1_max/y2_max*1.1);

  h2->GetXaxis()->SetTitle("DX, mm");
  h2->DrawNormalized("",h1->GetEntries()*h1_bin/h2_bin);

  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1);
  h1->SetMarkerColor(kBlue+3);
  h1->Draw("ep,sames");


  c1->Print(plot_file_name);
}


//-----------------------------------------------------------------------------
// trk_1/dy: 
//-----------------------------------------------------------------------------
int plot_trk_1_dy() {

  const char* hist_name      = "trk_1/dy";
  const char* plot_file_name = "trk_1_dy.eps";

  TH1F* h1  = gh1(Fn1,"TrackAna",hist_name);
  TH1F* h2  = gh1(Fn2,"TrackAna",hist_name);

  double h1_bin = h1->GetBinWidth(1);
  double h2_bin = h2->GetBinWidth(1);

  printf("%10.4f %10.4f\n",h1->GetBinWidth(1),h2->GetBinWidth(1));

  h2->SetFillColor(2);
  h2->SetFillStyle(3001);

  //  h2->GetXaxis()->SetRangeUser(10.,160.);

  double y1_max = h1->GetMaximum()/h1->Integral();
  double y2_max = h2->GetMaximum()/h2->Integral()*(h2_bin/h1_bin);

  if (y1_max > y2_max) h2->SetMaximum(h2->GetMaximum()*y1_max/y2_max*1.1);

  h2->GetXaxis()->SetTitle("DY, mm");
  h2->DrawNormalized("",h1->GetEntries()*h1_bin/h2_bin);

  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1);
  h1->SetMarkerColor(kBlue+3);
  h1->Draw("ep,sames");


  c1->Print(plot_file_name);
}



//-----------------------------------------------------------------------------
// trk_1/dy: 
//-----------------------------------------------------------------------------
int plot_trk_1_dt() {

  const char* hist_name      = "trk_1/dt";
  const char* plot_file_name = "trk_1_dt.eps";

  TH1F* h1  = gh1(Fn1,"TrackAna",hist_name);
  TH1F* h2  = gh1(Fn2,"TrackAna",hist_name);

  double h1_bin = h1->GetBinWidth(1);
  double h2_bin = h2->GetBinWidth(1);

  printf("%10.4f %10.4f\n",h1->GetBinWidth(1),h2->GetBinWidth(1));

  h2->SetFillColor(2);
  h2->SetFillStyle(3001);

  h2->GetXaxis()->SetRangeUser(-10.,10.);

  double y1_max = h1->GetMaximum()/h1->Integral();
  double y2_max = h2->GetMaximum()/h2->Integral()*(h2_bin/h1_bin);

  if (y1_max > y2_max) h2->SetMaximum(h2->GetMaximum()*y1_max/y2_max*1.1);

  h2->GetXaxis()->SetTitle("DT, ns");
  h2->DrawNormalized("",h1->GetEntries()*h1_bin/h2_bin);

  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1);
  h1->SetMarkerColor(kBlue+3);
  h1->Draw("ep,sames");


  c1->Print(plot_file_name);
}

//-----------------------------------------------------------------------------
// trk_1/nactv: 
//-----------------------------------------------------------------------------
int plot_trk_1_nactv() {

  const char* hist_name      = "trk_1/nactv";
  const char* plot_file_name = "trk_1_nactv.eps";
  const char* xtitle         = "N(active hits)";

  TH1F* h1  = gh1(Fn1,"TrackAna",hist_name);
  TH1F* h2  = gh1(Fn2,"TrackAna",hist_name);

  double h1_bin = h1->GetBinWidth(1);
  double h2_bin = h2->GetBinWidth(1);

  printf("%10.4f %10.4f\n",h1->GetBinWidth(1),h2->GetBinWidth(1));

  h2->SetFillColor(2);
  h2->SetFillStyle(3001);

  h2->GetXaxis()->SetRangeUser(0.,120.);

  double y1_max = h1->GetMaximum()/h1->Integral();
  double y2_max = h2->GetMaximum()/h2->Integral()*(h2_bin/h1_bin);

  printf("y1_max, y2_max: %10.4f %10.4f\n",y1_max,y2_max);

  if (y1_max > y2_max) h2->SetMaximum(h2->GetMaximum()*y1_max/y2_max*1.1);

  h2->GetXaxis()->SetTitle(xtitle);

  h2->DrawNormalized("",h1->GetEntries()*h1_bin/h2_bin);

  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1);
  h1->SetMarkerColor(kBlue+3);
  h1->Draw("ep,sames");

  c1->Print(plot_file_name);
}



//-----------------------------------------------------------------------------
// trk_1/nwrng: 
//-----------------------------------------------------------------------------
int plot_trk_1_nwrng() {

  const char* hist_name      = "trk_1/nwrng";
  const char* plot_file_name = "trk_1_nwrng.eps";
  const char* xtitle         = "N(non-active hits)";

  TH1F* h1  = gh1(Fn1,"TrackAna",hist_name);
  TH1F* h2  = gh1(Fn2,"TrackAna",hist_name);

  double h1_bin = h1->GetBinWidth(1);
  double h2_bin = h2->GetBinWidth(1);

  printf("%10.4f %10.4f\n",h1->GetBinWidth(1),h2->GetBinWidth(1));

  h2->SetFillColor(2);
  h2->SetFillStyle(3001);


  double y1_max = h1->GetMaximum()/h1->Integral();
  double y2_max = h2->GetMaximum()/h2->Integral()*(h2_bin/h1_bin);

  printf("y1_max, y2_max: %10.4f %10.4f\n",y1_max,y2_max);

  if (y1_max > y2_max) h2->SetMaximum(h2->GetMaximum()*y1_max/y2_max*1.1);

  h2->GetXaxis()->SetTitle(xtitle);
  //  h2->GetXaxis()->SetRangeUser(0.,120.);
  h2->DrawNormalized("",h1->GetEntries()*h1_bin/h2_bin);

  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1);
  h1->SetMarkerColor(kBlue+3);
  h1->Draw("ep,sames");

  c1->Print(plot_file_name);
}



//-----------------------------------------------------------------------------
// trk_1/nactv: 
//-----------------------------------------------------------------------------
int plot_trk_1_momerr() {

  const char* hist_name      = "trk_1/momerr";
  const char* plot_file_name = "trk_1_momerr.eps";
  const char* xtitle         = "#sigma_{P}, MeV/c";

  TH1F* h1  = gh1(Fn1,"TrackAna",hist_name);
  TH1F* h2  = gh1(Fn2,"TrackAna",hist_name);

  double h1_bin = h1->GetBinWidth(1);
  double h2_bin = h2->GetBinWidth(1);

  printf("%10.4f %10.4f\n",h1->GetBinWidth(1),h2->GetBinWidth(1));

  h2->SetFillColor(2);
  h2->SetFillStyle(3001);

  //  h2->GetXaxis()->SetRangeUser(0.,120.);

  double y1_max = h1->GetMaximum()/h1->Integral();
  double y2_max = h2->GetMaximum()/h2->Integral()*(h2_bin/h1_bin);

  printf("y1_max, y2_max: %10.4f %10.4f\n",y1_max,y2_max);

  if (y1_max > y2_max) h2->SetMaximum(h2->GetMaximum()*y1_max/y2_max*1.1);

  h2->GetXaxis()->SetTitle(xtitle);

  h2->DrawNormalized("",h1->GetEntries()*h1_bin/h2_bin);

  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1);
  h1->SetMarkerColor(kBlue+3);
  h1->Draw("ep,sames");


  c1->Print(plot_file_name);
}

//-----------------------------------------------------------------------------
// trk_1/dpf: 
//-----------------------------------------------------------------------------
int plot_trk_1_dpf() {

  const char* hist_name      = "trk_1/dpf";
  const char* plot_file_name = "trk_1_dpf.eps";
  const char* xtitle         = "P_{rec} -P(tracker front), MeV/c";

  TH1F* h1  = gh1(Fn1,"TrackAna",hist_name);
  TH1F* h2  = gh1(Fn2,"TrackAna",hist_name);

  double h1_bin = h1->GetBinWidth(1);
  double h2_bin = h2->GetBinWidth(1);

  printf("%10.4f %10.4f\n",h1->GetBinWidth(1),h2->GetBinWidth(1));

  h2->SetFillColor(2);
  h2->SetFillStyle(3001);

  h2->GetXaxis()->SetRangeUser(-2.,2.);

  double y1_max = h1->GetMaximum()/h1->Integral();
  double y2_max = h2->GetMaximum()/h2->Integral()*(h2_bin/h1_bin);

  printf("y1_max, y2_max: %10.4f %10.4f\n",y1_max,y2_max);

  if (y1_max > y2_max) h2->SetMaximum(h2->GetMaximum()*y1_max/y2_max*1.1);

  h2->GetXaxis()->SetTitle(xtitle);

  h2->DrawNormalized("",h1->GetEntries()*h1_bin/h2_bin);

  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1);
  h1->SetMarkerColor(kBlue+3);
  h1->Draw("ep,sames");


  c1->Print(plot_file_name);
}

//-----------------------------------------------------------------------------
// trk_1/dpfst: 
//-----------------------------------------------------------------------------
int plot_trk_1_pstout() {

  const char* hist_name      = "trk_1/pstout";
  const char* plot_file_name = "trk_1_pstout.eps";
  const char* xtitle         = "P(End of the Stopping Target), MeV/c";

  TH1F* h1  = gh1(Fn1,"TrackAna",hist_name);
  TH1F* h2  = gh1(Fn2,"TrackAna",hist_name);

  double h1_bin = h1->GetBinWidth(1);
  double h2_bin = h2->GetBinWidth(1);

  printf("%10.4f %10.4f\n",h1->GetBinWidth(1),h2->GetBinWidth(1));

  h2->SetFillColor(2);
  h2->SetFillStyle(3001);

  double y1_max = h1->GetMaximum()/h1->Integral();
  double y2_max = h2->GetMaximum()/h2->Integral()*(h2_bin/h1_bin);

  printf("y1_max, y2_max: %10.4f %10.4f\n",y1_max,y2_max);

  if (y1_max > y2_max) h2->SetMaximum(h2->GetMaximum()*y1_max/y2_max*1.1);

  h2->GetXaxis()->SetTitle(xtitle);
  //  h2->GetXaxis()->SetRangeUser(0.,120.);
  h2->DrawNormalized("",h1->GetEntries()*h1_bin/h2_bin);

  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1);
  h1->SetMarkerColor(kBlue+3);
  h1->Draw("ep,sames");


  c1->Print(plot_file_name);
}

//-----------------------------------------------------------------------------
// trk_1/dpfst: 
//-----------------------------------------------------------------------------
int plot_trk_1_dpfst() {

  const char* hist_name      = "trk_1/dpfst";
  const char* plot_file_name = "trk_1_dpfst.eps";
  const char* xtitle         = "P_{1} -P_{prod}, MeV/c";

  TH1F* h1  = gh1(Fn1,"TrackAna",hist_name);
  TH1F* h2  = gh1(Fn2,"TrackAna",hist_name);

  double h1_bin = h1->GetBinWidth(1);
  double h2_bin = h2->GetBinWidth(1);

  printf("%10.4f %10.4f\n",h1->GetBinWidth(1),h2->GetBinWidth(1));

  h2->SetFillColor(2);
  h2->SetFillStyle(3001);

  //  h2->GetXaxis()->SetRangeUser(0.,120.);

  double y1_max = h1->GetMaximum()/h1->Integral();
  double y2_max = h2->GetMaximum()/h2->Integral()*(h2_bin/h1_bin);

  printf("y1_max, y2_max: %10.4f %10.4f\n",y1_max,y2_max);

  if (y1_max > y2_max) h2->SetMaximum(h2->GetMaximum()*y1_max/y2_max*1.1);

  h2->GetXaxis()->SetTitle(xtitle);

  h2->DrawNormalized("",h1->GetEntries()*h1_bin/h2_bin);

  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1);
  h1->SetMarkerColor(kBlue+3);
  h1->Draw("ep,sames");


  c1->Print(plot_file_name);
}


//-----------------------------------------------------------------------------
// evt_6/ncl: N(clusters) for events with Set C tracks
//-----------------------------------------------------------------------------
int plot_evt_6_ncl() {

  const char* hist_name      = "evt_6/ncl";
  const char* plot_file_name = "evt_6_ncl.eps";
  const char* xtitle         = "N(clusters)";

  TH1F* h1  = gh1(Fn1,"TrackAna",hist_name);
  TH1F* h2  = gh1(Fn2,"TrackAna",hist_name);

  double h1_bin = h1->GetBinWidth(1);
  double h2_bin = h2->GetBinWidth(1);

  printf("%10.4f %10.4f\n",h1->GetBinWidth(1),h2->GetBinWidth(1));

  h2->SetFillColor(2);
  h2->SetFillStyle(3001);

  h2->GetXaxis()->SetRangeUser(0.,10.);

  double y1_max = h1->GetMaximum()/h1->Integral();
  double y2_max = h2->GetMaximum()/h2->Integral()*(h2_bin/h1_bin);

  printf("y1_max, y2_max: %10.4f %10.4f\n",y1_max,y2_max);

  if (y1_max > y2_max) h2->SetMaximum(h2->GetMaximum()*y1_max/y2_max*1.1);

  h2->GetXaxis()->SetTitle(xtitle);

  h2->DrawNormalized("",h1->GetEntries()*h1_bin/h2_bin);

  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1);
  h1->SetMarkerColor(kBlue+3);
  h1->Draw("ep,sames");

  c1->Print(plot_file_name);
}



//-----------------------------------------------------------------------------
// evt_6/zv: Z(vertex)
//-----------------------------------------------------------------------------
int plot_evt_6_zv() {

  const char* hist_name      = "evt_6/zv";
  const char* plot_file_name = "evt_6_zv.eps";
  const char* xtitle         = "Z(vertex), mm";

  TH1F* h1  = gh1(Fn1,"TrackAna",hist_name);
  TH1F* h2  = gh1(Fn2,"TrackAna",hist_name);

  double h1_bin = h1->GetBinWidth(1);
  double h2_bin = h2->GetBinWidth(1);

  printf("%10.4f %10.4f\n",h1->GetBinWidth(1),h2->GetBinWidth(1));

  h2->SetFillColor(2);
  h2->SetFillStyle(3001);

  //  h2->GetXaxis()->SetRangeUser(0.,10.);

  double y1_max = h1->GetMaximum()/h1->Integral();
  double y2_max = h2->GetMaximum()/h2->Integral()*(h2_bin/h1_bin);

  printf("y1_max, y2_max: %10.4f %10.4f\n",y1_max,y2_max);

  if (y1_max > y2_max) h2->SetMaximum(h2->GetMaximum()*y1_max/y2_max*1.1);

  h2->GetXaxis()->SetTitle(xtitle);

  h2->DrawNormalized("",h1->GetEntries()*h1_bin/h2_bin);

  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1);
  h1->SetMarkerColor(kBlue+3);
  h1->Draw("ep,sames");

  c1->Print(plot_file_name);
}

//-----------------------------------------------------------------------------
// evt_6/ecal: Z(vertex)
//-----------------------------------------------------------------------------
int plot_evt_6_ecal() {

  const char* hist_name      = "evt_6/ecal";
  const char* plot_file_name = "evt_6_ecal.eps";
  const char* xtitle         = "E(cal), MeV";

  TH1F* h1  = gh1(Fn1,"TrackAna",hist_name);
  TH1F* h2  = gh1(Fn2,"TrackAna",hist_name);

  double h1_bin = h1->GetBinWidth(1);
  double h2_bin = h2->GetBinWidth(1);

  printf("%10.4f %10.4f\n",h1->GetBinWidth(1),h2->GetBinWidth(1));

  h2->SetFillColor(2);
  h2->SetFillStyle(3001);

  //  h2->GetXaxis()->SetRangeUser(0.,10.);

  double y1_max = h1->GetMaximum()/h1->Integral();
  double y2_max = h2->GetMaximum()/h2->Integral()*(h2_bin/h1_bin);

  printf("y1_max, y2_max: %10.4f %10.4f\n",y1_max,y2_max);

  if (y1_max > y2_max) h2->SetMaximum(h2->GetMaximum()*y1_max/y2_max*1.1);

  h2->GetXaxis()->SetTitle(xtitle);

  h2->DrawNormalized("",h1->GetEntries()*h1_bin/h2_bin);

  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1);
  h1->SetMarkerColor(kBlue+3);
  h1->Draw("ep,sames");

  c1->Print(plot_file_name);
}





//-----------------------------------------------------------------------------
plot_validation_histograms() {
  plot_trk_1_p ();
  plot_trk_1_momerr();

  plot_trk_1_ep();

  plot_trk_1_dpf();
  plot_trk_1_dpfst();
					// MC: momentum on STOUT virtual detector
  plot_trk_1_pstout();
//-----------------------------------------------------------------------------
// track-to-cluster residuals
//-----------------------------------------------------------------------------
  plot_trk_1_dx();
  plot_trk_1_dy();
  plot_trk_1_dt();

//-----------------------------------------------------------------------------
// clusters in events with good (now - Set C) tracks
//-----------------------------------------------------------------------------
  plot_evt_6_ncl();
  plot_cls_2_energy();
  plot_evt_6_zv();
}
