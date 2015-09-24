//

// 1. plot

namespace {

  const char   *Fn1 = "~/hist/mu2e/v5_4/egun_stnmaker.track_ana.hist";
  const char   *Fn2 = "~/hist/mu2e/v4_2_6/e00s0160.track_ana.hist";

  int Initialized;

};

//-----------------------------------------------------------------------------
void print_canvas_with_date(TCanvas* C, const char* Name) {
  TDatime d;
  C->Print(Form("%i-%02i-%02i-%s.eps",d.GetYear(),d.GetMonth(),d.GetDay(),Name));
}


//-----------------------------------------------------------------------------
void move_stat_box(TH1F* Hist, double X1, double Y1, double X2, double Y2) {
  TPaveStats *ps = (TPaveStats*) Hist->GetListOfFunctions()->FindObject("stats");

  ps->SetX1NDC(X1); ps->SetY1NDC(Y1);
  ps->SetX2NDC(X2); ps->SetY2NDC(Y2);
}
		   
//-----------------------------------------------------------------------------
// cls_2/energy: cluster energy for events with "Set C" track
//-----------------------------------------------------------------------------
int plot_cls_2_energy() {

  const char* hist_name      = "cls_2/energy";
  const char* plot_file_name = "cls_2_energy";
  const char* title          = "E(cluster)";
  const char* xtitle         = "E, MeV";

  TCanvas* c1 = new TCanvas(Form("c_%s",plot_file_name),plot_file_name,1200,800);

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

  h2->SetTitle(title);
  h2->GetXaxis()->SetTitle(xtitle);
  h2->GetXaxis()->SetRangeUser(0.,150.);
  h2->GetXaxis()->SetTitle(xtitle);
  h2->DrawNormalized("",h1->GetEntries()*h1_bin/h2_bin);

  c1->Update(); // move_stat_box(h2,0.75,0.65,0.98,0.95);
  
  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1);
  h1->SetMarkerColor(kBlue+3);
  h1->Draw("ep,sames");

  c1->Update(); move_stat_box(h1,0.75,0.30,0.98,0.60);

  print_canvas_with_date(c1,plot_file_name);
}


//-----------------------------------------------------------------------------
// evt_6/ncl: N(clusters) for events with Set C tracks
//-----------------------------------------------------------------------------
int plot_evt_6_ncl() {

  const char* hist_name      = "evt_6/ncl";
  const char* plot_file_name = "evt_6_ncl";
  const char* xtitle         = "N(clusters)";

  TCanvas* c1 = new TCanvas(Form("c_%s",plot_file_name),plot_file_name,1200,800);

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

  print_canvas_with_date(c1,plot_file_name);
}

//-----------------------------------------------------------------------------
// evt_6/zv: Z(vertex)
//-----------------------------------------------------------------------------
int plot_evt_6_zv() {

  const char* hist_name      = "evt_6/zv";
  const char* plot_file_name = "evt_6_zv";
  const char* xtitle         = "Z(vertex), mm";

  TCanvas* c1 = new TCanvas(Form("c_%s",plot_file_name),plot_file_name,1200,800);

  TH1F* h1  = gh1(Fn1,"TrackAna",hist_name);
  TH1F* h2  = gh1(Fn2,"TrackAna",hist_name);

  double h1_bin = h1->GetBinWidth(1);
  double h2_bin = h2->GetBinWidth(1);

  printf("%10.4f %10.4f\n",h1->GetBinWidth(1),h2->GetBinWidth(1));

  double y1_max = h1->GetMaximum()/h1->Integral();
  double y2_max = h2->GetMaximum()/h2->Integral()*(h2_bin/h1_bin);

  printf("y1_max, y2_max: %10.4f %10.4f\n",y1_max,y2_max);

  if (y1_max > y2_max) h2->SetMaximum(h2->GetMaximum()*y1_max/y2_max*1.1);

  h2->SetFillColor(2);
  h2->SetFillStyle(3001);
  h2->GetXaxis()->SetRangeUser(5400.,6400.);

  h2->GetXaxis()->SetTitle(xtitle);

  h2->DrawNormalized("",h1->GetEntries()*h1_bin/h2_bin);

  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1);
  h1->SetMarkerColor(kBlue+3);
  h1->Draw("ep,sames");

  print_canvas_with_date(c1,plot_file_name);
}

//-----------------------------------------------------------------------------
int plot_trk_1_p() {

  const char* hist_name      = "trk_1/p";
  const char* plot_file_name = "trk_1_p";
  const char* title          = "Reconstructed Track Momentum";
  const char* xtitle         = "P, MeV/c";

  TCanvas* c1 = new TCanvas(Form("c_%s",plot_file_name),plot_file_name,1200,800);

  TH1F* h1  = gh1(Fn1,"TrackAna",hist_name);
  TH1F* h2  = gh1(Fn2,"TrackAna",hist_name);

  double h1_bin = h1->GetBinWidth(1);
  double h2_bin = h2->GetBinWidth(1);

  printf("%10.4f %10.4f\n",h1->GetBinWidth(1),h2->GetBinWidth(1));

  h2->SetTitle(title);
  h2->SetFillColor(2);
  h2->SetFillStyle(3001);
  h2->GetXaxis()->SetTitle(xtitle);

  h2->GetXaxis()->SetRangeUser(90.,110.);

  h2->DrawNormalized("",h1->GetEntries()*h1_bin/h2_bin);

  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1);
  h1->SetMarkerColor(kBlue+3);
  h1->Draw("ep,sames");

  print_canvas_with_date(c1,plot_file_name);
}


//-----------------------------------------------------------------------------
int plot_trk_1_ep() {

  const char* hist_name      = "trk_1/ep";
  const char* plot_file_name = "trk_1_ep";
  const char* title          = "E(cluster)/P(track)";
  const char* xtitle         = "E/P";

  TCanvas* c1 = new TCanvas(Form("c_%s",plot_file_name),plot_file_name,1200,800);

  TH1F* h1  = gh1(Fn1,"TrackAna",hist_name);
  TH1F* h2  = gh1(Fn2,"TrackAna",hist_name);

  double h1_bin = h1->GetBinWidth(1);
  double h2_bin = h2->GetBinWidth(1);

  printf("%10.4f %10.4f\n",h1->GetBinWidth(1),h2->GetBinWidth(1));

  double y1_max = h1->GetMaximum()/h1->GetEntries();
  double y2_max = h2->GetMaximum()/h2->GetEntries()*(h2_bin/h1_bin);

  if (y1_max > y2_max) h2->SetMaximum(h2->GetMaximum()*y1_max/y2_max*1.1);

  h2->SetTitle(title);
  h2->GetXaxis()->SetTitle(xtitle);

  h2->GetXaxis()->SetRangeUser(0.,1.2);

  h2->SetFillColor(2);
  h2->SetFillStyle(3001);
  h2->DrawNormalized("",h1->GetEntries()*h1_bin/h2_bin);

  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1);
  h1->SetMarkerColor(kBlue+3);
  h1->Draw("ep,sames");

  print_canvas_with_date(c1,plot_file_name);
}

//-----------------------------------------------------------------------------
// cls_2/energy: cluster energy for events with "Set C" track
//-----------------------------------------------------------------------------
int plot_trk_1_dx() {

  const char* hist_name      = "trk_1/dx";
  const char* plot_file_name = "trk_1_dx";
  const char* title          = "#Delta X = X_{track} - X_{cluster}";
  const char* xtitle         = "#Delta X, mm";

  TCanvas* c1 = new TCanvas(Form("c_%s",plot_file_name),plot_file_name,1200,800);

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

  h2->SetTitle(title);
  h2->GetXaxis()->SetTitle(xtitle);
  // h2->GetXaxis()->SetRangeUser(10.,160.);

  h2->DrawNormalized("",h1->GetEntries()*h1_bin/h2_bin);

  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1);
  h1->SetMarkerColor(kBlue+3);
  h1->Draw("ep,sames");

  c1->Update(); move_stat_box(h1,0.75,0.30,0.98,0.60);

  print_canvas_with_date(c1,plot_file_name);
}


//-----------------------------------------------------------------------------
// trk_1/dy: 
//-----------------------------------------------------------------------------
int plot_trk_1_dy() {

  const char* hist_name      = "trk_1/dy";
  const char* plot_file_name = "trk_1_dy";
  const char* title          = "#Delta Y = Y_{track} - Y_{cluster}";
  const char* xtitle         = "#Delta Y, mm";

  TCanvas* c1 = new TCanvas(Form("c_%s",plot_file_name),plot_file_name,1200,800);

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

  h2->SetTitle(title);
  h2->GetXaxis()->SetTitle(xtitle);
  //  h2->GetXaxis()->SetRangeUser(10.,160.);
  h2->DrawNormalized("",h1->GetEntries()*h1_bin/h2_bin);

  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1);
  h1->SetMarkerColor(kBlue+3);
  h1->Draw("ep,sames");

  c1->Update(); move_stat_box(h1,0.75,0.30,0.98,0.60);

  print_canvas_with_date(c1,plot_file_name);
}



//-----------------------------------------------------------------------------
// trk_1/dy: 
//-----------------------------------------------------------------------------
int plot_trk_1_dt() {

  const char* hist_name      = "trk_1/dt";
  const char* plot_file_name = "trk_1_dt";
  const char* xtitle         = "Track-Cluster #Delta{Y}, mm";

  TCanvas* c1 = new TCanvas(Form("c_%s",plot_file_name),plot_file_name,1200,800);
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

  h2->GetXaxis()->SetTitle(xtitle);
  h2->DrawNormalized("",h1->GetEntries()*h1_bin/h2_bin);

  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1);
  h1->SetMarkerColor(kBlue+3);
  h1->Draw("ep,sames");

  print_canvas_with_date(c1,plot_file_name);
}

//-----------------------------------------------------------------------------
// trk_1/nactv: 
//-----------------------------------------------------------------------------
int plot_trk_1_nactv() {

  const char* hist_name      = "trk_1/nactv";
  const char* plot_file_name = "trk_1_nactv";
  const char* title          = "Number of active hits (hits used by Kalman Fit)";
  const char* xtitle         = "N(active hits)";

  TCanvas* c1 = new TCanvas(Form("c_%s",plot_file_name),plot_file_name,1200,800);

  TH1F* h1  = gh1(Fn1,"TrackAna",hist_name);
  TH1F* h2  = gh1(Fn2,"TrackAna",hist_name);

  double h1_bin = h1->GetBinWidth(1);
  double h2_bin = h2->GetBinWidth(1);

  printf("%10.4f %10.4f\n",h1->GetBinWidth(1),h2->GetBinWidth(1));

  double y1_max = h1->GetMaximum()/h1->Integral();
  double y2_max = h2->GetMaximum()/h2->Integral()*(h2_bin/h1_bin);

  printf("y1_max, y2_max: %10.4f %10.4f\n",y1_max,y2_max);

  if (y1_max > y2_max) h2->SetMaximum(h2->GetMaximum()*y1_max/y2_max*1.1);

  h2->SetFillColor(2);
  h2->SetFillStyle(3001);

  h2->SetTitle(title);
  h2->GetXaxis()->SetRangeUser(0.,120.);
  h2->GetXaxis()->SetTitle(xtitle);
  h2->DrawNormalized("",h1->GetEntries()*h1_bin/h2_bin);

  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1);
  h1->SetMarkerColor(kBlue+3);
  h1->Draw("ep,sames");

  c1->Update();
  move_stat_box(h1,0.75,0.30,0.98,0.60);

  print_canvas_with_date(c1,plot_file_name);
}



//-----------------------------------------------------------------------------
// trk_1/nwrng: 
//-----------------------------------------------------------------------------
int plot_trk_1_nwrng() {

  const char* hist_name      = "trk_1/nwrng";
  const char* plot_file_name = "trk_1_nwrng";
  const char* xtitle         = "N(non-active hits)";

  TCanvas* c1 = new TCanvas(Form("c_%s",plot_file_name),plot_file_name,1200,800);

  TH1F* h1  = gh1(Fn1,"TrackAna",hist_name);
  TH1F* h2  = gh1(Fn2,"TrackAna",hist_name);

  double h1_bin = h1->GetBinWidth(1);
  double h2_bin = h2->GetBinWidth(1);

  printf("%10.4f %10.4f\n",h1->GetBinWidth(1),h2->GetBinWidth(1));

  double y1_max = h1->GetMaximum()/h1->Integral();
  double y2_max = h2->GetMaximum()/h2->Integral()*(h2_bin/h1_bin);

  if (y1_max > y2_max) h2->SetMaximum(h2->GetMaximum()*y1_max/y2_max*1.1);

  h2->SetFillColor(2);
  h2->SetFillStyle(3001);
  h2->GetXaxis()->SetTitle(xtitle);
  //  h2->GetXaxis()->SetRangeUser(0.,120.);
  h2->DrawNormalized("",h1->GetEntries()*h1_bin/h2_bin);

  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1);
  h1->SetMarkerColor(kBlue+3);
  h1->Draw("ep,sames");

  print_canvas_with_date(c1,plot_file_name);
}

//-----------------------------------------------------------------------------
// trk_1/nactv: 
//-----------------------------------------------------------------------------
int plot_trk_1_momerr() {

  const char* hist_name      = "trk_1/momerr";
  const char* plot_file_name = "trk_1_momerr";
  const char* title          = "track momentum error";
  const char* xtitle         = "#sigma_{P}, MeV/c";

  TCanvas* c1 = new TCanvas(Form("c_%s",plot_file_name),plot_file_name,1200,800);

  TH1F* h1  = gh1(Fn1,"TrackAna",hist_name);
  TH1F* h2  = gh1(Fn2,"TrackAna",hist_name);

  double h1_bin = h1->GetBinWidth(1);
  double h2_bin = h2->GetBinWidth(1);

  printf("%10.4f %10.4f\n",h1->GetBinWidth(1),h2->GetBinWidth(1));

  double y1_max = h1->GetMaximum()/h1->Integral();
  double y2_max = h2->GetMaximum()/h2->Integral()*(h2_bin/h1_bin);

  if (y1_max > y2_max) h2->SetMaximum(h2->GetMaximum()*y1_max/y2_max*1.1);

  h2->SetFillColor(2);
  h2->SetFillStyle(3001);

  h2->SetTitle(title);
  h2->GetXaxis()->SetTitle(xtitle);
  h2->GetXaxis()->SetRangeUser(0.,0.4);

  h2->DrawNormalized("",h1->GetEntries()*h1_bin/h2_bin);

  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1);
  h1->SetMarkerColor(kBlue+3);
  h1->Draw("ep,sames");

  c1->Update();

  move_stat_box(h1,0.75,0.30,0.98,0.60);
  
  c1->Modified();
  c1->Update();
  print_canvas_with_date(c1,plot_file_name);
}

//-----------------------------------------------------------------------------
// trk_1/dpf: 
//-----------------------------------------------------------------------------
int plot_trk_1_dpf() {

  const char* hist_name      = "trk_1/dpf";
  const char* plot_file_name = "trk_1_dpf";
  const char* xtitle         = "P_{rec} -P(tracker front), MeV/c";

  TCanvas* c1 = new TCanvas(Form("c_%s",plot_file_name),plot_file_name,1200,800);

  TH1F* h1  = gh1(Fn1,"TrackAna",hist_name);
  TH1F* h2  = gh1(Fn2,"TrackAna",hist_name);

  double h1_bin = h1->GetBinWidth(1);
  double h2_bin = h2->GetBinWidth(1);

  printf("%10.4f %10.4f\n",h1->GetBinWidth(1),h2->GetBinWidth(1));

  double y1_max = h1->GetMaximum()/h1->Integral();
  double y2_max = h2->GetMaximum()/h2->Integral()*(h2_bin/h1_bin);

  if (y1_max > y2_max) h2->SetMaximum(h2->GetMaximum()*y1_max/y2_max*1.1);

  h2->SetFillColor(2);
  h2->SetFillStyle(3001);

  h2->GetXaxis()->SetRangeUser(-2.,2.);

  h2->GetXaxis()->SetTitle(xtitle);
  h2->DrawNormalized("",h1->GetEntries()*h1_bin/h2_bin);

  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1);
  h1->SetMarkerColor(kBlue+3);
  h1->Draw("ep,sames");

  print_canvas_with_date(c1,plot_file_name);
}

//-----------------------------------------------------------------------------
// trk_1/dpfst: 
//-----------------------------------------------------------------------------
int plot_trk_1_pstout() {

  const char* hist_name      = "trk_1/pstout";
  const char* plot_file_name = "trk_1_pstout";
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

  print_canvas_with_date(c1,plot_file_name);
}

//-----------------------------------------------------------------------------
// trk_1/dpfst: 
//-----------------------------------------------------------------------------
int plot_trk_1_dpfst() {

  const char* hist_name      = "trk_1/dpfst";
  const char* plot_file_name = "trk_1_dpfst";
  const char* xtitle         = "P_{1} -P_{prod}, MeV/c";

  TH1F* h1  = gh1(Fn1,"TrackAna",hist_name);
  TH1F* h2  = gh1(Fn2,"TrackAna",hist_name);

  double h1_bin = h1->GetBinWidth(1);
  double h2_bin = h2->GetBinWidth(1);

  printf("%10.4f %10.4f\n",h1->GetBinWidth(1),h2->GetBinWidth(1));

  h2->SetFillColor(2);
  h2->SetFillStyle(3001);

  h2->GetXaxis()->SetRangeUser(-5,1.);

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

  print_canvas_with_date(c1,plot_file_name);
}

//-----------------------------------------------------------------------------
plot_validation_histograms() {
  plot_trk_1_p ();
  plot_trk_1_momerr();
  plot_trk_1_ep();
  plot_trk_1_nactv();

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
  plot_cls_2_energy();
  plot_evt_6_ncl();
  plot_evt_6_zv();
}
