///////////////////////////////////////////////////////////////////////////////
// compare histograms in file 'Fn' against reference histograms in file 'FnRef'
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
void print_canvas_with_date(TCanvas* C, const char* Name) {
  TDatime d;
  C->Print(Form("%i-%02i-%02i-%s.eps",d.GetYear(),d.GetMonth(),d.GetDay(),Name));
}

//-----------------------------------------------------------------------------
void plot_validation_c01(const char* FnRef, const char* Fn) {
  const char* oname="plot_validation_c01";

  TCanvas* c1 = new TCanvas("c01","c01",1200,900);
  c1->Divide(2,2);
//-----------------------------------------------------------------------------
// track momentum
//-----------------------------------------------------------------------------
  c1->cd(1);

  TH1* href = gh1(FnRef,"TrackAna","trk_1/p_2");
  TH1* h    = gh1(Fn   ,"TrackAna","trk_1/p_2");

  href->SetLineColor(1);
  href->SetFillColor(623);
  href->SetFillStyle(3001);

  href->GetXaxis()->SetRangeUser(70,120);

  h->SetMarkerStyle(20);
  h->SetMarkerSize (1.);
  h->SetMarkerColor(kBlue+3);

  href->Draw();
  h->Draw("sames,ep");
//-----------------------------------------------------------------------------
// npoints
//-----------------------------------------------------------------------------
  c1->cd(2);

  TH1* href = gh1(FnRef,"TrackAna","trk_1/nactv");
  TH1* h    = gh1(Fn   ,"TrackAna","trk_1/nactv");

  href->SetLineColor(1);
  href->SetFillColor(623);
  href->SetFillStyle(3001);
  href->GetXaxis()->SetRangeUser(0,119.9);

  h->SetMarkerStyle(20);
  h->SetMarkerSize (1.);
  h->SetMarkerColor(kBlue+3);

  href->Draw();
  h->Draw("sames,ep");
//-----------------------------------------------------------------------------
// chi2
//-----------------------------------------------------------------------------
  c1->cd(3);

  TH1* href = gh1(FnRef,"TrackAna","trk_1/chi2d");
  TH1* h    = gh1(Fn   ,"TrackAna","trk_1/chi2d");

  href->SetLineColor(1);
  href->SetFillColor(623);
  href->SetFillStyle(3001);
  href->GetXaxis()->SetRangeUser(0,4.999);

  h->SetMarkerStyle(20);
  h->SetMarkerSize (1.);
  h->SetMarkerColor(kBlue+3);

  href->Draw();
  h->Draw("sames,ep");
//-----------------------------------------------------------------------------
// tan(dip)
//-----------------------------------------------------------------------------
  c1->cd(4);

  TH1* href = gh1(FnRef,"TrackAna","trk_1/tdip");
  TH1* h    = gh1(Fn   ,"TrackAna","trk_1/tdip");

  href->SetLineColor(1);
  href->SetFillColor(623);
  href->SetFillStyle(3001);
  //  href->GetXaxis()->SetRangeUser(0,4.999);

  h->SetMarkerStyle(20);
  h->SetMarkerSize (1.);
  h->SetMarkerColor(kBlue+3);

  href->Draw();
  h->Draw("sames,ep");

  print_canvas_with_date(c1,oname);
}


//-----------------------------------------------------------------------------
void plot_validation_c02(const char* FnRef, const char* Fn) {

  const char* oname="plot_validation_c02";

  TCanvas* c1 = new TCanvas("c02","c02",1200,900);
  c1->Divide(2,2);
//-----------------------------------------------------------------------------
// DU
//-----------------------------------------------------------------------------
  c1->cd(1);

  TH1* href = gh1(FnRef,"TrackAna","trk_1/du");
  TH1* h    = gh1(Fn   ,"TrackAna","trk_1/du");

  href->SetLineColor(1);
  href->SetFillColor(623);
  href->SetFillStyle(3001);

  //  href->GetXaxis()->SetRangeUser(70,120);

  h->SetMarkerStyle(20);
  h->SetMarkerSize (1.);
  h->SetMarkerColor(kBlue+3);

  href->Draw();
  h->Draw("sames,ep");
//-----------------------------------------------------------------------------
// DV
//-----------------------------------------------------------------------------
  c1->cd(2);

  TH1* href = gh1(FnRef,"TrackAna","trk_1/dv");
  TH1* h    = gh1(Fn   ,"TrackAna","trk_1/dv");

  href->SetLineColor(1);
  href->SetFillColor(623);
  href->SetFillStyle(3001);
  //  href->GetXaxis()->SetRangeUser(0,119.9);

  h->SetMarkerStyle(20);
  h->SetMarkerSize (1.);
  h->SetMarkerColor(kBlue+3);

  href->Draw();
  h->Draw("sames,ep");
//-----------------------------------------------------------------------------
// DT
//-----------------------------------------------------------------------------
  c1->cd(3);

  TH1* href = gh1(FnRef,"TrackAna","trk_1/dt");
  TH1* h    = gh1(Fn   ,"TrackAna","trk_1/dt");

  href->SetLineColor(1);
  href->SetFillColor(623);
  href->SetFillStyle(3001);
  //  href->GetXaxis()->SetRangeUser(0,4.999);

  h->SetMarkerStyle(20);
  h->SetMarkerSize (1.);
  h->SetMarkerColor(kBlue+3);

  href->Draw();
  h->Draw("sames,ep");
//-----------------------------------------------------------------------------
// EP
//-----------------------------------------------------------------------------
  c1->cd(4);

  TH1* href = gh1(FnRef,"TrackAna","trk_1/ep");
  TH1* h    = gh1(Fn   ,"TrackAna","trk_1/ep");

  href->SetLineColor(1);
  href->SetFillColor(623);
  href->SetFillStyle(3001);
  //  href->GetXaxis()->SetRangeUser(0,4.999);

  h->SetMarkerStyle(20);
  h->SetMarkerSize (1.);
  h->SetMarkerColor(kBlue+3);

  href->Draw();
  h->Draw("sames,ep");

  print_canvas_with_date(c1,oname);
}

//-----------------------------------------------------------------------------
void plot_validation_c03(const char* FnRef, const char* Fn) {

  const char* oname="plot_validation_c03";

  TCanvas* c1 = new TCanvas("c03","c03",1200,900);
  c1->Divide(2,2);
//-----------------------------------------------------------------------------
// DU
//-----------------------------------------------------------------------------
  c1->cd(1);

  TH1* href = gh1(FnRef,"TrackAna","trk_1/d0");
  TH1* h    = gh1(Fn   ,"TrackAna","trk_1/d0");

  href->SetLineColor(1);
  href->SetFillColor(623);
  href->SetFillStyle(3001);

  //  href->GetXaxis()->SetRangeUser(70,120);

  h->SetMarkerStyle(20);
  h->SetMarkerSize (1.);
  h->SetMarkerColor(kBlue+3);

  href->Draw();
  h->Draw("sames,ep");
//-----------------------------------------------------------------------------
// DV
//-----------------------------------------------------------------------------
  c1->cd(2);

  TH1* href = gh1(FnRef,"TrackAna","trk_1/z0");
  TH1* h    = gh1(Fn   ,"TrackAna","trk_1/z0");

  href->SetLineColor(1);
  href->SetFillColor(623);
  href->SetFillStyle(3001);
  //  href->GetXaxis()->SetRangeUser(0,119.9);

  h->SetMarkerStyle(20);
  h->SetMarkerSize (1.);
  h->SetMarkerColor(kBlue+3);

  href->Draw();
  h->Draw("sames,ep");
//-----------------------------------------------------------------------------
// DT
//-----------------------------------------------------------------------------
  c1->cd(3);

  TH1* href = gh1(FnRef,"TrackAna","trk_1/t0");
  TH1* h    = gh1(Fn   ,"TrackAna","trk_1/t0");

  href->SetLineColor(1);
  href->SetFillColor(623);
  href->SetFillStyle(3001);
  //  href->GetXaxis()->SetRangeUser(0,4.999);

  h->SetMarkerStyle(20);
  h->SetMarkerSize (1.);
  h->SetMarkerColor(kBlue+3);

  href->Draw();
  h->Draw("sames,ep");
//-----------------------------------------------------------------------------
// EP
//-----------------------------------------------------------------------------
  c1->cd(4);

  TH1* href = gh1(FnRef,"TrackAna","trk_1/path");
  TH1* h    = gh1(Fn   ,"TrackAna","trk_1/path");

  href->SetLineColor(1);
  href->SetFillColor(623);
  href->SetFillStyle(3001);
  //  href->GetXaxis()->SetRangeUser(0,4.999);

  h->SetMarkerStyle(20);
  h->SetMarkerSize (1.);
  h->SetMarkerColor(kBlue+3);

  href->Draw();
  h->Draw("sames,ep");


  print_canvas_with_date(c1,oname);

}

