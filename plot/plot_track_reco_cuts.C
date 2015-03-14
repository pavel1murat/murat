//-----------------------------------------------------------------------------
// 1. N(active) > 20 cut
//-----------------------------------------------------------------------------
void plot_nactive_cut(const char* Filename, int PrintFlag) {

  TH1F *h0, *h1, *h4;

  h0 = gh1(Filename,"TrackRecoEffAna","eff_0/nactive_0");
  h1 = gh1(Filename,"TrackRecoEffAna","eff_0/nactive_1");
  h4 = gh1(Filename,"TrackRecoEffAna","eff_0/nactive_4");


  TCanvas* c = new TCanvas("c_eff_0_nactive","c_eff_0_nactive",0,0,1400,800);

  c->cd();

  //  gPad->SetLogy(1);

  h0->SetTitle("N(active) > 20");
  //  h0->GetXaxis()->SetRangeUser(50,109.99);
  h0->Draw();

  h1->SetFillStyle(3001);
  h1->SetFillColor(2);
  h1->SetLineColor(2);
  h1->Draw("sames");


  h4->SetFillStyle(1001);
  h4->SetFillColor(kBlue-6);
  h4->SetLineColor(kBlue-6);

  h4->Draw("sames");

  TArrow  *a1;
  a1 = new TArrow(20,50,20,20,0.015);
  a1->SetLineWidth(3);
  a1->Draw();

  c->Modified();
  c->Update();

  TPaveStats* s0 = (TPaveStats*) h0->GetListOfFunctions()->FindObject("stats");

  s0->SetX1NDC(0.70);
  s0->SetX2NDC(0.90);
  s0->SetY1NDC(0.75);
  s0->SetY2NDC(0.90);
  s0->Draw();

  TPaveStats* s1 = (TPaveStats*) h1->GetListOfFunctions()->FindObject("stats");

  s1->SetX1NDC(0.70);
  s1->SetX2NDC(0.90);
  s1->SetY1NDC(0.60);
  s1->SetY2NDC(0.75);
  s1->Draw();

  TPaveStats* s4 = (TPaveStats*) h4->GetListOfFunctions()->FindObject("stats");

  s4->SetX1NDC(0.70);
  s4->SetX2NDC(0.90);
  s4->SetY1NDC(0.45);
  s4->SetY2NDC(0.60);
  s4->Draw();

  if (PrintFlag) {
    c->Print("track_reco_cuts_01_nactive.eps");
  }
}

//-----------------------------------------------------------------------------
// 2. FitConsistency cut
//-----------------------------------------------------------------------------
void plot_fcons_cut(const char* Filename, int PrintFlag) {

  TH1F *h0, *h1, *h4;

  h0 = gh1(Filename,"TrackRecoEffAna","eff_0/fcons_0");
  h1 = gh1(Filename,"TrackRecoEffAna","eff_0/fcons_1");
  h4 = gh1(Filename,"TrackRecoEffAna","eff_0/fcons_4");


  TCanvas* c = new TCanvas("c_eff_0_fitcons","c_eff_0_fitcons",0,0,1400,800);

  c->cd();

  gPad->SetLogx(1);

  h0->SetTitle("fit probability > 2.e-3");
  //  h0->GetXaxis()->SetRangeUser(50,109.99);
  h0->Draw();

  h1->SetFillStyle(3001);
  h1->SetFillColor(2);
  h1->SetLineColor(2);
  h1->Draw("sames");


  h4->SetFillStyle(1001);
  h4->SetFillColor(kBlue-6);
  h4->SetLineColor(kBlue-6);

  h4->Draw("sames");

  TArrow  *a1;
  a1 = new TArrow(0.002,50,0.002,10,0.015);
  a1->SetLineWidth(3);
  a1->Draw();

  c->Modified();
  c->Update();

  TPaveStats* s0 = (TPaveStats*) h0->GetListOfFunctions()->FindObject("stats");

  s0->SetX1NDC(0.70);
  s0->SetX2NDC(0.90);
  s0->SetY1NDC(0.75);
  s0->SetY2NDC(0.90);
  s0->Draw();

  TPaveStats* s1 = (TPaveStats*) h1->GetListOfFunctions()->FindObject("stats");

  s1->SetX1NDC(0.70);
  s1->SetX2NDC(0.90);
  s1->SetY1NDC(0.60);
  s1->SetY2NDC(0.75);
  s1->Draw();

  TPaveStats* s4 = (TPaveStats*) h4->GetListOfFunctions()->FindObject("stats");

  s4->SetX1NDC(0.70);
  s4->SetX2NDC(0.90);
  s4->SetY1NDC(0.45);
  s4->SetY2NDC(0.60);
  s4->Draw();

  if (PrintFlag) {
    c->Print("track_reco_cuts_02_fcons.eps");
  }
}


//-----------------------------------------------------------------------------
// 3. FitMomErr cut
//-----------------------------------------------------------------------------
void plot_fitmomerr_cut(const char* Filename, int PrintFlag) {

  TH1F *h0, *h1, *h4;

  h0 = gh1(Filename,"TrackRecoEffAna","eff_0/fitmomerr_0");
  h1 = gh1(Filename,"TrackRecoEffAna","eff_0/fitmomerr_1");
  h4 = gh1(Filename,"TrackRecoEffAna","eff_0/fitmomerr_4");


  TCanvas* c = new TCanvas("c_eff_0_fitmomerr","c_eff_0_fitmomerr",0,0,1400,800);

  c->cd();

  gPad->SetLogy(1);

  h0->SetTitle("#sigma(P) < 0.25 MeV/c");
  //  h0->GetXaxis()->SetRangeUser(50,109.99);
  h0->Draw();

  h1->SetFillStyle(3001);
  h1->SetFillColor(2);
  h1->SetLineColor(2);
  h1->Draw("sames");


  h4->SetFillStyle(1001);
  h4->SetFillColor(kBlue-6);
  h4->SetLineColor(kBlue-6);

  h4->Draw("sames");

  TArrow  *a1;
  a1 = new TArrow(0.25,50,0.25,10,0.015);
  a1->SetLineWidth(3);
  a1->Draw();

  c->Modified();
  c->Update();

  TPaveStats* s0 = (TPaveStats*) h0->GetListOfFunctions()->FindObject("stats");

  s0->SetX1NDC(0.70);
  s0->SetX2NDC(0.90);
  s0->SetY1NDC(0.70);
  s0->SetY2NDC(0.90);
  s0->Draw();

  TPaveStats* s1 = (TPaveStats*) h1->GetListOfFunctions()->FindObject("stats");

  s1->SetX1NDC(0.70);
  s1->SetX2NDC(0.90);
  s1->SetY1NDC(0.50);
  s1->SetY2NDC(0.70);
  s1->Draw();

  TPaveStats* s4 = (TPaveStats*) h4->GetListOfFunctions()->FindObject("stats");

  s4->SetX1NDC(0.70);
  s4->SetX2NDC(0.90);
  s4->SetY1NDC(0.30);
  s4->SetY2NDC(0.50);
  s4->Draw();

  if (PrintFlag) {
    c->Print("track_reco_cuts_03_fitmomerr.eps");
  }
}

//-----------------------------------------------------------------------------
// 4. T0Err cut
//-----------------------------------------------------------------------------
void plot_t0err_cut(const char* Filename, int PrintFlag) {

  TH1F *h0, *h1, *h4;

  h0 = gh1(Filename,"TrackRecoEffAna","eff_0/t0err_0");
  h1 = gh1(Filename,"TrackRecoEffAna","eff_0/t0err_1");
  h4 = gh1(Filename,"TrackRecoEffAna","eff_0/t0err_4");


  TCanvas* c = new TCanvas("c_eff_0_t0err","c_eff_0_t0err",0,0,1400,800);

  c->cd();

  gPad->SetLogy(1);

  h0->SetTitle("#sigma(T0) < 0.9 ns");
  //  h0->GetXaxis()->SetRangeUser(50,109.99);
  h0->Draw();

  h1->SetFillStyle(3001);
  h1->SetFillColor(2);
  h1->SetLineColor(2);
  h1->Draw("sames");


  h4->SetFillStyle(1001);
  h4->SetFillColor(kBlue-6);
  h4->SetLineColor(kBlue-6);

  h4->Draw("sames");

  TArrow *a1;
  a1 = new TArrow(0.9.,50,0.9,10,0.015);
  a1->SetLineWidth(3);
  a1->Draw();

  c->Modified();
  c->Update();

  TPaveStats* s0 = (TPaveStats*) h0->GetListOfFunctions()->FindObject("stats");

  s0->SetX1NDC(0.70);
  s0->SetX2NDC(0.90);
  s0->SetY1NDC(0.70);
  s0->SetY2NDC(0.90);
  s0->Draw();

  TPaveStats* s1 = (TPaveStats*) h1->GetListOfFunctions()->FindObject("stats");

  s1->SetX1NDC(0.70);
  s1->SetX2NDC(0.90);
  s1->SetY1NDC(0.50);
  s1->SetY2NDC(0.70);
  s1->Draw();

  TPaveStats* s4 = (TPaveStats*) h4->GetListOfFunctions()->FindObject("stats");

  s4->SetX1NDC(0.70);
  s4->SetX2NDC(0.90);
  s4->SetY1NDC(0.30);
  s4->SetY2NDC(0.50);
  s4->Draw();

  if (PrintFlag) {
    c->Print("track_reco_cuts_04_t0err.eps");
  }
}

//-----------------------------------------------------------------------------
// 5. T0 cut
//-----------------------------------------------------------------------------
void plot_t0_cut(const char* Filename, int PrintFlag) {

  TH1F *h0, *h1, *h4;

  h0 = gh1(Filename,"TrackRecoEffAna","eff_0/t0_0");
  h1 = gh1(Filename,"TrackRecoEffAna","eff_0/t0_1");
  h4 = gh1(Filename,"TrackRecoEffAna","eff_0/t0_4");


  TCanvas* c = new TCanvas("c_eff_0_t0","c_eff_0_t0",0,0,1400,800);

  c->cd();

  // gPad->SetLogy(1);

  h0->SetTitle("T0 > 700 ns");
  //  h0->GetXaxis()->SetRangeUser(50,109.99);
  h0->Draw();

  h1->SetFillStyle(3001);
  h1->SetFillColor(2);
  h1->SetLineColor(2);
  h1->Draw("sames");


  h4->SetFillStyle(1001);
  h4->SetFillColor(kBlue-6);
  h4->SetLineColor(kBlue-6);

  h4->Draw("sames");

  TArrow *a1;
  a1 = new TArrow(700.,40,700,20,0.015);
  a1->SetLineWidth(3);
  a1->Draw();

  c->Modified();
  c->Update();

  TPaveStats* s0 = (TPaveStats*) h0->GetListOfFunctions()->FindObject("stats");

  s0->SetX1NDC(0.70);
  s0->SetX2NDC(0.90);
  s0->SetY1NDC(0.70);
  s0->SetY2NDC(0.90);
  s0->Draw();

  TPaveStats* s1 = (TPaveStats*) h1->GetListOfFunctions()->FindObject("stats");

  s1->SetX1NDC(0.70);
  s1->SetX2NDC(0.90);
  s1->SetY1NDC(0.50);
  s1->SetY2NDC(0.70);
  s1->Draw();

  TPaveStats* s4 = (TPaveStats*) h4->GetListOfFunctions()->FindObject("stats");

  s4->SetX1NDC(0.70);
  s4->SetX2NDC(0.90);
  s4->SetY1NDC(0.30);
  s4->SetY2NDC(0.50);
  s4->Draw();

  if (PrintFlag) {
    c->Print("track_reco_cuts_05_t0.eps");
  }
}

//-----------------------------------------------------------------------------
// 6. TanDIP cut
//-----------------------------------------------------------------------------
void plot_tandip_cut(const char* Filename, int PrintFlag) {

  TH1F *h0, *h1, *h4;

  h0 = gh1(Filename,"TrackRecoEffAna","eff_0/tandip_0");
  h1 = gh1(Filename,"TrackRecoEffAna","eff_0/tandip_1");
  h4 = gh1(Filename,"TrackRecoEffAna","eff_0/tandip_4");


  TCanvas* c = new TCanvas("c_eff_0_tandip","c_eff_0_tandip",0,0,1400,800);

  c->cd();

  //  gPad->SetLogy(1);

  h0->SetTitle("1 < tan(#lambda) < sqrt(3)");
  //  h0->GetXaxis()->SetRangeUser(50,109.99);
  h0->Draw();

  h1->SetFillStyle(3001);
  h1->SetFillColor(2);
  h1->SetLineColor(2);
  h1->Draw("sames");


  h4->SetFillStyle(1001);
  h4->SetFillColor(kBlue-6);
  h4->SetLineColor(kBlue-6);

  h4->Draw("sames");

  TArrow  *a1, *a2;
  a1 = new TArrow(1.,40,1,20,0.015);
  a1->SetLineWidth(3);
  a1->Draw();

  a2 = new TArrow(sqrt(3),40,sqrt(3.),20,0.015);
  a2->SetLineWidth(3);
  a2->Draw();

3  c->Modified();
  c->Update();

  TPaveStats* s0 = (TPaveStats*) h0->GetListOfFunctions()->FindObject("stats");

  s0->SetX1NDC(0.10);
  s0->SetX2NDC(0.30);
  s0->SetY1NDC(0.70);
  s0->SetY2NDC(0.90);
  s0->Draw();

  TPaveStats* s1 = (TPaveStats*) h1->GetListOfFunctions()->FindObject("stats");

  s1->SetX1NDC(0.10);
  s1->SetX2NDC(0.30);
  s1->SetY1NDC(0.50);
  s1->SetY2NDC(0.70);
  s1->Draw();

  TPaveStats* s4 = (TPaveStats*) h4->GetListOfFunctions()->FindObject("stats");

  s4->SetX1NDC(0.10);
  s4->SetX2NDC(0.30);
  s4->SetY1NDC(0.30);
  s4->SetY2NDC(0.50);
  s4->Draw();

  if (PrintFlag) {
    c->Print("track_reco_cuts_06_tandip.eps");
  }
}

//-----------------------------------------------------------------------------
// 7. D0 cut
//-----------------------------------------------------------------------------
void plot_d0_cut(const char* Filename, int PrintFlag) {

  TH1F *h0, *h1, *h4;

  h0 = gh1(Filename,"TrackRecoEffAna","eff_0/d0_0");
  h1 = gh1(Filename,"TrackRecoEffAna","eff_0/d0_1");
  h4 = gh1(Filename,"TrackRecoEffAna","eff_0/d0_4");


  TCanvas* c = new TCanvas("c_eff_0_d0","c_eff_0_d0",0,0,1400,800);

  c->cd();

  //  gPad->SetLogy(1);

  h0->SetTitle("-80 < R(max) < 105 mm");
  //  h0->GetXaxis()->SetRangeUser(50,109.99);
  h0->Draw();

  h1->SetFillStyle(3001);
  h1->SetFillColor(2);
  h1->SetLineColor(2);
  h1->Draw("sames");


  h4->SetFillStyle(1001);
  h4->SetFillColor(kBlue-6);
  h4->SetLineColor(kBlue-6);

  h4->Draw("sames");

  c->Modified();
  c->Update();

  TPaveStats* s0 = (TPaveStats*) h0->GetListOfFunctions()->FindObject("stats");

  s0->SetX1NDC(0.10);
  s0->SetX2NDC(0.30);
  s0->SetY1NDC(0.70);
  s0->SetY2NDC(0.90);
  s0->Draw();

  TPaveStats* s1 = (TPaveStats*) h1->GetListOfFunctions()->FindObject("stats");

  s1->SetX1NDC(0.10);
  s1->SetX2NDC(0.30);
  s1->SetY1NDC(0.50);
  s1->SetY2NDC(0.70);
  s1->Draw();

  TPaveStats* s4 = (TPaveStats*) h4->GetListOfFunctions()->FindObject("stats");

  s4->SetX1NDC(0.10);
  s4->SetX2NDC(0.30);
  s4->SetY1NDC(0.30);
  s4->SetY2NDC(0.50);
  s4->Draw();

  if (PrintFlag) {
    c->Print("track_reco_cuts_07_d0.eps");
  }
}


//-----------------------------------------------------------------------------
// 8. RMax cut
//-----------------------------------------------------------------------------
void plot_rmax_cut(const char* Filename, int PrintFlag) {

  TH1F *h0, *h1, *h4;

  h0 = gh1(Filename,"TrackRecoEffAna","eff_0/rmax_0");
  h1 = gh1(Filename,"TrackRecoEffAna","eff_0/rmax_1");
  h4 = gh1(Filename,"TrackRecoEffAna","eff_0/rmax_4");


  TCanvas* c = new TCanvas("c_eff_0_rmax","c_eff_0_rmax",0,0,1400,800);

  c->cd();

  //  gPad->SetLogy(1);

  h0->SetTitle("450 < R(max) < 680 mm");
  //  h0->GetXaxis()->SetRangeUser(50,109.99);
  h0->Draw();

  h1->SetFillStyle(3001);
  h1->SetFillColor(2);
  h1->SetLineColor(2);
  h1->Draw("sames");


  h4->SetFillStyle(1001);
  h4->SetFillColor(kBlue-6);
  h4->SetLineColor(kBlue-6);

  h4->Draw("sames");

  c->Modified();
  c->Update();

  TPaveStats* s0 = (TPaveStats*) h0->GetListOfFunctions()->FindObject("stats");

  s0->SetX1NDC(0.10);
  s0->SetX2NDC(0.30);
  s0->SetY1NDC(0.70);
  s0->SetY2NDC(0.90);
  s0->Draw();

  TPaveStats* s1 = (TPaveStats*) h1->GetListOfFunctions()->FindObject("stats");

  s1->SetX1NDC(0.10);
  s1->SetX2NDC(0.30);
  s1->SetY1NDC(0.50);
  s1->SetY2NDC(0.70);
  s1->Draw();

  TPaveStats* s4 = (TPaveStats*) h4->GetListOfFunctions()->FindObject("stats");

  s4->SetX1NDC(0.10);
  s4->SetX2NDC(0.30);
  s4->SetY1NDC(0.30);
  s4->SetY2NDC(0.50);
  s4->Draw();

  if (PrintFlag) {
    c->Print("track_reco_cuts_08_rmax.eps");
  }
}


//-----------------------------------------------------------------------------
// 9. Signal momentum window
//-----------------------------------------------------------------------------
void plot_track_mom_cut(const char* Filename, int PrintFlag = 0) {

  TH1F *h0, *h1, *h4;

  h0 = gh1(Filename,"TrackRecoEffAna","eff_0/p_0");
  h1 = gh1(Filename,"TrackRecoEffAna","eff_0/p_1");
  h4 = gh1(Filename,"TrackRecoEffAna","eff_0/p_4");

  TCanvas* c = new TCanvas("c_eff_0_p","c_eff_0_p",0,0,1400,800);

  c->cd();

  gPad->SetLogy(1);

  h0->SetTitle("103.5 < P < 105.0 MeV/c");
  h0->GetXaxis()->SetRangeUser(50,109.99);
  h0->Draw();

  h1->SetFillStyle(3001);
  h1->SetFillColor(2);
  h1->SetLineColor(2);
  h1->Draw("sames");


  h4->SetFillStyle(1001);
  h4->SetFillColor(kBlue-6);
  h4->SetLineColor(kBlue-6);

  h4->Draw("sames");

  c->Modified();
  c->Update();

  TPaveStats* s0 = (TPaveStats*) h0->GetListOfFunctions()->FindObject("stats");

  s0->SetX1NDC(0.10);
  s0->SetX2NDC(0.30);
  s0->SetY1NDC(0.70);
  s0->SetY2NDC(0.90);
  s0->Draw();

  TPaveStats* s1 = (TPaveStats*) h1->GetListOfFunctions()->FindObject("stats");

  s1->SetX1NDC(0.10);
  s1->SetX2NDC(0.30);
  s1->SetY1NDC(0.50);
  s1->SetY2NDC(0.70);
  s1->Draw();

  TPaveStats* s4 = (TPaveStats*) h4->GetListOfFunctions()->FindObject("stats");

  s4->SetX1NDC(0.10);
  s4->SetX2NDC(0.30);
  s4->SetY1NDC(0.30);
  s4->SetY2NDC(0.50);
  s4->Draw();

  if (PrintFlag) {
    c->Print("track_reco_cuts_09_reco_mom.eps");
  }
}


//-----------------------------------------------------------------------------
void plot_all(const char* Fn) {

  plot_nactive_cut(Fn,1);
  plot_fcons_cut(Fn,1);
  plot_fitmomerr_cut(Fn,1);
  plot_t0err_cut(Fn,1);
  plot_t0_cut(Fn,1);
  plot_tandip_cut(Fn,1);
  plot_d0_cut(Fn,1);
  plot_rmax_cut(Fn,1);
  plot_track_mom_cut(Fn,1);
}

