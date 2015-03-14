//

//-----------------------------------------------------------------------------
plot_p_setc(const char* Fn) {

  TCanvas* c = new TCanvas("c","c",900,600);
  
  TH1F* h_tpr  = gh1(Fn,"TrackAna","trk_62/p2");
  TH1F* h_cpr  = gh1(Fn,"TrackAna","trk_52/p2");
  
  //  h_ef->GetXaxis()->SetRangeUser(-15,14.99);
  h_tpr->GetXaxis()->SetTitle("P [MeV/c]");
  h_tpr->GetXaxis()->SetRangeUser(90.,110.);

  double binsize=h_tpr->GetBinWidth(1);

  h_tpr->GetYaxis()->SetTitle(Form("Entries /( %3.2f MeV/c)", binsize));
  h_tpr->SetName("TrkPatRec");
  h_tpr->SetTitle("Momentum distribution for cut Set C tracks");
  
  h_cpr->SetName("CalPatRec");
  h_cpr->SetFillStyle(3002);
  h_cpr->SetFillColor(kRed);
  h_cpr->SetLineColor(kRed);

  h_tpr->SetStats(0);
  h_cpr->SetStats(0);

  h_tpr->Draw();
  h_cpr->Draw("sames");


  TLegend *leg = new TLegend(0.2, 0.72, 0.39, 0.85, NULL, "brNDC");
  leg->SetFillStyle(1001);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(1);
  leg->SetShadowColor(0);

  leg->AddEntry(h_tpr,"TrkPatRec","f");
  leg->AddEntry(h_cpr,"CalPatRec","f");
  leg->Draw();
}


//-----------------------------------------------------------------------------
plot_chi2_all(const char* Fn) {

  TCanvas* c = new TCanvas("c","c",900,600);
  
  TH1F* h_tpr  = gh1(Fn,"TrackAna","trk_60/chi2d");
  TH1F* h_cpr  = gh1(Fn,"TrackAna","trk_50/chi2d");
  
  //  h_ef->GetXaxis()->SetRangeUser(-15,14.99);
  //  h_tpr->GetXaxis()->SetTitle("P [MeV/c]");
  //  h_tpr->GetXaxis()->SetRangeUser(90.,110.);

  double binsize=h_tpr->GetBinWidth(1);

  //  h_tpr->GetYaxis()->SetTitle(Form("Entries /( %3.2f MeV/c)", binsize));
  h_tpr->SetName("TrkPatRec");
  h_tpr->SetTitle("#chi^{2}/N_{DOF} for reconstructed tracks");
  
  h_cpr->SetName("CalPatRec");
  h_cpr->SetFillStyle(3002);
  h_cpr->SetFillColor(kRed);
  h_cpr->SetLineColor(kRed);

  h_tpr->SetStats(0);
  h_cpr->SetStats(0);

  h_tpr->Draw();
  h_cpr->Draw("sames");


  TLegend *leg = new TLegend(0.6, 0.72, 0.8, 0.85, NULL, "brNDC");
  leg->SetFillStyle(1001);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(1);
  leg->SetShadowColor(0);

  leg->AddEntry(h_tpr,"TrkPatRec","f");
  leg->AddEntry(h_cpr,"CalPatRec","f");
  leg->Draw();
}


//-----------------------------------------------------------------------------
plot_nactive_setc(const char* Fn) {

  TCanvas* c = new TCanvas("c","c",900,600);
  
  TH1F* h_tpr  = gh1(Fn,"TrackAna","trk_62/nactv");
  TH1F* h_cpr  = gh1(Fn,"TrackAna","trk_52/nactv");
  
  //  h_ef->GetXaxis()->SetRangeUser(-15,14.99);
  h_tpr->GetXaxis()->SetTitle("N(hits)");
  //  h_tpr->GetXaxis()->SetRangeUser(90.,110.);

  double binsize=h_tpr->GetBinWidth(1);

  h_tpr->GetYaxis()->SetTitle(Form("Entries /( %3.2f MeV/c)", binsize));
  h_tpr->SetName("TrkPatRec");
  h_tpr->SetTitle("Number of reconstructed hits on Set C track");
  
  h_cpr->SetName("CalPatRec");
  h_cpr->SetFillStyle(3002);
  h_cpr->SetFillColor(kRed);
  h_cpr->SetLineColor(kRed);

  h_tpr->SetStats(0);
  h_cpr->SetStats(0);

  h_tpr->Draw();
  h_cpr->Draw("sames");


  TLegend *leg = new TLegend(0.7, 0.72, 0.9, 0.85, NULL, "brNDC");
  leg->SetFillStyle(1001);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(1);
  leg->SetShadowColor(0);

  leg->AddEntry(h_tpr,"TrkPatRec","f");
  leg->AddEntry(h_cpr,"CalPatRec","f");
  leg->Draw();
}


//-----------------------------------------------------------------------------
// here we use several files
//-----------------------------------------------------------------------------
void plot_calpatrec_eff() {

  int const npt = 4;

  float x    [npt] = {0.    , 1.  , 2.  , 4   };
  float ratio[npt] = {1.0165, 1.11, 1.28, 1.90};
  //  float eff  [npt] = {0.1099, 0.1052, 0.1004, 0.0858};
  float eff  [npt] = {0.0987, 0.0945, 0.0903, 0.0791};

  float  ex  [npt] = { 0., 0., 0., 0};
  float  ey  [npt];

  float y[4];

  float eff0 = 0.1095;

  float sf   = 0.8829; // 0.874;

  for (int i=0; i<4; i++) {
    y[i] = eff[i];
    ey [i] = 0.001;
  }
  

  TGraphErrors* gr = new TGraphErrors(npt,x,y,ex,ey);

  gr->SetMarkerStyle(20);
  gr->SetMarkerSize (1);

  gr->SetTitle("");
  gr->GetXaxis()->SetTitle("Background / (nominal expected background)");


  gr->GetYaxis()->SetRangeUser(0.05,0.12);
  gr->Draw("alp");

  TText* txt = new TText(0.03,0.115,"combined TrkPatRec+CalPatRec reconstruction efficiency for CE");
  txt->SetTextSize(0.04);
  txt->SetTextFont(52);

  txt->Draw();
}
