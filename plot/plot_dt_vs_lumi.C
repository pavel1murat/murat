//
//-----------------------------------------------------------------------------
// trk_18: Set C + T0 > 700ns
//-----------------------------------------------------------------------------
void plot_dt_vs_lumi() {

  char hist_dir[200];

  sprintf(hist_dir,"%s/hist/mu2e",gSystem->Getenv("HOME"));

  TH1* h1 = gh1(Form("%s/v4_1_8/cnvs2402.track_ana.hist",hist_dir),"TrackAna","trk_18/dt");
  TH1* h2 = gh1(Form("%s/v4_2_1/cnvs1412.track_ana.hist",hist_dir),"TrackAna","trk_18/dt");
  TH1* h3 = gh1(Form("%s/v4_2_1/cnvs1512.track_ana.hist",hist_dir),"TrackAna","trk_18/dt");
  TH1* h4 = gh1(Form("%s/v4_2_1/cnvs1612.track_ana.hist",hist_dir),"TrackAna","trk_18/dt");

  h1->SetLineColor(2);
  h1->GetXaxis()->SetTitle("#DeltaT, ns");
  h1->DrawNormalized("",1);

  h2->DrawNormalized("same",1);

  //  h3->SetFillColor(kBlue-9);
  //  h3->SetFillStyle(3002);
  h3->DrawNormalized("same",1);

  h4->SetFillColor(kBlue-3);
  h4->SetFillStyle(3002);
  h4->DrawNormalized("same",1);

  TLegend* leg = new TLegend(0.15,0.6,0.4,0.8);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);

  leg->AddEntry(h1,"CE","f");
  leg->AddEntry(h2,"CE+MIXP3x1","f");
  leg->AddEntry(h3,"CE+MIXP3x2","f");
  leg->AddEntry(h4,"CE+MIXP3x4","f");
		
  leg->Draw();
}

//-----------------------------------------------------------------------------
void plot_sigmat_vs_lumi() {

  char hist_dir[200];

  TH1F* h[4];

  sprintf(hist_dir,"%s/hist/mu2e",gSystem->Getenv("HOME"));

  h[0] = gh1(Form("%s/v4_1_8/cnvs2402.track_ana.hist",hist_dir),"TrackAna","trk_18/dt");
  h[1] = gh1(Form("%s/v4_2_1/cnvs1412.track_ana.hist",hist_dir),"TrackAna","trk_18/dt");
  h[2] = gh1(Form("%s/v4_2_1/cnvs1512.track_ana.hist",hist_dir),"TrackAna","trk_18/dt");
  h[3] = gh1(Form("%s/v4_2_1/cnvs1612.track_ana.hist",hist_dir),"TrackAna","trk_18/dt");

  double x[4], ex[4], s[4], e[4];

  x[0] = 0;
  x[1] = 1;
  x[2] = 2;
  x[3] = 4;

  for (int i=0; i<4; i++) {
    h[i]->Fit("gaus","Q");
    s[i] = h[i]->GetFunction("gaus")->GetParameter(2);
    e[i] = h[i]->GetFunction("gaus")->GetParError (2);
    ex[i] = 0;
  }

  TH2F* h2 = new TH2F("st_vs_occ","timing resolution vs occupancy",1,-1,5,1,0,1);

  h2->GetXaxis()->SetTitle("background scale factor, x0 = CE-only");
  h2->GetYaxis()->SetTitle("#sigma_{T}, ns");
  h2->Draw();

  TGraphErrors* gr = new TGraphErrors(4,x,s,ex,e);

  gr->Draw("LP,same");

}
