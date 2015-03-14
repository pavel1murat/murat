//
//-----------------------------------------------------------------------------
// trk_18: Set C + T0 > 700ns
//-----------------------------------------------------------------------------
void plot_ep_vs_lumi() {

  char hist_dir[200];

  sprintf(hist_dir,"%s/hist/mu2e",gSystem->Getenv("HOME"));

  TH1* h1 = gh1(Form("%s/v4_1_8/cnvs2402.track_ana.hist",hist_dir),"TrackAna","trk_18/ep");
  TH1* h2 = gh1(Form("%s/v4_2_1/cnvs1412.track_ana.hist",hist_dir),"TrackAna","trk_18/ep");
  TH1* h3 = gh1(Form("%s/v4_2_1/cnvs1512.track_ana.hist",hist_dir),"TrackAna","trk_18/ep");
  TH1* h4 = gh1(Form("%s/v4_2_1/cnvs1612.track_ana.hist",hist_dir),"TrackAna","trk_18/ep");

  h1->SetLineColor(2);
  h1->Draw();

  h2->Draw("same");

  //  h3->SetFillColor(kBlue-9);
  //  h3->SetFillStyle(3002);
  h3->Draw("same");

  h4->SetFillColor(kBlue-3);
  h4->SetFillStyle(3002);
  h4->Draw("same");

  TLegend* leg = new TLegend(0.15,0.6,0.4,0.8);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);

  leg->AddEntry(h1,"CE","f");
  leg->AddEntry(h2,"CE+MIXP3x1","f");
  leg->AddEntry(h3,"CE+MIXP3x2","f");
  leg->AddEntry(h4,"CE+MIXP3x4","f");
		
  leg->Draw();
}
