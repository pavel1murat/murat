//

tdr_plot_delta_t() {

  TCanvas* c = new TCanvas("c","c",900,600);

  TH1F* h_e = gh1("hist/e00s0002.track_ana.hist","TrackAna","trk_1/dt");
  TH1F* h_m = gh1("hist/m00s0002.track_ana.hist","TrackAna","trk_1/dt");

  h_e->GetXaxis()->SetRangeUser(-15,14.99);
  h_e->GetXaxis()->SetTitle("#Delta T, ns");
  h_e->SetTitle("T(track)-T(cluster)");

  h_e->DrawNormalized("",1);
  h_m->SetFillStyle(3002);
  h_m->SetFillColor(kBlue-7);
  h_m->DrawNormalized("same",1);

  c->Draw();
}

//-----------------------------------------------------------------------------
// printing: 
// c->Print("/home/murat/figures/mu2e/tdr/trackana_trk_1_ep.eps")
//-----------------------------------------------------------------------------
tdr_plot_ep() {

  TCanvas* c = new TCanvas("c","c",900,600);

  TH1F* h_e = gh1("hist/e00s0002.track_ana.hist","TrackAna","trk_1/ep");
  TH1F* h_m = gh1("hist/m00s0002.track_ana.hist","TrackAna","trk_1/ep");

  h_e->GetXaxis()->SetRangeUser(0,1.199);
  h_e->GetXaxis()->SetTitle("E/P");
  h_e->SetTitle("E/P");

  h_e->SetMaximum(10000);
  h_e->DrawNormalized("",1);
  h_m->SetFillStyle(3002);
  h_m->SetFillColor(kBlue-7);
  h_m->DrawNormalized("same",1);

  c->Draw();
}
