//


#include "murat/scripts/plot_utilities.hh"
#include "murat/scripts/fit_cb2.C"

void plot_dpf_fit(int Print = 1) {

  TCanvas* c = new TCanvas("c","",1400,800);

  cb2_fit_crystal_ball("~/hist/mu2e/v6_0_0/e00s6000.track_comp_use_mva.hist",
		       "TrackComp",
		       "trk_3/dpf",0,-5,5);

  TH1F* h = (TH1F*) gROOT->FindObject("h_cb2_fit_crystal_ball");

  h->SetTitle("");
  h->GetXaxis()->SetTitle("#Delta P, MeV/c");

  draw_label_ndc("#Delta P = P_{reco} - P_{MC}",0.15,0.8,0.06);

  c->SetLogy(1);
  c->Modified();
  c->Update();

  if (Print != 0) {
    c->Print("resolution_at_the_tracker_front_cb2_fit.eps");
  }
}
