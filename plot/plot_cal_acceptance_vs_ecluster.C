//
void plot_cal_acceptance_vs_ecluster(const char* Fn) {
  TH1F* h = gh1(Fn,"TrackAna","trk_1/ep");

  double eff, err;

  double qtot = h->Integral();

  TH1F* heff = (TH1F*) h->Clone("heff");

  heff->Reset();

  // I forgot to include underflows to EP hist, so normalization comes from here

  double q0 = gh1(Fn,"TrackAna","trk_1/p")->GetEntries();

  int nb = h->GetNbinsX();

  for (int i=1; i<nb; i++) {
    double qn = h->Integral(1,i);

    eff = (qtot-qn)/q0;
    err = sqrt(eff*(1-eff)/q0);

    heff->SetBinContent(i,eff);
    heff->SetBinError  (i,err);
  }

  heff->SetTitle("Calorimeter Efficiency vs E/P");
  heff->GetXaxis()->SetTitle("E_{cluster}/P_{track}");
  heff->GetXaxis()->SetRangeUser(0.,0.999);
  heff->Draw();


}

