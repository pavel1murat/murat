//

namespace {
  int const  kNDatasets = 4;

  char* cut_label1[kNDatasets] = {
    "cnvs0202"   , "cnvs1402", "cnvs1502", "cnvs1602"
  };

  char* cut_label[kNDatasets] = {
    "CE"   , "CE+MIXP3", "CE+MIXP3-x2", "CE+MIXP3-x4"
  };

  const char* filename [kNDatasets] = {
    "/home/murat/hist/mu2e/v4_1_8/cnvs2402.track_ana.hist",
    "/home/murat/hist/mu2e/v4_1_8/cnvs1402.track_ana.hist",
    "/home/murat/hist/mu2e/v4_1_8/cnvs1502.track_ana.hist",
    "/home/murat/hist/mu2e/v4_1_8/cnvs1602.track_ana.hist"
  };

};

//-----------------------------------------------------------------------------
void plot_pid_eff_vs_dataset() {
  

  TH1F* h1;

  float v1, v2, qn0, qn1, err;

  TH1F* h_eff = new TH1F("heff","Efficiency vs dataset",kNDatasets,0,kNDatasets);

  for (int i=0; i<kNDatasets; i++) {

    h1 = gh1(filename[i],"TrackAna","trk_19/llhr");

    qn0 = h1->Integral();
    qn1 = h1->Integral(0,100);

    v1 = 1-qn1/qn0;
    v2 = sqrt(qn0-qn1)/qn0;

    h_eff->SetBinContent(i+1,v1);
    h_eff->SetBinError  (i+1,v2);
  }

  for (int i=1; i<=kNDatasets; i++) {
    h_eff->GetXaxis()->SetBinLabel(i,cut_label[i-1]);
  }

  h_eff->SetTitle(Form("PID Effficiency vs Dataset"));
  h_eff->SetStats(0);
  h_eff->SetMaximum(1.1);
  h_eff->SetMinimum(0.);
  h_eff->Draw();
  h_eff->Draw("same,text45");
}
