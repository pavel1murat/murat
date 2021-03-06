//

namespace {
  int const  kNDatasets = 4;

  char* cut_label1[kNDatasets] = {
    "cnvs0202"   , "cnvs1402", "cnvs1502", "cnvs1602"
  };

  char* cut_label[kNDatasets] = {
    "CE"   , "CE+MIXP3", "CE+MIXP3-x2", "CE+MIXP3-x4"
  };

  char* cut_label_cd2[kNDatasets] = {
    "CE"   , "CE+MIXP3", "CE+MIXP3-x2", "CE+MIXP3-x3"
  };

  const int qevents[kNDatasets] = {
    200000*1.4, 20000, 20000, 20000
  };

  const char* filename [kNDatasets] = {
    "/home/murat/hist/mu2e/v4_1_8/cnvs0202.track_ana.hist",
    "/home/murat/hist/mu2e/v4_2_1/cnvs1412.track_ana.hist",
    "/home/murat/hist/mu2e/v4_2_1/cnvs1512.track_ana.hist",
    "/home/murat/hist/mu2e/v4_2_1/cnvs1612.track_ana.hist"
  };

};

//-----------------------------------------------------------------------------
void plot_trk_eff_vs_dataset() {
  

  TH1F* h_mpr = new TH1F("h_mpr","MergePatRec Efficiency",4,0,4);
  TH1F* h_tpr = new TH1F("h_tpr","TrkPatRec   Efficiency",4,0,4);
  TH1F* hh    = h_tpr->Clone("hh");

  float dat_mpr[10], dat_tpr[10], mpr_over_tpr[10];
  float v1, v2, err;

  for (int i=0; i<kNDatasets; i++) {

    dat_mpr[0] = qevents[i];
    dat_tpr[0] = dat_mpr[0];

    dat_mpr[9] = gh1(filename[i],"TrackAna","evt_18/ce_costh")->GetEntries();
    dat_tpr[9] = gh1(filename[i],"TrackAna","evt_28/ce_costh")->GetEntries();

    v1 = dat_mpr[9]/(dat_mpr[0]+0.);
    v2 = dat_tpr[9]/(dat_tpr[0]+0.);

    mpr_over_tpr[i] = v1/v2;

    err = v1/sqrt(qevents[i]*v1);
      
    h_mpr->SetBinContent(i+1,v1);
    //    h_mpr->SetBinError  (i+1,err);
    h_tpr->SetBinContent(i+1,v2);
    hh->SetBinContent(i+1,v2);
    hh->SetBinError  (i+1,err);
  }

  for (int i=1; i<=kNDatasets; i++) {
    h_mpr->GetXaxis()->SetBinLabel(i,cut_label[i-1]);
  }

  h_mpr->SetTitle(Form("MergePatRec+CalPatRec efficiency"));
  h_mpr->SetStats(0);
  h_mpr->SetMaximum(0.16);
  h_mpr->SetMinimum(0.);
  h_mpr->Draw();
  h_mpr->Draw("same,text45");

  h_tpr->SetFillStyle(3001);
  h_tpr->SetFillColor(kBlue-7);
  h_tpr->Draw("same");

  hh->SetMarkerStyle(20);
  hh->SetMarkerSize(1.);
  hh->Draw("same,ep");


  TLegend* leg = new TLegend(0.7,0.75,0.9,0.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);

  leg->AddEntry(h_mpr,"TrkPatRec+CalPatRec","f");
  leg->AddEntry(h_tpr,"TrkPatRec","f");

  leg->Draw();
}


//-----------------------------------------------------------------------------
void plot_trk_eff_ratio_vs_dataset() {
  

  TH1F* h_rat = new TH1F("h_mpr","MergePatRec/TrkPatRec Ratio",4,0,4);
  TH1F* hh    = h_rat->Clone("hh");

  float dat_mpr[10], dat_tpr[10], mpr_over_tpr[10];
  float v1, v2, err;

  for (int i=0; i<kNDatasets; i++) {

    dat_mpr[0] = qevents[i];
    dat_tpr[0] = dat_mpr[0];

    dat_mpr[9] = gh1(filename[i],"TrackAna","evt_18/ce_costh")->GetEntries();
    dat_tpr[9] = gh1(filename[i],"TrackAna","evt_28/ce_costh")->GetEntries();

    v1 = dat_mpr[9]/(dat_mpr[0]+0.);
    v2 = dat_tpr[9]/(dat_tpr[0]+0.);

    mpr_over_tpr[i] = v1/v2;

    printf(" %10.5f\n", mpr_over_tpr[i]);

    err = v1/sqrt(qevents[i]*v1);
      
    hh->SetBinContent(i+1,mpr_over_tpr[i]);

  }

  for (int i=1; i<=kNDatasets; i++) {
    hh->GetXaxis()->SetBinLabel(i,cut_label[i-1]);
  }

  hh->SetTitle(Form("#epsilon(TrkPatRec+CalPatRec)/#epsilon(TrkPatRec)"));
  hh->SetMinimum(0.5);
  hh->SetMaximum(2.5);
  hh->SetStats(0);
  hh->Draw();
  hh->Draw("same,text45");
}


//-----------------------------------------------------------------------------
void plot_trk_eff_vs_dataset_cd2() {
  
  float dat_mpr[4] = { 0.1177, 0.1116, 0.1047, 0.0978};
  float dat_tpr[4] = { 0.1140, 0.104 , 0.092 , 0.075 };

  TH1F* h_mpr = new TH1F("h_mpr","MergePatRec Efficiency",4,0,4);
  TH1F* h_tpr = new TH1F("h_tpr","TrkPatRec   Efficiency",4,0,4);
  TH1F* hh    = h_tpr->Clone("hh");

  float dat_mpr[10], dat_tpr[10], mpr_over_tpr[10];
  float v1, v2, err;

  for (int i=0; i<kNDatasets; i++) {

    v1 = dat_mpr[i];
    v2 = dat_tpr[i];

    mpr_over_tpr[i] = v1/v2;

    err = 0.005;
      
    h_mpr->SetBinContent(i+1,v1);
    h_mpr->SetBinError  (i+1,err);
    h_tpr->SetBinContent(i+1,v2);
    hh->SetBinContent(i+1,v2);
    hh->SetBinError  (i+1,0.001);
  }

  for (int i=1; i<=kNDatasets; i++) {
    h_mpr->GetXaxis()->SetBinLabel(i,cut_label_cd2[i-1]);
  }

  h_mpr->SetTitle(Form("MergePatRec+CalPatRec efficiency"));
  h_mpr->SetStats(0);
  h_mpr->SetMaximum(0.16);
  h_mpr->SetMinimum(0.);
  h_mpr->SetMarkerStyle(24);
  h_mpr->SetMarkerSize(1.5);
  h_mpr->Draw("ep");
  h_mpr->Draw("same,text45");

  h_tpr->SetFillStyle(3001);
  h_tpr->SetFillColor(kBlue-7);
  h_tpr->Draw("same");

  hh->SetMarkerStyle(20);
  hh->SetMarkerSize(2.);
  hh->Draw("same,ep");


  TLegend* leg = new TLegend(0.7,0.75,0.9,0.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);

  leg->AddEntry(h_mpr,"TrkPatRec+CalPatRec","f");
  leg->AddEntry(h_tpr,"TrkPatRec","f");

  leg->Draw();
}

