//

namespace {
  char* cut_label[10] = {
    "DetFilter"   , "N(MC)>20", "Pf>100 MeV/c", "CE Pitch", "Track Fit", 
    "Fit Qual", "T0>700"  , "Reco Pitch", "cosmics", "103.5 < P < 105 MeV/c"
  };

};

//-----------------------------------------------------------------------------
// TrkAlg = "TrkPatRec", or "CalPatRec", or "MergePatRec"
//-----------------------------------------------------------------------------
void plot_daves_ladder(const char* Filename, const char* TrkAlg, int NEvents) {

  TH1F* h_mpr = new TH1F("h_mpr","MergePatRec Efficiency",10,0,10);
  TH1F* h_tpr = new TH1F("h_tpr","TrkPatRec   Efficiency",10,0,10);

  float dat_mpr[10], dat_tpr[10];	// for TrackCompAna
  int   iset;

  if      (strcmp(TrkAlg,"TrkPatRec") == 0) iset = 10;
  else if (strcmp(TrkAlg,"CalPatRec") == 0) iset = 20;
  
					// total N(events) in the dataset ( <= NEvents)
  dat_mpr[0] = gh1(Filename,"TrackComp","evt_0/ce_costh")->GetEntries();
  //  dat_tpr[0] = dat_mpr[0];
					// CE N(straw hits) > 20
  dat_mpr[1] = gh1(Filename,"TrackComp",Form("evt_%i/ce_costh",iset))->GetEntries();
  //  dat_tpr[1] = dat_mpr[1];
					// CE P > 100
  dat_mpr[2] = gh1(Filename,"TrackComp",Form("evt_%i/ce_costh",iset+1))->GetEntries();
  //  dat_tpr[2] = dat_mpr[2];
					// CE pitch
  dat_mpr[3] = gh1(Filename,"TrackComp",Form("evt_%i/ce_costh",iset+2))->GetEntries();
  //  dat_tpr[3] = dat_mpr[3];
					// commonality ends, track reconstructed and fit
  dat_mpr[4] = gh1(Filename,"TrackComp",Form("evt_%i/ce_costh",iset+3))->GetEntries();
  //  dat_tpr[4] = gh1(Filename,"TrackComp","evt_24/ce_costh")->GetEntries();

					// TrkQual
  dat_mpr[5] = gh1(Filename,"TrackComp",Form("evt_%i/ce_costh",iset+4))->GetEntries();
  //  dat_tpr[5] = gh1(Filename,"TrackComp","evt_25/ce_costh")->GetEntries();

  dat_mpr[6] = gh1(Filename,"TrackComp",Form("evt_%i/ce_costh",iset+5))->GetEntries();
  //  dat_tpr[6] = gh1(Filename,"TrackComp","evt_26/ce_costh")->GetEntries();

  dat_mpr[7] = gh1(Filename,"TrackComp",Form("evt_%i/ce_costh",iset+7))->GetEntries();
  //  dat_tpr[7] = gh1(Filename,"TrackComp","evt_27/ce_costh")->GetEntries();

  dat_mpr[8] = gh1(Filename,"TrackComp",Form("evt_%i/ce_costh",iset+8))->GetEntries();
  //  dat_tpr[8] = gh1(Filename,"TrackComp","evt_27/ce_costh")->GetEntries();

  dat_mpr[9] = gh1(Filename,"TrackComp",Form("evt_%i/ce_costh",iset+9))->GetEntries();
  //  dat_tpr[9] = gh1(Filename,"TrackComp","evt_27/ce_costh")->GetEntries();

  float v1, v2;

  for (int i=1; i<=10; i++) {
    h_mpr->GetXaxis()->SetBinLabel(i,cut_label[i-1]);

    v1 = dat_mpr[i-1]/(NEvents+0.);
    //    v2 = dat_tpr[i-1]/(dat_tpr[0]+0.);

    h_mpr->SetBinContent(i,v1);
    //    h_tpr->SetBinContent(i,v2);
  }

  h_mpr->SetTitle(Form("MergePatRec+CalPatRec efficiency: %s",Filename));
  h_mpr->SetStats(0);
  h_mpr->SetMinimum(0.);
  h_mpr->Draw();
  h_mpr->Draw("same,text45");

//   h_tpr->SetFillStyle(3003);
//   h_tpr->SetFillColor(kBlue-7);
//   h_tpr->Draw("same");


  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);

  leg->AddEntry(h_mpr,"TrkPatRec+CalPatRec","f");
  //  leg->AddEntry(h_tpr,"TrkPatRec","f");

  leg->Draw();
}


//-----------------------------------------------------------------------------
a2_mpr() {

  TH1F* h_cnv00202_mergepatrec = new TH1F("h_cnv00202_mergepatrec","Efficiency",10,0,10);


    float val;
    for (int i=1; i<=10; i++) {
    h_cnv00202_mergepatrec->GetXaxis()->SetBinLabel(i,cut_label[i-1]);

    if (i == 1) val = 1;
    else val = cnv00202_mergepatrec[i-1]/(cnv00202_mergepatrec[i-2]+0.);

    h_cnv00202_mergepatrec->SetBinContent(i,val);
  }

    h_cnv00202_mergepatrec->SetMaximum(1.5);
  h_cnv00202_mergepatrec->SetTitle("Mergepatrec+CalPatRec efficiency");
  h_cnv00202_mergepatrec->SetStats(0);
  h_cnv00202_mergepatrec->SetMinimum(0.);
  h_cnv00202_mergepatrec->Draw();
  h_cnv00202_mergepatrec->Draw("same,text45");


  TH1F* h_cnv00302_mergepatrec = new TH1F("h_cnv00302_mergepatrec","Efficiency",10,0,10);


  for (int i=1; i<=10; i++) {
    h_cnv00302_mergepatrec->GetXaxis()->SetBinLabel(i,cut_label[i-1]);

    if (i == 1) val = 1;
    else val = cnv00302_mergepatrec[i-1]/(cnv00302_mergepatrec[i-2]+0.);

    h_cnv00302_mergepatrec->SetBinContent(i,val);
  }

  h_cnv00302_mergepatrec->SetFillStyle(3003);
  h_cnv00302_mergepatrec->SetFillColor(kBlue-7);
  h_cnv00302_mergepatrec->Draw("same");


  TH1F* h_cnv00402_mergepatrec = new TH1F("h_cnv00402_mergepatrec","Efficiency",10,0,10);


  for (int i=1; i<=10; i++) {
    h_cnv00402_mergepatrec->GetXaxis()->SetBinLabel(i,cut_label[i-1]);

    if (i == 1) val = 1;
    else val = cnv00402_mergepatrec[i-1]/(cnv00402_mergepatrec[i-2]+0.);

    h_cnv00402_mergepatrec->SetBinContent(i,val);
  }
  h_cnv00402_mergepatrec->SetFillStyle(3001);
  h_cnv00402_mergepatrec->SetFillColor(kBlue-7);
  h_cnv00402_mergepatrec->Draw("same");


  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);

  leg->AddEntry(h_cnv00202_mergepatrec,"CE only","f");
  leg->AddEntry(h_cnv00302_mergepatrec,"CE+MIXP2","f");
  leg->AddEntry(h_cnv00402_mergepatrec,"CE+MIXP3","f");

  leg->Draw();
}
