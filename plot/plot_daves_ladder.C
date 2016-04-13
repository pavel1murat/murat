//

namespace {
  char* cut_label[10] = {
    "DetFilter"   , "Nsh(CE)>20", "P(CE TF)>100 MeV/c", "CE Pitch", "Track Fit", 
    "Fit Qual", "T0>700"  , "Reco Pitch", "cosmics", "103.5 < P < 105 MeV/c"
  };

  const char* e42s5721_track_ana  = "~/hist/mu2e/v5_7_0/e42s5721.track_ana.hist" ; // CE+BGR  , matcorr in CalPatRec
  const char* e42s5721_track_comp = "~/hist/mu2e/v5_7_0/e42s5721.track_comp.hist"; // CE+BGR  , matcorr in CalPatRec

};


//-----------------------------------------------------------------------------
void make_eff_ladder_hist(const char* Filename, const char* TrkAlg, int NEvents, TH1F*& Hist) {

  if (Hist != 0) delete Hist;
  
  Hist = new TH1F(Form("h_eff_ladder_%s",TrkAlg),"Efficiency",10,0,10);
  
  float dat[10];	// for TrackCompAna circa Apr'2016
  int   iset;

  if      (strcmp(TrkAlg,"TrkPatRec") == 0) iset = 10;
  else if (strcmp(TrkAlg,"CalPatRec") == 0) iset = 20;
  
  // total N(events) in the dataset ( <= NEvents)
  dat[0] = gh1(Filename,"TrackComp","evt_0/ce_costh")->GetEntries();

  // CE N(straw hits) > 20
  dat[1] = gh1(Filename,"TrackComp",Form("evt_%i/ce_costh",iset))->GetEntries();

  // CE P > 100
  dat[2] = gh1(Filename,"TrackComp",Form("evt_%i/ce_costh",iset+1))->GetEntries();

  // CE pitch
  dat[3] = gh1(Filename,"TrackComp",Form("evt_%i/ce_costh",iset+2))->GetEntries();

  // commonality ends, track reconstructed and fit
  dat[4] = gh1(Filename,"TrackComp",Form("evt_%i/ce_costh",iset+3))->GetEntries();

  // TrkQual
  dat[5] = gh1(Filename,"TrackComp",Form("evt_%i/ce_costh",iset+4))->GetEntries();

  dat[6] = gh1(Filename,"TrackComp",Form("evt_%i/ce_costh",iset+5))->GetEntries();

  dat[7] = gh1(Filename,"TrackComp",Form("evt_%i/ce_costh",iset+7))->GetEntries();

  dat[8] = gh1(Filename,"TrackComp",Form("evt_%i/ce_costh",iset+8))->GetEntries();

  dat[9] = gh1(Filename,"TrackComp",Form("evt_%i/ce_costh",iset+9))->GetEntries();

  for (int i=1; i<=10; i++) {
    Hist->GetXaxis()->SetBinLabel(i,cut_label[i-1]);
    float v1 = dat[i-1]/(NEvents+0.);
    Hist->SetBinContent(i,v1);
  }
}

//-----------------------------------------------------------------------------
// TrkAlg = "TrkPatRec", or "CalPatRec", or "MergePatRec"
//-----------------------------------------------------------------------------
void plot_daves_ladder(const char* Filename, const char* TrkAlg, int NEvents) {

  TH1F* hist(0);
  
  make_eff_ladder_hist(Filename,TrkAlg,NEvents,hist);

  hist->SetTitle(Form("%s efficiency: %s",TrkAlg,Filename));
  hist->SetStats(0);
  hist->SetMaximum(0.8);
  hist->SetMinimum(0.);
  hist->Draw();
  hist->Draw("same,text45");

  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);

  leg->AddEntry(hist,TrkAlg,"f");

  leg->Draw();
}


//-----------------------------------------------------------------------------
// TrkAlg = "TrkPatRec", or "CalPatRec", or "MergePatRec"
//-----------------------------------------------------------------------------
void plot_trkpatrec_vs_calpatrec_ladders(const char* Filename, int NEvents) {

  TH1F  *h_tpr(0), *h_cpr(0);
  
  make_eff_ladder_hist(Filename,"TrkPatRec",NEvents,h_tpr);
  make_eff_ladder_hist(Filename,"CalPatRec",NEvents,h_cpr);

  h_tpr->SetTitle("");
  h_tpr->SetStats(0);
  h_tpr->SetMaximum(0.8);
  h_tpr->SetMinimum(0.);
  h_tpr->SetLineColor(2);
  h_tpr->SetLineWidth(2);
  h_tpr->Draw();
  h_tpr->Draw("same,text45");

  h_cpr->Draw("same");

  TText* label = new TText(1.5,0.75,Form("efficiency, TrkPatRec vs CalPatRec: %s",Filename));
  label->SetTextSize(0.035);
  label->SetTextFont(52);
  label->Draw();
  
  TLegend* leg = new TLegend(0.7,0.6,0.9,0.75);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);

  leg->AddEntry(h_tpr,"TrkPatRec","f");
  leg->AddEntry(h_cpr,"CalPatRec","f");

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
