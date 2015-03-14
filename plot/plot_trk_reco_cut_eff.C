///////////////////////////////////////////////////////////////////////////////
// use StrawHitAna module to study the track reco efficiency
///////////////////////////////////////////////////////////////////////////////

namespace {
  TH1F      *h1, *h2;
  TCanvas   *c;
};

namespace {
  char* cut_label[10] = {
    "Total"   , "N(MC)>20", "Pf>100 MeV/c", "CE Pitch", "Track Fit", 
    "Fit Qual", "T0>700"  , "Reco Pitch", "cosmics", "P > 103.5 MeV/c"
  };

};

//-----------------------------------------------------------------------------
// 1. N(CE hits) >= 20
//-----------------------------------------------------------------------------
void plot_eff_nhs(const char* Fn) {

  h1 = gh1(Fn,"StrawHitAna","evt_0/nhs");
  h2 = gh1(Fn,"StrawHitAna","evt_1/nhs");

  TCanvas* c = new TCanvas("c_nhs","N Hits Signal",0,0,1000,700);

  c->cd();
  gPad->SetLogy(1);

  h2->SetLineColor(kBlue-6);
  h2->SetFillColor(kBlue-6);
  h2->SetFillStyle(3001);

  h1->Draw();

  h2->Draw("sames");

  c->Modified();
  c->Update();
  
  TPaveStats* s1 = (TPaveStats*) h1->GetListOfFunctions()->FindObject("stats");

  s1->SetY1NDC(0.65);
  s1->SetY2NDC(0.90);
  s1->Draw();

  TPaveStats* s2 = (TPaveStats*) h2->GetListOfFunctions()->FindObject("stats");

  s2->SetY1NDC(0.40);
  s2->SetY2NDC(0.65);
  s2->Draw();
}

//-----------------------------------------------------------------------------
// 2. MC Momentum (Tracker Front) > 100.
//-----------------------------------------------------------------------------
void plot_eff_mom_tf(const char* Fn) {

  h1 = gh1(Fn,"StrawHitAna","evt_1/mom_tf");
  h2 = gh1(Fn,"StrawHitAna","evt_2/mom_tf");

  TCanvas* c = new TCanvas("c_mom_tf","Momentum (TF) > 100 MeV/c",0,0,1000,700);

  c->cd();
  gPad->SetLogy(1);

  h2->SetLineColor(kBlue-6);
  h2->SetFillColor(kBlue-6);
  h2->SetFillStyle(3001);

  h1->Draw();

  h2->Draw("sames");

  c->Modified();
  c->Update();
  
  TPaveStats* s1 = (TPaveStats*) h1->GetListOfFunctions()->FindObject("stats");

  s1->SetX1NDC(0.1);
  s1->SetX2NDC(0.30);
  s1->SetY1NDC(0.65);
  s1->SetY2NDC(0.90);
  s1->Draw();

  TPaveStats* s2 = (TPaveStats*) h2->GetListOfFunctions()->FindObject("stats");

  s2->SetX1NDC(0.1);
  s2->SetX2NDC(0.30);
  s2->SetY1NDC(0.40);
  s2->SetY2NDC(0.65);
  s2->Draw();
}

//-----------------------------------------------------------------------------
// 3. < tan(dip) MC < sqrt(3.)
//-----------------------------------------------------------------------------
void plot_eff_tdip_mc(const char* Fn) {

  h1 = gh1(Fn,"StrawHitAna","evt_2/tdip_mc");
  h2 = gh1(Fn,"StrawHitAna","evt_3/tdip_mc");

  TCanvas* c = new TCanvas("c_tdip_mc","1 < tan(dip) < sqrt(3)",0,0,1000,700);

  c->cd();
  //  gPad->SetLogy(1);

  h2->SetLineColor(kBlue-6);
  h2->SetFillColor(kBlue-6);
  h2->SetFillStyle(3001);

  h1->Draw();

  h2->Draw("sames");

  c->Modified();
  c->Update();
  
  TPaveStats* s1 = (TPaveStats*) h1->GetListOfFunctions()->FindObject("stats");

  s1->SetX1NDC(0.1);
  s1->SetX2NDC(0.30);
  s1->SetY1NDC(0.65);
  s1->SetY2NDC(0.90);
  s1->Draw();

  TPaveStats* s2 = (TPaveStats*) h2->GetListOfFunctions()->FindObject("stats");

  s2->SetX1NDC(0.1);
  s2->SetX2NDC(0.30);
  s2->SetY1NDC(0.40);
  s2->SetY2NDC(0.65);
  s2->Draw();
}


//-----------------------------------------------------------------------------
// 4. N(reconstructed tracks) > 0
//-----------------------------------------------------------------------------
void plot_eff_ntrk(const char* Fn) {

  h1 = gh1(Fn,"StrawHitAna","evt_3/ntrk");
  h2 = gh1(Fn,"StrawHitAna","evt_4/ntrk");

  TCanvas* c = new TCanvas("c_ntrk","N(tracks) > 0",0,0,1000,700);

  c->cd();
  gPad->SetLogy(1);

  h2->SetLineColor(kBlue-6);
  h2->SetFillColor(kBlue-6);
  h2->SetFillStyle(3001);

  h1->GetXaxis()->SetRangeUser(0.,4.99);
  h1->Draw();

  h2->Draw("sames");

  c->Modified();
  c->Update();
  
  TPaveStats* s1 = (TPaveStats*) h1->GetListOfFunctions()->FindObject("stats");

  s1->SetX1NDC(0.70);
  s1->SetX2NDC(0.90);
  s1->SetY1NDC(0.65);
  s1->SetY2NDC(0.90);
  s1->Draw();

  TPaveStats* s2 = (TPaveStats*) h2->GetListOfFunctions()->FindObject("stats");

  s2->SetX1NDC(0.70);
  s2->SetX2NDC(0.90);
  s2->SetY1NDC(0.40);
  s2->SetY2NDC(0.65);
  s2->Draw();
}


//-----------------------------------------------------------------------------
// 5. Track quality cuts passed
//-----------------------------------------------------------------------------
void plot_eff_track_quality(const char* Fn) {

  h1 = gh1(Fn,"StrawHitAna","evt_4/ntrk");
  h2 = gh1(Fn,"StrawHitAna","evt_5/ntrk");

  TCanvas* c = new TCanvas("c_quality_cuts","Track Quality Cuts Passed",0,0,1000,700);

  c->cd();
  gPad->SetLogy(1);

  h2->SetLineColor(kBlue-6);
  h2->SetFillColor(kBlue-6);
  h2->SetFillStyle(3001);

  h1->GetXaxis()->SetRangeUser(0.,4.99);
  h1->Draw();

  h2->Draw("sames");

  c->Modified();
  c->Update();
  
  TPaveStats* s1 = (TPaveStats*) h1->GetListOfFunctions()->FindObject("stats");

  s1->SetX1NDC(0.70);
  s1->SetX2NDC(0.90);
  s1->SetY1NDC(0.65);
  s1->SetY2NDC(0.90);
  s1->Draw();

  TPaveStats* s2 = (TPaveStats*) h2->GetListOfFunctions()->FindObject("stats");

  s2->SetX1NDC(0.70);
  s2->SetX2NDC(0.90);
  s2->SetY1NDC(0.40);
  s2->SetY2NDC(0.65);
  s2->Draw();
}

//-----------------------------------------------------------------------------
// 6. Track T0 > 700 ns
//-----------------------------------------------------------------------------
void plot_eff_t0(const char* Fn) {

  h1 = gh1(Fn,"StrawHitAna","evt_5/t0");
  h2 = gh1(Fn,"StrawHitAna","evt_6/t0");

  TCanvas* c = new TCanvas("c_t0","Track T0",0,0,1000,700);

  c->cd();
  //  gPad->SetLogy(1);

  h2->SetLineColor(kBlue-6);
  h2->SetFillColor(kBlue-6);
  h2->SetFillStyle(3001);

  //  h1->GetXaxis()->SetRangeUser(0.,4.99);
  h1->Draw();

  h2->Draw("sames");

  c->Modified();
  c->Update();
  
  TPaveStats* s1 = (TPaveStats*) h1->GetListOfFunctions()->FindObject("stats");

  s1->SetX1NDC(0.70);
  s1->SetX2NDC(0.90);
  s1->SetY1NDC(0.65);
  s1->SetY2NDC(0.90);
  s1->Draw();

  TPaveStats* s2 = (TPaveStats*) h2->GetListOfFunctions()->FindObject("stats");

  s2->SetX1NDC(0.70);
  s2->SetX2NDC(0.90);
  s2->SetY1NDC(0.40);
  s2->SetY2NDC(0.65);
  s2->Draw();
}

//-----------------------------------------------------------------------------
// 7. 1. < Reco pitch < sqrt(3.)
//-----------------------------------------------------------------------------
void plot_eff_reco_tdip(const char* Fn) {

  h1 = gh1(Fn,"StrawHitAna","evt_6/tdip");
  h2 = gh1(Fn,"StrawHitAna","evt_7/tdip");

  TCanvas* c = new TCanvas("c_tdip","1 < Reco(Pitch) < sqrt(3.)",0,0,1000,700);

  c->cd();
  //  gPad->SetLogy(1);

  h2->SetLineColor(kBlue-6);
  h2->SetFillColor(kBlue-6);
  h2->SetFillStyle(3001);

  h1->GetXaxis()->SetRangeUser(0.6,1.99);
  h1->Draw();

  h2->Draw("sames");

  c->Modified();
  c->Update();
  
  TPaveStats* s1 = (TPaveStats*) h1->GetListOfFunctions()->FindObject("stats");

  s1->SetX1NDC(0.70);
  s1->SetX2NDC(0.90);
  s1->SetY1NDC(0.65);
  s1->SetY2NDC(0.90);
  s1->Draw();

  TPaveStats* s2 = (TPaveStats*) h2->GetListOfFunctions()->FindObject("stats");

  s2->SetX1NDC(0.70);
  s2->SetX2NDC(0.90);
  s2->SetY1NDC(0.40);
  s2->SetY2NDC(0.65);
  s2->Draw();
}

//-----------------------------------------------------------------------------
// 8. cosmics rejection
//-----------------------------------------------------------------------------
void plot_eff_cosmics(const char* Fn) {

  h1 = gh1(Fn,"StrawHitAna","evt_7/ntrk");
  h2 = gh1(Fn,"StrawHitAna","evt_8/ntrk");

  TCanvas* c = new TCanvas("c_d0","cosmics: Rmin < A, RMax < B",0,0,1000,700);

  c->cd();
  gPad->SetLogy(1);

  h2->SetLineColor(kBlue-6);
  h2->SetFillColor(kBlue-6);
  h2->SetFillStyle(3001);

  h1->GetXaxis()->SetRangeUser(0.6,4.99);
  h1->Draw();

  h2->Draw("sames");

  c->Modified();
  c->Update();
  
  TPaveStats* s1 = (TPaveStats*) h1->GetListOfFunctions()->FindObject("stats");

  s1->SetX1NDC(0.70);
  s1->SetX2NDC(0.90);
  s1->SetY1NDC(0.65);
  s1->SetY2NDC(0.90);
  s1->Draw();

  TPaveStats* s2 = (TPaveStats*) h2->GetListOfFunctions()->FindObject("stats");

  s2->SetX1NDC(0.70);
  s2->SetX2NDC(0.90);
  s2->SetY1NDC(0.40);
  s2->SetY2NDC(0.65);
  s2->Draw();
}

//-----------------------------------------------------------------------------
// 9. sigmal momentum window
//-----------------------------------------------------------------------------
void plot_eff_mom(const char* Fn) {

  h1 = gh1(Fn,"StrawHitAna","evt_8/p");
  h2 = gh1(Fn,"StrawHitAna","evt_9/p");

  TCanvas* c = new TCanvas("c_p","10.5 < Reconstructed momentum < 105.0",0,0,1000,700);

  c->cd();
  gPad->SetLogy(1);

  h2->SetLineColor(kBlue-6);
  h2->SetFillColor(kBlue-6);
  h2->SetFillStyle(3001);

  h1->GetXaxis()->SetRangeUser(60,109.99);
  h1->Draw();

  h2->Draw("sames");

  c->Modified();
  c->Update();
  
  TPaveStats* s1 = (TPaveStats*) h1->GetListOfFunctions()->FindObject("stats");

  s1->SetX1NDC(0.10);
  s1->SetX2NDC(0.30);
  s1->SetY1NDC(0.65);
  s1->SetY2NDC(0.90);
  s1->Draw();

  TPaveStats* s2 = (TPaveStats*) h2->GetListOfFunctions()->FindObject("stats");

  s2->SetX1NDC(0.10);
  s2->SetX2NDC(0.30);
  s2->SetY1NDC(0.40);
  s2->SetY2NDC(0.65);
  s2->Draw();
}


//-----------------------------------------------------------------------------
void plot_ladder(const char* Filename) {


  TH1F* h_mpr = new TH1F("h_mpr","MergePatRec Efficiency",10,0,10);
  //  TH1F* h_tpr = new TH1F("h_tpr","TrkPatRec   Efficiency",10,0,10);

  float dat_mpr[10], dat_tpr[10];


  dat_mpr[0] = gh1(Filename,"StrawHitAna","evt_0/nhs")->GetEntries();
  dat_mpr[1] = gh1(Filename,"StrawHitAna","evt_1/nhs")->GetEntries();
  dat_mpr[2] = gh1(Filename,"StrawHitAna","evt_2/nhs")->GetEntries();
  dat_mpr[3] = gh1(Filename,"StrawHitAna","evt_3/nhs")->GetEntries();
  dat_mpr[4] = gh1(Filename,"StrawHitAna","evt_4/nhs")->GetEntries();
  dat_mpr[5] = gh1(Filename,"StrawHitAna","evt_5/nhs")->GetEntries();
  dat_mpr[6] = gh1(Filename,"StrawHitAna","evt_6/nhs")->GetEntries();
  dat_mpr[7] = gh1(Filename,"StrawHitAna","evt_7/nhs")->GetEntries();
  dat_mpr[8] = gh1(Filename,"StrawHitAna","evt_8/nhs")->GetEntries();
  dat_mpr[9] = gh1(Filename,"StrawHitAna","evt_9/nhs")->GetEntries();


  float v1, v2;

  for (int i=1; i<=10; i++) {
    h_mpr->GetXaxis()->SetBinLabel(i,cut_label[i-1]);

    v1 = dat_mpr[i-1]/(dat_mpr[0]+0.);

    h_mpr->SetBinContent(i,v1);
  }

  h_mpr->SetTitle(Form("MergePatRec+CalPatRec efficiency: %s",Filename));
  h_mpr->SetStats(0);
  h_mpr->SetMinimum(0.);
  h_mpr->Draw();
  h_mpr->Draw("same,text45");

  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);

  leg->AddEntry(h_mpr,"TrkPatRec+CalPatRec","f");

  leg->Draw();
}


