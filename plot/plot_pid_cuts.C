///////////////////////////////////////////////////////////////////////////////
// use TrackAna module histograms
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
// 1. P
//-----------------------------------------------------------------------------
void plot_eff_p(const char* Fn) {

  h1 = gh1(Fn,"TrackAna","trk_1/p_2");
  h2 = gh1(Fn,"TrackAna","trk_13/p_2");

  TCanvas* c = new TCanvas("c_p","Track Momentum",0,0,1100,700);

  c->cd();
  gPad->SetLogy(1);

  h2->SetLineColor(kBlue-6);
  h2->SetFillColor(kBlue-6);
  h2->SetFillStyle(3001);

  h1->GetXaxis()->SetRangeUser(60,119.999);
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
// 2. E/P
//-----------------------------------------------------------------------------
void plot_eff_ep(const char* Fn) {

  h1 = gh1(Fn,"TrackAna","trk_13/ep");
  h2 = gh1(Fn,"TrackAna","trk_29/ep");

  TCanvas* c = new TCanvas("c_ep","Track Momentum",0,0,1100,700);

  c->cd();
  gPad->SetLogy(1);

  h2->SetLineColor(kBlue-6);
  h2->SetFillColor(kBlue-6);
  h2->SetFillStyle(3001);

  h1->GetXaxis()->SetRangeUser(60,119.999);
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
// 3. DT
//-----------------------------------------------------------------------------
void plot_eff_dt(const char* Fn) {

  h1 = gh1(Fn,"TrackAna","trk_29/dt");
  h2 = gh1(Fn,"TrackAna","trk_32/dt");

  TCanvas* c = new TCanvas("c_dt","Delta T",0,0,1100,700);

  c->cd();
  gPad->SetLogy(1);

  h2->SetLineColor(kBlue-6);
  h2->SetFillColor(kBlue-6);
  h2->SetFillStyle(3001);

  h1->GetXaxis()->SetRangeUser(60,119.999);
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
// 4. chi2(TCM)
//-----------------------------------------------------------------------------
void plot_eff_chi2tcm(const char* Fn) {

  h1 = gh1(Fn,"TrackAna","trk_32/chi2tcm");
  h2 = gh1(Fn,"TrackAna","trk_25/chi2tcm");

  TCanvas* c = new TCanvas("c_chi2tcm","Chi2 TCM",0,0,1100,700);

  c->cd();
  gPad->SetLogy(1);

  h2->SetLineColor(kBlue-6);
  h2->SetFillColor(kBlue-6);
  h2->SetFillStyle(3001);

  //  h1->GetXaxis()->SetRangeUser(60,119.999);
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



