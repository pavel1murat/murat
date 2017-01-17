///////////////////////////////////////////////////////////////////////////////
// use TrackAna module histograms
///////////////////////////////////////////////////////////////////////////////
#include <string>

namespace {
  TH1F      *h1, *h2;
  TCanvas   *c;
  
  std::string figures_dir = ".";
};

namespace {
  std::string cut_label[10] = {
    "Total"   , "N(MC)>20", "Pf>100 MeV/c", "CE Pitch", "Track Fit", 
    "Fit Qual", "T0>700"  , "Reco Pitch", "cosmics", "P > 103.5 MeV/c"
  };

};

//-----------------------------------------------------------------------------
// redraw the stat box, use the histogram line color
//-----------------------------------------------------------------------------
void plot_stat_box(TH1* Hist, double X1, double Y1, double X2, double Y2) {
  TPaveStats* s = (TPaveStats*) Hist->GetListOfFunctions()->FindObject("stats");
  if (s != NULL) {
    s->SetLineColor(Hist->GetLineColor());
    s->SetTextColor(Hist->GetLineColor());
    s->SetX1NDC(X1);
    s->SetY1NDC(Y1);
    s->SetX2NDC(X2);
    s->SetY2NDC(Y2);
    s->Draw();
  }
  else {
    printf("ERROR: stat box is not defined for %s\n",Hist->GetName());
  }
}

//-----------------------------------------------------------------------------
// 1. P
//-----------------------------------------------------------------------------
void plot_eff_p(const char* Fn, int Print = 0) {

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

  printf(" efficiency = %10.5f\n",h2->Integral()/h1->GetEntries());
//-----------------------------------------------------------------------------
// .eps files are written into figures/eps subdirectory
//-----------------------------------------------------------------------------
  if (Print == 1) c->Print(Form("%s/eps/pid_cut_eff_p.eps",figures_dir.data()));
}

//-----------------------------------------------------------------------------
// 2. 0.1 < E/P < 1.15
//-----------------------------------------------------------------------------
void plot_eff_ep(const char* Fn, int Print = 0) {

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
  
  printf(" efficiency = %10.5f\n",h2->Integral()/h1->GetEntries());
//-----------------------------------------------------------------------------
// .eps files are written into figures/eps subdirectory
//-----------------------------------------------------------------------------
  if (Print == 1) c->Print(Form("%s/eps/pid_cut_eff_ep.eps",figures_dir.data()));
}


//-----------------------------------------------------------------------------
// 3. chi2(TCM)
//-----------------------------------------------------------------------------
void plot_eff_chi2tcm(const char* Fn, int Print = 0) {

  h1 = gh1(Fn,"TrackAna","trk_29/chi2tcm");   // P+E/P cuts
  h2 = gh1(Fn,"TrackAna","trk_32/chi2tcm");   // P+E/P+tcmchi2 cuts

  TCanvas* c = new TCanvas("c_chi2tcm","Chi2(TCM)",0,0,1100,700);

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
  
  printf(" efficiency = %10.5f\n",h2->Integral()/h1->GetEntries());
//-----------------------------------------------------------------------------
// .eps files are written into figures/eps subdirectory
//-----------------------------------------------------------------------------
  if (Print == 1) c->Print(Form("%s/eps/pid_cut_eff_chi2tcm.eps",figures_dir.data()));
}


//-----------------------------------------------------------------------------
// 4. DT
//-----------------------------------------------------------------------------
void plot_eff_dt(const char* Fn, int Print = 0) {

  h1 = gh1(Fn,"TrackAna","trk_32/dt"); // P+E/P+TCMCHI2 cuts
  h2 = gh1(Fn,"TrackAna","trk_25/dt"); // P+E/P+TCMCHI2+DT cuts

  TCanvas* c = new TCanvas("c_dt","Chi2 TCM",0,0,1100,700);

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
  
  printf(" efficiency = %10.5f\n",h2->Integral()/h1->GetEntries());
//-----------------------------------------------------------------------------
// .eps files are written into figures/eps subdirectory
//-----------------------------------------------------------------------------
  if (Print == 1) c->Print(Form("%s/eps/pid_cut_eff_dt.eps",figures_dir.data()));
}



