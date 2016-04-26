///////////////////////////////////////////////////////////////////////////////
// fit efficiency vs straw hit emax
// data come from e31s5730 dataset, each fileset of which is 100K events 
// generated with a given EMAX threshold (000002 = emax = 0.002 MeV)
///////////////////////////////////////////////////////////////////////////////

#include "murat/scripts/datasets.hh"
#include "murat/scripts/fitters.hh"

const char* figures_dir =  ".";

//-----------------------------------------------------------------------------
// redraw the stat box, witht eh color of the histogram itself
//-----------------------------------------------------------------------------
void plot_stat_box(TH1* Hist, double X1, double Y1, double X2, double Y2) {
  TPaveStats* s = (TPaveStats*) Hist->GetListOfFunctions()->FindObject("stats");
  s->SetLineColor(Hist->GetLineColor());
  s->SetTextColor(Hist->GetLineColor());
  s->SetX1NDC(X1); s->SetY1NDC(Y1); s->SetX2NDC(X2); s->SetY2NDC(Y2);
  s->Draw();
}

//-----------------------------------------------------------------------------
void draw_label_ndc(const char* Text, double X1, double Y1, double FontSize, double Font = 52) {
  TLatex* label = new TLatex(X1,Y1,Text);
  label->SetNDC();
  label->SetTextSize(0.035);
  label->SetTextFont(52);
  label->Draw();
}

//-----------------------------------------------------------------------------
double f(double* X, double* P) {
  double x = X[0];

  double fun = P[0]/(1+(P[2]*P[2])*(x-P[1])*(x-P[1])/(x));

  return fun;
}

//-----------------------------------------------------------------------------
void fit_eff_vs_emax(int Print = 0) {

  int const npt(8);

  double emax [npt] = { 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009 };

  double eff  [npt] = { 0.10330, 0.11036, 0.11265, 0.11200, 0.11246, 0.11075, 0.10913, 0.10855};

  double ex[npt], ey[npt];


  for (int i=0; i<npt; i++) {
    ex[i] = 0;
    ey[i] = eff[i]/100.;
  }

  TCanvas* c = new TCanvas("c_fit","fit eff vs emax",1400,800);

  TH1F* h = new TH1F("h1","",100,0,0.01);
  h->SetStats(0);
  h->GetXaxis()->SetTitle("straw hit E(max), MeV");
  h->GetYaxis()->SetRangeUser(0.10,0.12);
  h->Draw();
  
  TGraphErrors* gr = new TGraphErrors(npt,emax,eff,ex,ey);

  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(1.0);
  
  gr->Draw("same,p");

  TF1* f = new TF1("f",f,0.001,0.01,3);

  f->SetParameters(0.11,0.005,0.1);

  gr->Fit(f,"","",0.002,0.009);

  draw_label_ndc("TrackAna/trk_33/p efficiency (103,105) MeV",0.15,0.85,20,52);

  if (Print == 1) c->Print(Form("%s/%s",figures_dir,Form("efficiency_vs_straw_hit_emax_v573.eps")));
}
