// -*- mode:c++ -*-
// it makes sens to have test script names including teh name of the tested class explicitly
//

#include "TCanvas.h"
#include "TRandom3.h"
#include "TH1D.h"
#include "TMath.h"

#include "murat/alg/TLogLHR.hh"

TLogLHR*          g_llhr(nullptr);

TH1D *h_lhr_001_pu, *h_lhr_001_pb;
//-----------------------------------------------------------------------------
// show how a measurement of N events changes the background expectation
//-----------------------------------------------------------------------------
void lhr_001(double CL, double MuB, double  MuS, double NObs) {
  TLogLHR*   llhr(nullptr);
  
  if (llhr == nullptr) llhr = new TLogLHR("lhr_005",CL);
  else                 llhr->SetCL(CL);

  llhr->fDebugLevel = 1;
  
  int nmax = 50;
  double pu[nmax], pb[nmax];

  llhr->InitPoissonDist(MuB,MuS,  -1,pu,nmax);
  llhr->InitPoissonDist(MuB,MuS,NObs,pb,nmax);

  //-----------------------------------------------------------------------------
  // plottting
  //-----------------------------------------------------------------------------
  h_lhr_001_pu = new TH1D("h_pu","h p(unbiased)",nmax,0,nmax);
  h_lhr_001_pb = new TH1D("h_pb","p(biased)"    ,nmax,0,nmax);

  printf(" lhr_001: booked hists\n");

  for (int i=0; i<nmax; i++) {
    h_lhr_001_pu->SetBinContent(i+1,pu[i]);
    h_lhr_001_pu->SetBinError  (i+1,0.001);
    h_lhr_001_pb->SetBinContent(i+1,pb[i]);
    h_lhr_001_pb->SetBinError  (i+1,0.001);
  }

  printf(" lhr_001: filled hists\n");

  h_lhr_001_pu->SetMarkerStyle(20);
  h_lhr_001_pu->SetMarkerColor(kRed);
  h_lhr_001_pb->SetMarkerStyle(20);
  h_lhr_001_pb->SetMarkerColor(kBlue+2);

  TCanvas* c = new TCanvas("c_lhr_001","LLHR 001",1100,700);

  h_lhr_001_pu->Draw();
  h_lhr_001_pb->Draw("sames");
}

//-----------------------------------------------------------------------------
// for "5-sigma" consistency, set CL = -1
//-----------------------------------------------------------------------------
void lhr_005(double CL, double B=0.1, double  SMin = 0, double SMax = 0, double NPoints = 1) {
  TLogLHR*   llhr;
  
  if (llhr == nullptr) llhr = new TLogLHR("lhr_005",CL);
  else                 llhr->SetCL(CL);

  g_llhr = llhr;
  
  double prob[10000], mean_llhr[10000];

  llhr->ConfInterval(B,SMin,SMax,NPoints,prob,mean_llhr);

  TH1D *h_prob, *h_llhr;

  double step;
  if (NPoints > 1) {
    step   = (SMax-SMin)/(NPoints-1);
    h_prob = new TH1D("h_prob","lhr_005 Prob"  ,NPoints,SMin-step/2,SMax+step/2);
    h_llhr = new TH1D("h_llhr","lhr_005 <LLHR>",NPoints,SMin-step/2,SMax+step/2);
  }
  else {
    h_prob = new TH1D("h_prob","lhr_005 Prob",  1, SMin-0.5, SMin+0.5);
    h_llhr = new TH1D("h_llhr","lhr_005 <LLHR>",1, SMin-0.5, SMin+0.5);
  }
  
  for (int i=0; i<NPoints; i++) {
    printf("i, prob[i] : %3i %12.5e\n",i, prob[i]);
    h_prob->SetBinContent(i+1,prob[i]);
    h_llhr->SetBinContent(i+1,mean_llhr[i]);
  }

  TCanvas* c = new TCanvas("c_lhr_005","FC 005",1400,700);
  c->Divide(2,1);
  if (! c->GetShowEventStatus()) c->ToggleEventStatus();

  c->cd(1);
  // gPad->SetLogy(1);
  gPad->SetGridy(1);
  gPad->SetCrosshair(1);

  h_prob->GetYaxis()->SetRangeUser(1.e-2,1);
  h_prob->Draw();
  
  c->cd(2);
  h_llhr->Draw();
}

//-----------------------------------------------------------------------------
// experiment with the expected background  N
// scan a range of signal hypotheses [SMin,SMax] with a step 'Step'
// determine CL confidence interval around the measurement of N events
//-----------------------------------------------------------------------------
void lhr_006(double CL, double MuB=0.1, double N = 5, double  SMin = 0, double SMax = 0, int NPoints = 1) {
  TLogLHR*   llhr(nullptr);

  if (llhr == nullptr) llhr = new TLogLHR("lhr_006",CL);
  else                 llhr->SetCL(CL);

  g_llhr = llhr;
  
  double prob[NPoints], mean_llhr[NPoints];

  llhr->MeasInterval(MuB,N,SMin,SMax,NPoints,prob,mean_llhr);

  TH1D *h_prob, *h_llhr;
  
  double step;
  if (NPoints > 1) {
    step = (SMax-SMin)/(NPoints-1);
    // h_prob = new TH1D("h_prob","lhr_006 Prob"  ,NPoints,SMin-step/2,SMax+step/2);
    h_llhr = new TH1D("h_llhr","lhr_006 <LLHR>",NPoints,SMin-step/2,SMax+step/2);
  }
  else {
    // h_prob = new TH1D("h_prob","lhr_006 Prob",  1, SMin-0.5, SMin+0.5);
    h_llhr = new TH1D("h_llhr","lhr_006 <LLHR>",1, SMin-0.5, SMin+0.5);
  }
  
  for (int i=0; i<NPoints; i++) {
    // printf("i, prob[i], <llhr[i]> : %3i %12.5e %12.5e\n",i, prob[i],mean_llhr[i]);
    // h_prob->SetBinContent(i+1,prob     [i]);
    h_llhr->SetBinContent(i+1,mean_llhr[i]);
  }

  TCanvas* c = new TCanvas("c_lhr_006","FC 006",1400,600);
  //  gPad->SetLogy(1);
  // c->Divide(2,1);
  if (! c->GetShowEventStatus()) c->ToggleEventStatus();
  c->SetCrosshair(1);

  // c->cd(1);
  // gPad->SetGridx(1);
  // gPad->SetGridy(1);
  // h_prob->GetYaxis()->SetRangeUser(1.e-5,1);
  // h_prob->Draw();

  // c->cd(2);
  h_llhr->SetTitle(Form("CL:%5.3f MuB: %6.3f Nobs:%3f SMin: %6.3f SMax: %6.3f NPoints:%3i",
			CL,MuB,N,SMin,SMax,NPoints));
  
  h_llhr->GetYaxis()->SetRangeUser(-20,1);
  h_llhr->Draw();
}

//-----------------------------------------------------------------------------
// experiment with the expected background  N
// scan a range of signal hypotheses [SMin,SMax] with a step 'Step'
// determine CLb propability of a CL-level observation of a signal
//-----------------------------------------------------------------------------
void lhr_007(double CL, double MuB=0.1, double  SMin = 0, double SMax = 0, int NPoints = 1) {
  TLogLHR*   llhr(nullptr);

  if (llhr == nullptr) llhr = new TLogLHR("lhr_007",CL);
  else                 llhr->SetCL(CL);

  g_llhr = llhr;
  
  double  x[NPoints], prob[NPoints];
  
  llhr->DiscoveryProbCLb(MuB,SMin,SMax,NPoints,x,prob);

  TCanvas* c = new TCanvas("c_lhr_007","FC 006",1100,800);
  if (! c->GetShowEventStatus()) c->ToggleEventStatus();
  gPad->SetCrosshair(1);

  TGraph* gr = new TGraph(NPoints,x,prob);

  gr->SetName(Form("gr_bgr_%06i_lhr_007",int(MuB*1000)));
  gr->SetTitle(Form("CLb discovery prob for bgr=%5.3f events",MuB));

  gr->SetMarkerStyle(20);
  gr->Draw("alp");

  gPad->SetGridy(1);
  gPad->SetGridx(1);

  gPad->Modified();
  gPad->Update();
}

//-----------------------------------------------------------------------------
// experiment with the expected background  N
// scan a range of signal hypotheses [SMin,SMax] with a step 'Step'
// determine CLs propability of a CL-level observation of a signal
//-----------------------------------------------------------------------------
void lhr_008(double CL, double MuB=0.1, double  SMin = 0, double SMax = 0, int NPoints = 1) {
  TLogLHR*   llhr;

  if (llhr == nullptr) llhr = new TLogLHR("lhr_007",CL);
  else                 llhr->SetCL(CL);

  g_llhr = llhr;
  
  double  x[NPoints], prob[NPoints];
  
  llhr->DiscoveryProbCLs(MuB,SMin,SMax,NPoints,x,prob);

  double step = (NPoints > 1) ? (SMax-SMin)/(NPoints-1) : 0;
  for (int i=0; i<NPoints; i++) {
    x[i] = SMin+i*step;
  }
  
  TCanvas* c = new TCanvas("c_lhr_008","FC 008",1100,800);
  if (! c->GetShowEventStatus()) c->ToggleEventStatus();
  gPad->SetCrosshair(1);

  TGraph* gr = new TGraph(NPoints,x,prob);

  gr->SetName(Form("gr_bgr_%06i_lhr_008",int(MuB*1000)));
  gr->SetTitle(Form("CLs discovery prob for bgr=%5.3f events",MuB));

  gr->SetMarkerStyle(20);
  gr->Draw("alp");

  gPad->SetGridy(1);
  gPad->SetGridx(1);

  gPad->Modified();
  gPad->Update();
}


void lhr_011(double CL, double MuB=0.1, double N = 5) {
  TLogLHR*   llhr;

  if (llhr == nullptr) llhr = new TLogLHR("lhr_011",CL);
  else                 llhr->SetCL(CL);

  g_llhr = llhr;
  
  double smin = 0;
  double smax = 4*N;

  int NPoints = (smax-smin)/0.01 + 1;
  
  double prob[100000], mean_llhr[100000];

  llhr->MeasInterval(MuB,N,smin,smax,NPoints,prob,mean_llhr);

  TH1D *h_prob, *h_llhr;

  double step;
  if (NPoints > 1) {
    step = (smax-smin)/(NPoints-1);

    h_prob = new TH1D("h_prob","lhr_011 Prob"  ,NPoints,smin-step/2,smax+step/2);
    h_llhr = new TH1D("h_llhr","lhr_011 <LLHR>",NPoints,smin-step/2,smax+step/2);
  }
  else {
    h_prob = new TH1D("h_prob","lhr_011 Prob",  1, smin-0.5, smin+0.5);
    h_llhr = new TH1D("h_llhr","lhr_011 <LLHR>",1, smin-0.5, smin+0.5);
  }
  
  for (int i=0; i<NPoints; i++) {
    // printf("i, prob[i], <llhr[i]> : %3i %12.5e %12.5e\n",i, prob[i],mean_llhr[i]);
    h_prob->SetBinContent(i+1,prob     [i]);
    h_llhr->SetBinContent(i+1,mean_llhr[i]);
  }

  TCanvas* c = new TCanvas("c_lhr_011","FC 011",1400,600);
  gPad->SetLogy(1);

  c->Divide(2,1);
  if (! c->GetShowEventStatus()) c->ToggleEventStatus();

  c->cd(1);
   gPad->SetGridx(1);
   gPad->SetGridy(1);
   gPad->SetCrosshair(1);

  h_prob->GetYaxis()->SetRangeUser(1.e-5,1);
  h_prob->Draw();

  c->cd(2);
  h_llhr->GetYaxis()->SetRangeUser(-15,1);
  h_llhr->Draw();
}
