// -*- mode:c++ -*-

#include "TCanvas.h"
#include "TRandom3.h"
#include "TH1D.h"
#include "TMath.h"

#include "murat/alg/TFeldmanCousinsA.hh"
#include "murat/alg/TLogLHR.hh"

TFeldmanCousinsA* fc;
TLogLHR*          llhr(nullptr);

//-----------------------------------------------------------------------------
// for "5-sigma" consistency, set CL = -1
//-----------------------------------------------------------------------------
void fc_005(double CL, double B=0.1, double  SMin = 0, double SMax = 0, double NPoints = 1) {

  if (!gInterpreter->IsLoaded("lib/libmurat_alg.so")) gInterpreter->Load("lib/libmurat_alg.so");
  
  if (llhr == nullptr) llhr = new TLogLHR("fc_005",CL);
  else                 llhr->SetCL(CL);

  double prob[10000], mean_llhr[10000];

  llhr->ConfInterval(B,SMin,SMax,NPoints,prob,mean_llhr);

  TH1D *h_prob, *h_llhr;
  
  double step;
  if (NPoints > 1) {
    step   = (SMax-SMin)/(NPoints-1);
    h_prob = new TH1D("h_prob","fc_005 Prob"  ,NPoints,SMin-step/2,SMax+step/2);
    h_llhr = new TH1D("h_llhr","fc_005 <LLHR>",NPoints,SMin-step/2,SMax+step/2);
  }
  else {
    h_prob = new TH1D("h_prob","fc_005 Prob",  1, SMin-0.5, SMin+0.5);
    h_llhr = new TH1D("h_llhr","fc_005 <LLHR>",1, SMin-0.5, SMin+0.5);
  }
  
  for (int i=0; i<NPoints; i++) {
    printf("i, prob[i] : %3i %12.5e\n",i, prob[i]);
    h_prob->SetBinContent(i+1,prob[i]);
    h_llhr->SetBinContent(i+1,mean_llhr[i]);
  }

  TCanvas* c = new TCanvas("c_fc_005","FC 005",1400,700);
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
void fc_006(double CL, double MuB=0.1, double N = 5, double  SMin = 0, double SMax = 0, int NPoints = 1) {

  if (!gInterpreter->IsLoaded("lib/libmurat_alg.so")) gInterpreter->Load("lib/libmurat_alg.so");
  
  if (llhr == nullptr) llhr = new TLogLHR("fc_006",CL);
  else                 llhr->SetCL(CL);

  double prob[NPoints], mean_llhr[NPoints];

  llhr->MeasInterval(MuB,N,SMin,SMax,NPoints,prob,mean_llhr);

  TH1D *h_prob, *h_llhr;
  
  double step;
  if (NPoints > 1) {
    step = (SMax-SMin)/(NPoints-1);
    h_prob = new TH1D("h_prob","fc_006 Prob"  ,NPoints,SMin-step/2,SMax+step/2);
    h_llhr = new TH1D("h_llhr","fc_006 <LLHR>",NPoints,SMin-step/2,SMax+step/2);
  }
  else {
    h_prob = new TH1D("h_prob","fc_006 Prob",  1, SMin-0.5, SMin+0.5);
    h_llhr = new TH1D("h_llhr","fc_006 <LLHR>",1, SMin-0.5, SMin+0.5);
  }
  
  for (int i=0; i<NPoints; i++) {
    // printf("i, prob[i], <llhr[i]> : %3i %12.5e %12.5e\n",i, prob[i],mean_llhr[i]);
    h_prob->SetBinContent(i+1,prob     [i]);
    h_llhr->SetBinContent(i+1,mean_llhr[i]);
  }

  TCanvas* c = new TCanvas("c_fc_006","FC 006",1400,600);
  //  gPad->SetLogy(1);
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

//-----------------------------------------------------------------------------
// experiment with the expected background  N
// scan a range of signal hypotheses [SMin,SMax] with a step 'Step'
// determine CLb propability of a CL-level observation of a signal
//-----------------------------------------------------------------------------
void fc_007(double CL, double MuB=0.1, double  SMin = 0, double SMax = 0, int NPoints = 1) {

  if (!gInterpreter->IsLoaded("lib/libmurat_alg.so")) gInterpreter->Load("lib/libmurat_alg.so");
  
  if (llhr == nullptr) llhr = new TLogLHR("fc_007",CL);
  else                 llhr->SetCL(CL);

  double  x[NPoints], prob[NPoints];
  
  llhr->DiscoveryProbCLb(MuB,SMin,SMax,NPoints,x,prob);

  TCanvas* c = new TCanvas("c_fc_007","FC 006",1100,800);
  if (! c->GetShowEventStatus()) c->ToggleEventStatus();
  gPad->SetCrosshair(1);

  TGraph* gr = new TGraph(NPoints,x,prob);

  gr->SetName(Form("gr_bgr_%06i_fc_007",int(MuB*1000)));
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
void fc_008(double CL, double MuB=0.1, double  SMin = 0, double SMax = 0, int NPoints = 1) {

  if (!gInterpreter->IsLoaded("lib/libmurat_alg.so")) gInterpreter->Load("lib/libmurat_alg.so");
  
  if (llhr == nullptr) llhr = new TLogLHR("fc_007",CL);
  else                 llhr->SetCL(CL);

  double  x[NPoints], prob[NPoints];
  
  llhr->DiscoveryProbCLs(MuB,SMin,SMax,NPoints,x,prob);

  double step = (NPoints > 1) ? (SMax-SMin)/(NPoints-1) : 0;
  for (int i=0; i<NPoints; i++) {
    x[i] = SMin+i*step;
  }
  
  TCanvas* c = new TCanvas("c_fc_008","FC 008",1100,800);
  if (! c->GetShowEventStatus()) c->ToggleEventStatus();
  gPad->SetCrosshair(1);

  TGraph* gr = new TGraph(NPoints,x,prob);

  gr->SetName(Form("gr_bgr_%06i_fc_008",int(MuB*1000)));
  gr->SetTitle(Form("CLs discovery prob for bgr=%5.3f events",MuB));

  gr->SetMarkerStyle(20);
  gr->Draw("alp");

  gPad->SetGridy(1);
  gPad->SetGridx(1);

  gPad->Modified();
  gPad->Update();
}


