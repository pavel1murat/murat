// -*- mode:c++ -*-

#include "TCanvas.h"
#include "TRandom3.h"
#include "TH1D.h"
#include "TMath.h"

#include "murat/alg/TFeldmanCousinsA.hh"
#include "murat/alg/TLogLHR.hh"

TFeldmanCousinsA* fc;
TLogLHR*          llhr;

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
void fc_004(double B=0.1, int N = 5.0, double  SMin = 0, double SMax = 20, double NSteps = 100) {

  if (!gInterpreter->IsLoaded("lib/libmurat_alg.so")) gInterpreter->Load("lib/libmurat_alg.so");
  
  llhr = new TLogLHR("a",0.9);

  double prob[1000], mean_llhr[1000];

  TH1D *h_prob, *h_llhr;

  if (NSteps > 1) {
    h_prob = new TH1D("h_prob","fc_004 Prob"  ,NSteps,SMin,SMax);
    h_llhr = new TH1D("h_llhr","fc_004 <LLHR>",NSteps,SMin,SMax);
  }
  else {
    h_prob = new TH1D("h_prob","fc_004 Prob",  1, SMin-0.5, SMin+0.5);
    h_llhr = new TH1D("h_llhr","fc_004 <LLHR>",1, SMin-0.5, SMin+0.5);
  }
  
  llhr->ConfInterval(B,N,SMin,SMax,NSteps-1,prob,mean_llhr);

  for (int i=0; i<NSteps; i++) {
    h_prob->SetBinContent(i+1,prob[i]);
    h_llhr->SetBinContent(i+1,mean_llhr[i]);
  }

  TCanvas* c = new TCanvas("c_fc_004","FC 004",1600,900);
  c->Divide(2,1);

  c->cd(1); 
  h_prob->Draw();

  c->cd(2);
  h_llhr->Draw();
}


//-----------------------------------------------------------------------------
// try upper limit - set nsteps = 1, SMin = 0;
//-----------------------------------------------------------------------------
void fc_005(double CL, double B=0.1, int N = 5, double  SMin = 0, double SMax = 0, double NSteps = 1) {

  if (!gInterpreter->IsLoaded("lib/libmurat_alg.so")) gInterpreter->Load("lib/libmurat_alg.so");
  
  llhr = new TLogLHR("fc_005",CL);

  double prob[1000], mean_llhr[1000];

  llhr->ConfInterval(B,N,SMin,SMax,NSteps,prob,mean_llhr);

  TH1D *h_prob, *h_llhr;
  
  double step;
  if (NSteps > 1) {
    step   = (SMax-SMin)/(NSteps-1);
    h_prob = new TH1D("h_prob","fc_005 Prob"  ,NSteps,SMin-step/2,SMax+step/2);
    h_llhr = new TH1D("h_llhr","fc_005 <LLHR>",NSteps,SMin-step/2,SMax+step/2);
  }
  else {
    h_prob = new TH1D("h_prob","fc_005 Prob",  1, SMin-0.5, SMin+0.5);
    h_llhr = new TH1D("h_llhr","fc_005 <LLHR>",1, SMin-0.5, SMin+0.5);
  }
  
  for (int i=0; i<NSteps; i++) {
    printf("i, prob[i] : %3i %12.5e\n",i, prob[i]);
    h_prob->SetBinContent(i+1,prob[i]);
    h_llhr->SetBinContent(i+1,mean_llhr[i]);
  }

  TCanvas* c = new TCanvas("c_fc_005","FC 005",1400,700);
  c->Divide(2,1);

  c->cd(1);
  gPad->SetLogy(1);
  gPad->SetGridy(1);
  gPad->SetCrosshair(1);

  h_prob->GetYaxis()->SetRangeUser(1.e-2,1);
  h_prob->Draw();
  
  c->cd(2);
  h_llhr->Draw();
}


//-----------------------------------------------------------------------------
// try measurement
//-----------------------------------------------------------------------------
void fc_006(double CL, double B=0.1, int N = 5, double  SMin = 0, double SMax = 0, double NSteps = 1) {

  if (!gInterpreter->IsLoaded("lib/libmurat_alg.so")) gInterpreter->Load("lib/libmurat_alg.so");
  
  llhr = new TLogLHR("fc_006",CL);

  double prob[1000], mean_llhr[1000];

  llhr->MeasInterval(B,N,SMin,SMax,NSteps,prob,mean_llhr);

  TH1D *h_prob, *h_llhr;
  
  double step;
  if (NSteps > 1) {
    step = (SMax-SMin)/(NSteps-1);
    h_prob = new TH1D("h_prob","fc_006 Prob"  ,NSteps,SMin-step/2,SMax+step/2);
    h_llhr = new TH1D("h_llhr","fc_006 <LLHR>",NSteps,SMin-step/2,SMax+step/2);
  }
  else {
    h_prob = new TH1D("h_prob","fc_006 Prob",  1, SMin-0.5, SMin+0.5);
    h_llhr = new TH1D("h_llhr","fc_006 <LLHR>",1, SMin-0.5, SMin+0.5);
  }
  
  for (int i=0; i<NSteps; i++) {
    printf("i, prob[i], <llhr[i]> : %3i %12.5e %12.5e\n",i, prob[i],mean_llhr[i]);
    h_prob->SetBinContent(i+1,prob[i]);
    h_llhr->SetBinContent(i+1,mean_llhr[i]);
  }

  TCanvas* c = new TCanvas("c_fc_006","FC 006",1400,600);
  //  gPad->SetLogy(1);
  c->Divide(2,1);
  c->cd(1);
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  gPad->SetCrosshair(1);
  h_prob->GetYaxis()->SetRangeUser(1.e-5,1);
  h_prob->Draw();

  c->cd(2);
  h_llhr->Draw();
}


