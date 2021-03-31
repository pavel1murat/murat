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

  double prob[1000];

  TH1D* hist;

  if (NSteps > 1) hist = new TH1D("hist","hist",NSteps,SMin,SMax);
  else            hist = new TH1D("hist","hist",    1, SMin-0.5, SMin+0.5);
  
  llhr->ConfInterval(B,N,SMin,SMax,NSteps-1,prob);

  for (int i=0; i<NSteps; i++) {
    hist->SetBinContent(i+1,prob[i]);
  }

  TCanvas* c = new TCanvas("c_fc_004","FC 004",1200,900);
  hist->Draw();
}


//-----------------------------------------------------------------------------
// try upper limit - set nsteps = 1, SMin = 0;
//-----------------------------------------------------------------------------
void fc_005(double CL, double B=0.1, int N = 5, double  SMin = 0, double SMax = 0, double NSteps = 1) {

  if (!gInterpreter->IsLoaded("lib/libmurat_alg.so")) gInterpreter->Load("lib/libmurat_alg.so");
  
  llhr = new TLogLHR("fc_005",CL);

  double prob[1000];

  llhr->ConfInterval(B,N,SMin,SMax,NSteps,prob);

  TH1D* hist;
  
  double step;
  if (NSteps > 1) {
    step = (SMax-SMin)/(NSteps-1);
    hist = new TH1D("hist","hist",NSteps,SMin-step/2,SMax+step/2);
  }
  else {
    hist = new TH1D("hist","hist",    1, SMin-0.5, SMin+0.5);
  }
  
  for (int i=0; i<NSteps; i++) {
    printf("i, prob[i] : %3i %12.5e\n",i, prob[i]);
    hist->SetBinContent(i+1,prob[i]);
  }

  TCanvas* c = new TCanvas("c_fc_005","FC 005",1000,600);
  gPad->SetLogy(1);
  gPad->SetGridy(1);
  gPad->SetCrosshair(1);
  hist->GetYaxis()->SetRangeUser(1.e-2,1);
  hist->Draw();
}


//-----------------------------------------------------------------------------
// try measurement
//-----------------------------------------------------------------------------
void fc_006(double CL, double B=0.1, int N = 5, double  SMin = 0, double SMax = 0, double NSteps = 1) {

  if (!gInterpreter->IsLoaded("lib/libmurat_alg.so")) gInterpreter->Load("lib/libmurat_alg.so");
  
  llhr = new TLogLHR("fc_006",CL);

  double prob[1000];

  llhr->MeasInterval(B,N,SMin,SMax,NSteps,prob);

  TH1D* hist;
  
  double step;
  if (NSteps > 1) {
    step = (SMax-SMin)/(NSteps-1);
    hist = new TH1D("h_fc_006","hist",NSteps,SMin-step/2,SMax+step/2);
  }
  else {
    hist = new TH1D("h_fc_006","hist",    1, SMin-0.5, SMin+0.5);
  }
  
  for (int i=0; i<NSteps; i++) {
    printf("i, prob[i] : %3i %12.5e\n",i, prob[i]);
    hist->SetBinContent(i+1,prob[i]);
  }

  TCanvas* c = new TCanvas("c_fc_006","FC 006",1000,600);
  //  gPad->SetLogy(1);
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  gPad->SetCrosshair(1);
  hist->GetYaxis()->SetRangeUser(1.e-2,1);
  hist->Draw();
}


