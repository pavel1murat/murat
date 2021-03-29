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
// plot B and (B+S) distributions of the number of expected events
//-----------------------------------------------------------------------------
void fc_001() {
  
  fc = new TFeldmanCousinsA("a",-1,1);

  double b = 0.5;
  double s = 3.0;
  
  fc->Init(b,s);

  TCanvas* c = new TCanvas("c","c",800,600);
  
  fc->fHist.fBgProb->Draw();
  fc->fHist.fBgProb->SetLineColor(kRed+2);
  
  fc->fHist.fBsProb->Draw("sames");

  c->Modified();
  c->Update();
}

void fc_003(double b=0.1, double s = 5.0 ) {
  //  double b(0.4), s);

  if ( !gInterpreter->IsLoaded("lib/libmurat_alg.so")) {
    gInterpreter->Load("lib/libmurat_alg.so");
  }
  
  fc = new TFeldmanCousinsA("a",0.9);
  
  fc->Init(b,s);

  TRandom3 r3;
  int nbx(1000);
  
  TH1D* h20 = new TH1D("h20_fc_003","X(B+S)"          ,50,0,50);
  TH1D* h31 = new TH1D("h31_fc_003","LH(B|B+S)"       ,nbx,0, 1);
  TH1D* h32 = new TH1D("h32_fc_003","LLH(B|B)"        ,nbx,0, 1);
  TH1D* h33 = new TH1D("h33_fc_003","LLH(B-(B+S)|B+S)",nbx,-100, 100);
  TH1D* h40 = new TH1D("h40_fc_003","prob"            ,nbx,-100, 100);
  
  long int nexp = 10000000;
  for (long int i=0; i<nexp; i++) {

    int xbs = r3.Poisson(b+s);          // B+S
    h20->Fill(xbs);
					// H0 likelihood
    double lh_b = fc->fBgProb[xbs];
    h31->Fill(lh_b);
					// H1 likelihood
    double lh_bs = fc->fBsProb[xbs];  
    h32->Fill(lh_bs);

    double llh_bs = TMath::Log(lh_b/lh_bs);
    h33->Fill(llh_bs);
  }

  for (int i=1; i<=nbx; i++) {
    double prob = 1-h33->Integral(i+1,nbx)/((double) nexp);
    // printf("i = %3i prob=%12.5e\n",i,prob);
    h40->SetBinContent(i,prob);
  }
  
  TCanvas* c = new TCanvas("c_fc_003","FC 003",1200,900);

  c->Divide(2,2);
  c->cd(1);
  h20->Draw();
  
  c->cd(2);
  h31->SetLineColor(kRed+2);
  h31->SetFillColor(kRed+2);
  h31->SetFillStyle(3001);
  h31->Draw("");
  h32->Draw("sames");

  c->cd(3);
  h33->Draw("");

  c->cd(4);
  h40->Draw("");

  c->Modified();
  c->Update();
}

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
void fc_005(double B=0.1, int N = 5, double  SMin = 0, double SMax = 0, double NSteps = 1) {

  if (!gInterpreter->IsLoaded("lib/libmurat_alg.so")) gInterpreter->Load("lib/libmurat_alg.so");
  
  llhr = new TLogLHR("fc_005",-1);

  double prob[1000];

  TH1D* hist;

  if (NSteps > 1) hist = new TH1D("hist","hist",NSteps,SMin,SMax);
  else            hist = new TH1D("hist","hist",    1, SMin-0.5, SMin+0.5);
  
  llhr->ConfInterval(B,N,SMin,SMax,NSteps-1,prob);

  for (int i=0; i<NSteps; i++) {
    hist->SetBinContent(i+1,prob[i]);
  }

  TCanvas* c = new TCanvas("c_fc_005","FC 005",1000,600);
  gPad->SetLogy(1);
  hist->Draw();
}
