///////////////////////////////////////////////////////////////////////////////
// 2016-03-18 P.Murat: default crystal ball parameterization 
//
// - gaussian   for dx > -|abs_alpha|
// - power tail for dx < -|abs_alpha|
//
//   no positive tails
///////////////////////////////////////////////////////////////////////////////
/* --

Double_t CrystalBall(Double_t *x,Double_t *par) {

  Double_t t = (x[0]-par[2])/par[3];
  if (par[0] < 0) t = -t;

  Double_t absAlpha = fabs((Double_t)par[0]);

  if (t >= -absAlpha) {
    return par[4]*exp(-0.5*t*t);
  }
  else {
    Double_t a =  TMath::Power(par[1]/absAlpha,par[1])*exp(-0.5*absAlpha*absAlpha);
    Double_t b= par[1]/absAlpha - absAlpha; 

    return par[4]*(a/TMath::Power(b - t, par[1]));
  }
}

-- */
#include "TH1F.h"
#include "TF1.h"
#include "TMath.h"
#include "Stntuple/val/stntuple_val_functions.hh"

//-----------------------------------------------------------------------------
// A1*exp((x-x0)^2/2sig^2 for (x-x0)/sig > -alp
// A1*A2*(B-(x-x0)/sig)^{-n}
// A2 = 
// P[0] = N
// P[1] = x[0]
// P[2] = sigma
// P[3] = alpha  (> 0)
// P[4] = n
//-----------------------------------------------------------------------------
double f_crystal_ball(double* X, double* P) {
  double f, alpha, abs_alpha, n, a, b;

  double dx = (X[0]-P[1])/P[2];

  alpha     = P[3];
  abs_alpha = fabs(alpha);
  n         = P[4];

  if (dx > -abs_alpha) {
    f = P[0]*TMath::Exp(-dx*dx/2.);
  }
  else {
    a = TMath::Power((n/abs_alpha),n)*TMath::Exp(-(alpha*alpha)/2);
    b = n/abs_alpha-abs_alpha;
    f = P[0]*a*TMath::Power((b-dx),-n);
  }

  return f;
}

namespace {
  TH1F* h; 
  TF1*  f; 
};


//-----------------------------------------------------------------------------
void cb_init_parameters(TF1* F, double P0, double P1, double P2, double P3, double P4) 
{
  F->SetParameter(0,P0);
  F->SetParameter(1,P1);
  F->SetParameter(2,P2);
  F->SetParameter(3,P3);
  F->SetParameter(4,P4);
}

//-----------------------------------------------------------------------------
void cb_fit(double XMin, double XMax, const char* FitOpt = "") {
  h->GetXaxis()->SetRangeUser(XMin,XMax);
  h->Fit(f,FitOpt,"",XMin,XMax);
}

//-----------------------------------------------------------------------------
void fit_crystal_ball(const char* File, const char* Module, const char* Hist, 
		      double X0, double XMin, double XMax) {

  h = (TH1F*) gh1(File,Module,Hist)->Clone("h_fit");

  f = new TF1("f_crystal_ball",f_crystal_ball,XMin,XMax,5);

  f->SetParName(0,"anorm");
  f->SetParName(1,"mean");
  f->SetParName(2,"sigma");
  f->SetParName(3,"alpha");
  f->SetParName(4,"N");

  f->SetParLimits(2,0.,5.);

  double anorm = h->GetEntries()/10;
  cb_init_parameters(f,anorm,X0,0.02,0.05,3.);
  cb_fit(XMin,XMax);

  // h->Draw();
  // f->Draw("same");
}

//-----------------------------------------------------------------------------
// fit a histogram, normalized to a unit area
//-----------------------------------------------------------------------------
void fit_crystal_ball_norm(const char* File, const char* Module, const char* Hist, 
			   double X0, double XMin, double XMax) {

  h = (TH1F*) gh1(File,Module,Hist)->Clone("h_fit");
  h->Scale(1./h->Integral());

  f = new TF1("f_crystal_ball",f_crystal_ball,XMin,XMax,7);

  f->SetParName(0,"anorm");
  f->SetParName(1,"mean");
  f->SetParName(2,"sigma");
  f->SetParName(3,"alpha");
  f->SetParName(4,"N");

  f->SetParLimits(2,0.,10.);

  double anorm = h->Integral()/10;
  cb_init_parameters(f,anorm,X0,2.,3.,1.);
  cb_fit(XMin,XMax);
}


//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
void cb_fit_ep(TH1* Hist, double X0, double XMin, double XMax, const char* FitOpt) {

  h      = (TH1F*) Hist;
  double anorm = h->GetEntries()/15;


  f = new TF1("f_crystal_ball",f_crystal_ball,XMin,XMax,5);

  f->SetParName(0,"anorm");
  f->SetParName(1,"mean");
  f->SetParName(2,"sigma");
  f->SetParName(3,"alpha");
  f->SetParName(4,"N");


  cb_init_parameters(f,anorm,X0,0.02,0.2,10.);
  cb_fit            (XMin,XMax,FitOpt);
}

