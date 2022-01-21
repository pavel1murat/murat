///////////////////////////////////////////////////////////////////////////////
// cb4: lognormal
// make it a C++ class to simplify fitting multiple functions
//
// usage:
// ------
// .L murat/scripts/fit_lognormal.C
// cbx = new fit_lognormal("x");
// cbx->fit(file,....)
//
// or
//
// TH1F*hist = ...
// cbx->fit(...)
//
// to find the momentum value corresponding to the half-maximum,
// scan a range of momenta with 
// 
// cbx->GetFunc()->Eval(104.5) 
///////////////////////////////////////////////////////////////////////////////
#include "TMath.h"
#include "TNamed.h"
#include "TF1.h"

//-----------------------------------------------------------------------------
class fit_lognormal: public TNamed { 
public:

  TH1F* fHist;    // not owned by the fitter
  TH1F* fHGaus;   // owned by the fitter

  TF1*  fFunc;
  TF1*  fGaus;

  fit_lognormal(const char* Name = "x");
  ~fit_lognormal();

  static double  core(double* X, double* P);      // lognormal core

  void           init_parameters    (double P0, double P1, double P2, double P3);

  void           fit_histogram      (double XMin, double XMax);
//-----------------------------------------------------------------------------
// two functions doing the same, but with different interfaces
//-----------------------------------------------------------------------------
  void           fit(TH1* Hist, double X0, double XMin, double XMax);

  void           fit(const char* File, const char* Module, const char* Hist, 
		     double X0, double XMin, double XMax);

  TF1* GetFunc() { return fFunc; }
  TH1* GetHist() { return (TH1*) fHist; }

};


//-----------------------------------------------------------------------------
fit_lognormal::fit_lognormal(const char* Name): TNamed(Name,Name) {
  fHist  = nullptr;
  fHGaus = nullptr;
  
  fFunc = new TF1(Form("%s_f_fit_lognormal"  ,GetName()),fit_lognormal::f_core ,-1,1,6);
  fFunc->SetParName(0,"anorm");
  fFunc->SetParName(1,"x0");
  fFunc->SetParName(2,"scale");
  fFunc->SetParName(3,"sigm");

  // fFunc->SetParLimits(2,0.,1.);
  // fFunc->SetParLimits(3,0.,1.);
  //  F->SetParLimits(5,0.,5.);
  //  F->SetParLimits(6,0.,5.);

  fFunc->SetLineWidth(1);

  //  fGaus = new TF1(Form("%s_f_agaus",GetName()),cb4::agauss,-1,1,4);
}

//-----------------------------------------------------------------------------
fit_lognormal::~fit_lognormal() {
  // delete fGaus;
}

//-----------------------------------------------------------------------------
// dx = (x-P[3])/P[2]
// p[0]*exp(p[1] * (dx − exp(dx)) )
// default Landau distribution for energy losses: P[2] = 2
// decant fit: Landau (-x) constant=260000 mpv=0.5 sigma=0.18
//-----------------------------------------------------------------------------
double fit_lognormal::core(double* X, double* P) {

  double x0    = P[1];
  double s1    = P[2]; 
  double dx    = (X[0]-x0);
  double logdx = log(dx/P[2]);
  double sig   = P[3];
  double sig2  = sig*sig;

  double f     = P[0]*exp(-(logdx*logdx)/(2*sig*sig))/dx;

  return f;
}

//-----------------------------------------------------------------------------
double fit_lognormal::f_cb4(double* X, double* P) {
  double x0, abs_dx0, x1, f, sigm, beta1, abs_beta1, beta2, abs_beta2, n1, n2, a, b, dx0, dx1, xdx;

  x0         = P[1];
  sigm       = P[2];
  a          = P[3];
  beta1      = P[4];   // x1/sigm
  beta2      = P[5];   // x2/sigm

  dx0        = (X[0]-x0)/sigm;
  abs_dx0    = fabs(dx0);
					// in units of sigma
  abs_beta1  = fabs(beta1);
  abs_beta2  = fabs(beta2);

  if (xdx < -abs_beta1) {
    n1 = abs_beta1*a*(1+exp(-abs_beta1));
    b  = exp((a*(abs_dx0-exp(-abs_dx0))))/pow(beta1,n1);
    f  = P[0]*b*pow(dx0,n1);
  }
  else if (xdx < abs_beta2) { 
					// gauss
    f = P[0]*exp(a*(dx0-exp(-dx0)));
  }
  else {
    n2 = beta2*a*(1+exp(-beta2));
    b  = exp((a*(dx0-exp(-dx0))))/pow(beta2,n2);
    f  = P[0]*b*pow(dx0,n2);
  }
  
  return f;
}

//-----------------------------------------------------------------------------
void fit_lognormal::fit_histogram(double XMin, double XMax) {
  fHist->GetXaxis()->SetRangeUser(XMin,XMax);
  fHist->Fit(fFunc,"L","",XMin,XMax);
}

//-----------------------------------------------------------------------------
void fit_lognormal::init_parameters(double P0, double P1, double P2, double P3, 
			  double P4, double P5) {
  fFunc->SetParameter(0,P0);
  fFunc->SetParameter(1,P1);
  fFunc->SetParameter(2,P2);
  fFunc->SetParameter(3,P3);
  fFunc->SetParameter(4,P4);
  fFunc->SetParameter(5,P5);
}

//-----------------------------------------------------------------------------
void fit_lognormal::fit(TH1* Hist, double X0, double XMin, double XMax, double Sigma = -1) {

  fHist      = (TH1F*) Hist->Clone("h_fit_cb4");
  double anorm = fHist->GetEntries()/50;

  printf("anorm = %12.5e\n",anorm);

  init_parameters(anorm,X0, 0.180,1., 1.,1.);

  fit_histogram(XMin,XMax);
  // 					// gaussian component of the fit, to display
  // fGaus->SetMinimum(XMin);
  // fGaus->SetMaximum(XMax);
  // fGaus->SetParameter(0,fFunc->GetParameter(0));
  // fGaus->SetParameter(1,fFunc->GetParameter(1));
  // fGaus->SetParameter(2,fFunc->GetParameter(2));
  // fGaus->SetParameter(3,fFunc->GetParameter(3));

  // fHGaus = (TH1F*) fHist->Clone(Form("%s_h_agaus",GetName()));

  // fHGaus->Reset();

  // for (int i=1; i<fHist->GetNbinsX(); i++) {
  //   float x = fHist->GetBinCenter(i);
  //   float y = fGaus->Eval(x);
  //   fHGaus->SetBinContent(i,y);
  //   fHGaus->SetBinError  (i,0);
  // }
  
  // fHGaus->SetLineColor(kBlue+3);
  // fHGaus->SetFillColor(kBlue+3);
  // fHGaus->SetFillStyle(3003);
  
  // fHGaus->Draw("same");
//-----------------------------------------------------------------------------
// calculate the background fraction: high-momentum side , above the gaussian
// use region from alpha2 to XMax
//-----------------------------------------------------------------------------
  // double total  = fHist->Integral(1,Hist->GetNbinsX());

  // double x0     = fFunc->GetParameter("mean"  );
  // double sigmaL = fFunc->GetParameter("Ls" );
  // double sigmaR = fFunc->GetParameter("Rs" );
  // double alpha2 = fFunc->GetParameter("alpha2");
  // double x1     = x0+fabs(sigmaR*alpha2);

  // printf("x0 = %10.4f, sigmaL = %10.4f sigmaR = %10.4f alpha2 = %10.4f\n",x0,sigmaL,sigmaR,alpha2);
 
  // double htail = (fFunc->Integral(x1,XMax)-fGaus->Integral(x1,XMax))/fHist->GetBinWidth(1);

  // printf("total = %12.5f, HTail = %12.5f, htail/total : %12.5e\n",total,htail,htail/total);
}

//-----------------------------------------------------------------------------
void fit_lognormal::fit(const char* File, const char* Module, const char* Hist, 
	      double X0, double XMin, double XMax, double Sigma = -1) {

  TH1* h = (TH1F*) gh1(File,Module,Hist);

  this->fit(h,X0,XMin,XMax,Sigma);
}

