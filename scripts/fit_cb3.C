///////////////////////////////////////////////////////////////////////////////
// cb3: pseudo Crystall Ball function : aymmetric gaussian with power-law tails on both sides
// make it a C++ class to simplify fitting multiple functions
//
// usage:
// ------
// .L murat/scripts/fit_cb2.C
// cbx = new cb2("x");
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
class cb3: public TNamed { 
public:

  TH1F* fHist;    // not owned by the fitter
  TH1F* fHGaus;   // owned by the fitter

  TF1*  fFunc;
  TF1*  fGaus;

  cb3(const char* Name = "x");
  ~cb3();

  static double  agauss(double* X, double* P);
  static double  f_cb3 (double* X, double* P);    // asymmetric double Crystall Ball itself

  void           init_parameters    (double P0, double P1, double P2, double P3, 
				     double P4, double P5, double P6, double P7);

  void           fit_histogram      (double XMin, double XMax);
//-----------------------------------------------------------------------------
// two functions doing the same, but with different interfaces
//-----------------------------------------------------------------------------
  void           fit(TH1* Hist, double X0, double XMin, double XMax, 
		     double SigmaL = -1., double SigmaR=-1.);

  void           fit(const char* File, const char* Module, const char* Hist, 
		     double X0, double XMin, double XMax, 
		     double SigmaL=-1., double SigmaR=-1.);

  TF1* GetFunc() { return fFunc; }
  TH1* GetHist() { return (TH1*) fHist; }

};


//-----------------------------------------------------------------------------
cb3::cb3(const char* Name): TNamed(Name,Name) {
  fHist  = nullptr;
  fHGaus = nullptr;
  
  fFunc = new TF1(Form("%s_f_cb3"  ,GetName()),cb3::f_cb3 ,-1,1,8);
  fFunc->SetParName(0,"anorm");
  fFunc->SetParName(1,"mean");
  fFunc->SetParName(2,"Ls");
  fFunc->SetParName(3,"Rs");
  fFunc->SetParName(4,"alpha1");
  fFunc->SetParName(5,"N1");
  fFunc->SetParName(6,"alpha2");
  fFunc->SetParName(7,"N2");

  fFunc->SetParLimits(2,0.,1.);
  fFunc->SetParLimits(3,0.,1.);
  //  F->SetParLimits(5,0.,5.);
  //  F->SetParLimits(6,0.,5.);

  fFunc->SetLineWidth(1);

  fGaus = new TF1(Form("%s_f_agaus",GetName()),cb3::agauss,-1,1,4);
}

//-----------------------------------------------------------------------------
cb3::~cb3() {
  delete fGaus;
}

//-----------------------------------------------------------------------------
double cb3::agauss(double* X, double* P) {
  double x0, f(0), sig, xdx;

  x0         = P[1];
  double dx  = X[0]-x0;

  if (dx < 0) sig = P[2];
  else        sig = P[3];

  xdx        = (X[0]-x0)/sig;

  f = P[0]*TMath::Exp(-xdx*xdx/2.);
  
  return f;
}

//-----------------------------------------------------------------------------
double cb3::f_cb3(double* X, double* P) {
  double x0, x1, f, sigl, sigr, alpha1, abs_alpha1, alpha2, abs_alpha2, n1, n2, a, b, dx0, dx1, xdx;

  x0         = P[1];
  sigl       = P[2];
  sigr       = P[3];
  alpha1     = P[4];
  n1         = P[5];
  alpha2     = P[6];
  n2         = P[7];

  dx0        = X[0]-x0;
  if (dx0 < 0) xdx = dx0/sigl;
  else         xdx = dx0/sigr;

  abs_alpha1 = fabs(alpha1);
  abs_alpha2 = fabs(alpha2);

  if (xdx < -abs_alpha1) {
    a = TMath::Power((n1/abs_alpha1),n1)*TMath::Exp(-(alpha1*alpha1)/2);
    b = n1/abs_alpha1-abs_alpha1;
    f = P[0]*a*TMath::Power((b-xdx),-n1);
  }
  else if (xdx < abs_alpha2) { 
					// gauss
    f = P[0]*TMath::Exp(-xdx*xdx/2.);
  }
  else {
    a = TMath::Power((n2/abs_alpha2),n2)*TMath::Exp(-(alpha2*alpha2)/2);
    b = -n2/abs_alpha2+abs_alpha2;
    f = P[0]*a*TMath::Power((xdx-b),-n2);
  }
  
  return f;
}

//-----------------------------------------------------------------------------
void cb3::fit_histogram(double XMin, double XMax) {
  fHist->GetXaxis()->SetRangeUser(XMin,XMax);
  fHist->Fit(fFunc,"L","",XMin,XMax);
}

//-----------------------------------------------------------------------------
void cb3::init_parameters(double P0, double P1, double P2, double P3, 
			  double P4, double P5, double P6, double P7) {
  fFunc->SetParameter(0,P0);
  fFunc->SetParameter(1,P1);
  fFunc->SetParameter(2,P2);
  fFunc->SetParameter(3,P3);
  fFunc->SetParameter(4,P4);
  fFunc->SetParameter(5,P5);
  fFunc->SetParameter(6,P6);
  fFunc->SetParameter(7,P7);
}

//-----------------------------------------------------------------------------
void cb3::fit(TH1* Hist, double X0, double XMin, double XMax, double SigmaL = -1, double SigmaR = -1.) {

  fHist      = (TH1F*) Hist->Clone("h_fit_cb3");
  double anorm = fHist->GetEntries()/10;

  printf("anorm = %12.5e\n",anorm);

  init_parameters(anorm,X0, 0.180, 0.180, 1,.4,1.,4);

  if (SigmaL > 0) {			// effectively, fix Sigma
    fFunc->SetParameter(2,SigmaL);
    fFunc->SetParLimits(2,SigmaL*(1+1.e-5),SigmaL*(1+1.e-5));
  }

  if (SigmaR > 0) {			// effectively, fix Sigma
    fFunc->SetParameter(3,SigmaR);
    fFunc->SetParLimits(3,SigmaR*(1+1.e-5),SigmaR*(1+1.e-5));
  }

  fit_histogram(XMin,XMax);
					// gaussian component of the fit, to display
  fGaus->SetMinimum(XMin);
  fGaus->SetMaximum(XMax);
  fGaus->SetParameter(0,fFunc->GetParameter(0));
  fGaus->SetParameter(1,fFunc->GetParameter(1));
  fGaus->SetParameter(2,fFunc->GetParameter(2));
  fGaus->SetParameter(3,fFunc->GetParameter(3));

  fHGaus = (TH1F*) fHist->Clone(Form("%s_h_agaus",GetName()));

  fHGaus->Reset();

  for (int i=1; i<fHist->GetNbinsX(); i++) {
    float x = fHist->GetBinCenter(i);
    float y = fGaus->Eval(x);
    fHGaus->SetBinContent(i,y);
    fHGaus->SetBinError  (i,0);
  }
  
  fHGaus->SetLineColor(kBlue+3);
  fHGaus->SetFillColor(kBlue+3);
  fHGaus->SetFillStyle(3003);
  
  fHGaus->Draw("same");
//-----------------------------------------------------------------------------
// calculate the background fraction: high-momentum side , above the gaussian
// use region from alpha2 to XMax
//-----------------------------------------------------------------------------
  double total  = fHist->Integral(1,Hist->GetNbinsX());

  double x0     = fFunc->GetParameter("mean"  );
  double sigmaL = fFunc->GetParameter("Ls" );
  double sigmaR = fFunc->GetParameter("Rs" );
  double alpha2 = fFunc->GetParameter("alpha2");
  double x1     = x0+fabs(sigmaR*alpha2);

  printf("x0 = %10.4f, sigmaL = %10.4f sigmaR = %10.4f alpha2 = %10.4f\n",x0,sigmaL,sigmaR,alpha2);
 
  double htail = (fFunc->Integral(x1,XMax)-fGaus->Integral(x1,XMax))/fHist->GetBinWidth(1);

  printf("total = %12.5f, HTail = %12.5f, htail/total : %12.5e\n",total,htail,htail/total);
}

//-----------------------------------------------------------------------------
void cb3::fit(const char* File, const char* Module, const char* Hist, 
	      double X0, double XMin, double XMax, double SigmaL = -1., double SigmaR = -1.) {

  TH1* h = (TH1F*) gh1(File,Module,Hist);

  this->fit(h,X0,XMin,XMax,SigmaL,SigmaR);
}

