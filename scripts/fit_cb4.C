///////////////////////////////////////////////////////////////////////////////
// parameterization by Shihua Huang: Landau core with the power law tails on both sides
// make it a C++ class to simplify fitting multiple functions
//
// usage:
// ------
// .L murat/scripts/fit_cb4.C
// cbx = new cb4("x");
// cbx->init_parameters(...)
//
// and then either
//
// cbx->fit(file,module,....)
//
// or
//
// TH1F* hist = ... 
// cbx->fHist = (TH1F*) hist->Clone("cbx_fhist")
// cbx->fit_histogram(hist,xmin,xmax)
//
// to find the momentum value corresponding to the half-maximum,
// scan a range of momenta with 
// 
// cbx->GetFunc()->Eval(104.5) 
// Minimizer is Minuit2 / Migrad
// Chi2                      =      121.367
// NDf                       =           67
// Edm                       =  9.73605e-05
// NCalls                    =           85
// anorm                     =      290.015   +/-   0           
// x0                        =      104.398   +/-   0           
// a                         =     0.666853   +/-   0           
// b                         =      3.08502   +/-   0           
// alp1                      =      1.27325   +/-   0           
// n1                        =       7.7983   +/-   0           
// alp2                      =        7.132   +/-   0           
// n2                        =            4   +/-   0           
//
// bpip4b0:  x3->init_parameters(860,69.05,0.45,3.3,2.16,6.8,0.9,2.1)
// cpos1b0:  x2->init_parameters(200,92.3,0.67,3.1,1.3,8,7,4)
// rpc04b0:  x3->init_parameters(1280,128.64,0.158,11.2,5,2,1.3,2.)

///////////////////////////////////////////////////////////////////////////////
#include "TMath.h"
#include "TNamed.h"
#include "TF1.h"

#include "Stntuple/val/stntuple_val_functions.hh"

//-----------------------------------------------------------------------------
class cb4: public TNamed { 
public:

  static int fgDebug;

  TH1F* fHist;    // not owned by the fitter
  TH1F* fHCore;   // owned by the fitter

  TF1*  fFunc;
  TF1*  fCore;

  TString  fFitOpt;

  cb4(const char* Name = "x");
  ~cb4();

  static double  core(double* X, double* P);      // Landau core

  static double  f_cb4 (double* X, double* P);    // Landau core + power tails

  void           init_parameters    (double P0, double P1, double P2, double P3, 
				     double P4, double P5, double P6, double P7);

  void           fit_histogram      (TH1F* Hist, double XMin, double XMax);
//-----------------------------------------------------------------------------
// two functions doing the same, but with different interfaces
//-----------------------------------------------------------------------------
  void           fit(double XMin, double XMax);

  void           fit(const char* File, const char* Module, const char* Hist, 
		     double XMin, double XMax);

  TF1* GetFunc() { return fFunc; }
  TH1* GetHist() { return (TH1*) fHist; }

};


int cb4::fgDebug (0);

//-----------------------------------------------------------------------------
cb4::cb4(const char* Name): TNamed(Name,Name) {
  fHist  = nullptr;
  fHCore = nullptr;

  fFitOpt = "";
  
  fFunc = new TF1(Form("%s_f_cb4"  ,GetName()),cb4::f_cb4 ,-1,1,8);
  fFunc->SetParName(0,"anorm");
  fFunc->SetParName(1,"x0");
  fFunc->SetParName(2,"a");
  fFunc->SetParName(3,"b");
  fFunc->SetParName(4,"alp1");
  fFunc->SetParName(5,"n1");
  fFunc->SetParName(6,"alp2");
  fFunc->SetParName(7,"n2");

  fFunc->SetLineWidth(1);

  fCore = new TF1(Form("%s_f_core",GetName()),cb4::core,-1,1,4);
}

//-----------------------------------------------------------------------------
cb4::~cb4() {
  // delete fGaus;
}

//-----------------------------------------------------------------------------
// dx = (x-P[3])/P[2]
// p[0]*exp(p[1] * (dx − exp(dx)) )
// default Landau distribution for energy losses: P[2] = 2
// decant fit: Landau (-x) constant=260000 mpv=0.5 sigma=0.18
//-----------------------------------------------------------------------------
double cb4::core(double* X, double* P) {
  double x0, f(0), sig, xdx;

  double dx  = (X[0]-P[1]);

  f = P[0]*TMath::Exp(P[2]*(P[3]*dx-TMath::Exp(P[3]*dx)));
  
  return f;
}

//-----------------------------------------------------------------------------
double cb4::f_cb4(double* X, double* P) {
  double x0, dx, abs_dx0, x1, f, sigm, beta1, abs_beta1, beta2, abs_beta2;
  double n1, n2, a, b, A1, B1, alp1, alp2;

  x0         = P[1];
  a          = P[2];
  b          = P[3];
  alp1       = P[4];
  n1         = P[5];
  alp2       = P[6];
  n2         = P[7];

  if (cb4::fgDebug > 0) {
    printf("P[0]=%10.3e x0=%10.3e a=%10.3e b=%10.3e alp1=%10.3e n1=%10.3e alp1=%10.3e n1=%10.3e\n",
           P[0],x0,a,b,alp1,n1,alp2,n2);
  }

  dx        = X[0]-x0;

  if (dx < -alp1) {
    double B1 = -alp1+n1/(a*b*(1-exp(-b*alp1)));
    double A1 = P[0]*exp(a*(-b*alp1-exp(-b*alp1)))*pow(B1+alp1,n1);

    f  = A1/pow(B1-dx,n1);

    if (cb4::fgDebug > 0) {
      printf("[1] : A1,B1,f : %10.3e %10.3e %10.3e \n", A1,B1,f );
    }
  }
  else if (dx < alp2) { 
					// gauss
    f = P[0]*exp(a*(b*dx-exp(b*dx)));

    if (cb4::fgDebug > 0) {
      printf("[2] : f : %10.3e\n", f );
    }
  }
  else {
    double B2 = -alp2+n2/(a*b*(exp(b*alp2)-1));
    double A2 = P[0]*exp(a*(b*alp2-exp(b*alp2)))*pow(B2+alp2,n2);

    f  = A2/pow(B2+dx,n2);
    
    if (cb4::fgDebug > 0) {
      printf("[3] : A2,B2,f : %10.3e %10.3e %10.3e \n", A2,B2,f );
    }
  }
  
  return f;
}

//-----------------------------------------------------------------------------
void cb4::init_parameters(double P0, double P1, double P2, double P3, 
			  double P4, double P5, double P6, double P7) {

  fFunc->SetParameter(0,P0);
  fFunc->SetParameter(1,P1);
  fFunc->SetParError (1,0.1);
  fFunc->SetParameter(2,P2);
  fFunc->SetParameter(3,P3);
  
  fFunc->SetParameter(4,P4);            // alp1
  fFunc->SetParLimits(4,0,100);
  fFunc->SetParError (4,0.1);

  fFunc->SetParameter(5,P5);            // n1
  fFunc->SetParLimits(5,0,100);
  fFunc->SetParError (5,0.1);

  fFunc->SetParameter(6,P6);            // alp2
  fFunc->SetParLimits(6,0,100);
  fFunc->SetParError (6,0.1);

  fFunc->SetParameter(7,P7);            // n2
  fFunc->SetParLimits(7,0,100);
  fFunc->SetParError (7,0.1);
}

//-----------------------------------------------------------------------------
void cb4::fit_histogram(TH1F* Hist, double XMin, double XMax) {
  fHist->GetXaxis()->SetRangeUser(XMin,XMax);
  fHist->Fit(fFunc,fFitOpt.Data(),"",XMin,XMax);
}

//-----------------------------------------------------------------------------
void cb4::fit(const char* File, const char* Module, const char* Hist, 
	      double XMin, double XMax) {

  TH1F* h = (TH1F*) gh1(File,Module,Hist);
  fHist  = (TH1F*) h->Clone("h_fit_cb4");

  this->fit_histogram(h,XMin,XMax);
}

//-----------------------------------------------------------------------------
// assume init_parameters is a separate step
//-----------------------------------------------------------------------------
void cb4::fit(double XMin, double XMax) {

  fit_histogram(fHist,XMin,XMax);
//-----------------------------------------------------------------------------
// core  component of the fit, to display
//-----------------------------------------------------------------------------
  fCore->SetMinimum(XMin);
  fCore->SetMaximum(XMax);
  fCore->SetParameter(0,fFunc->GetParameter(0));
  fCore->SetParameter(1,fFunc->GetParameter(1));
  fCore->SetParameter(2,fFunc->GetParameter(2));
  fCore->SetParameter(3,fFunc->GetParameter(3));

  fHCore = (TH1F*) fHist->Clone(Form("%s_h_core",GetName()));

  fHCore->Reset();

  for (int i=1; i<fHist->GetNbinsX(); i++) {
    float x = fHist->GetBinCenter(i);
    float y = fCore->Eval(x);
    fHCore->SetBinContent(i,y);
    fHCore->SetBinError  (i,0);
  }
  
  fHCore->SetLineColor(kBlue+3);
  fHCore->SetFillColor(kBlue+3);
  fHCore->SetFillStyle(3003);
  
  //  fHCore->Draw("same");
  fCore->SetLineColor(kBlue+3);
  fCore->SetFillColor(kBlue+3);
  fCore->SetFillStyle(3003);
  fCore->Draw("same");
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
