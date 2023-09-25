///////////////////////////////////////////////////////////////////////////////
// example:
/*
   TFitDioRcSpectrum rf;
   fr.main_pol_dio();
   fr.PrintSpectrum();
*/

#include "ana/TFitDioRcSpectrum.hh"
#include "TMath.h"
#include "TFile.h"

namespace{
//-----------------------------------------------------------------------------
// local for the file
// Mode = 0: dx is calculated starting from XMin
// Mode = 1: dx is calculated starting from XMax
//-----------------------------------------------------------------------------
  double range[] = { 
   1,   0.1,   2.1,
   1,   1.9,  12.1,
   1,  11.9,  16.1,
   1,  15.9,  40.1,
   1,  39.9,  50.1,
   1,  49.9,  56.1,
   1,  55.9,  61.1,
   1,  60.9,  67.1, 
   1,  66.9,  73.1,
   1,  72.9,  79.1,
   1,  78.9,  85.1,
   1,  84.9,  90.1,
   1,  89.9,  95.1,
   1,  94.9,  99.1,
   1,  98.9, 101.1,
   1, 100.9, 102.6,
   2, 102.4, 104.1,
   2, 103.9, 104.9,
   -1
  };
}

TFitDioRcSpectrum::Range_t*   TFitDioRcSpectrum::fRange(NULL);
double                        TFitDioRcSpectrum::fP   [100][10];
int                           TFitDioRcSpectrum::fgMode;
int                           TFitDioRcSpectrum::fNRanges;
int                           TFitDioRcSpectrum::fgDebugMode;
double                        TFitDioRcSpectrum::X0;

//-----------------------------------------------------------------------------
TFitDioRcSpectrum::TFitDioRcSpectrum() {
  
  for (int i=0; range[3*i] > 0; i++) fNRanges++;

  printf(" nranges = %3i\n",fNRanges);

  fRange = new Range_t[fNRanges];
  
  for (int i=0; i<fNRanges; i++) {
    fRange[i].fMode = range[3*i  ];
    fRange[i].fXMin = range[3*i+1];
    fRange[i].fXMax = range[3*i+2];
  }

  init_hist();

  fgDebugMode = 0;
}


//-----------------------------------------------------------------------------
TFitDioRcSpectrum::~TFitDioRcSpectrum() {
  delete fRange;
}

//-----------------------------------------------------------------------------
double TFitDioRcSpectrum::f_pol_dio(double* X, double* P) {

  double dx;

  if (fgMode == 1) dx = X[0]-X0;
  else             dx = X0-X[0];
  //  double dx    = fabs(X[0]-X0);

  double alpha = P[1];
  
  double f     = P[0] + TMath::Power(dx,alpha)*(P[2]+P[3]*dx+P[4]*dx*dx+P[5]*dx*dx*dx);

  // printf("X[0],P[0],P[1],P[2],P[3],f = %10.3f %12.5e %12.5e %12.5e %12.5e %12.5e\n",
  // 	 X[0],P[0],P[1],P[2],P[3],f);

  return f;
}



//-----------------------------------------------------------------------------
// Mode = 0: dx is calculated starting from XMin
// Mode = 1: dx is calculated starting from XMax
//-----------------------------------------------------------------------------
void TFitDioRcSpectrum::fit_pol_dio(int Mode, double XMin, double XMax) {

  fNPar = 6;
  fFunc = new TF1("f",f_pol_dio,XMin,XMax,fNPar);
  //  fFunc->SetNpx(1);

  double y0, y1, y2, y22, y12, dydx, dydx2, dydx3, dx, ddx;
  
  fgMode = Mode;
  if      (fgMode == 1) {
    X0 = XMin;

    dx = XMax-XMin;
    //    double xm = (XMax+XMin)/2;

    TSpline3* s = new TSpline3(Hist);
    
    y0   = s->Eval(X0);  

    ddx = 0.5;
    y1  = s->Eval(XMin);  
    y12 = s->Eval(XMin+ddx);  
    y2  = s->Eval(XMax);  
    y22 = s->Eval(XMax-ddx);  

    dydx  = (y2-y1)/dx; 
    dydx2 = ((y12-y1)/ddx - (y2-y22)/ddx)/dx;
    dydx3 = dydx2/dx/10;

    
    // printf("dx,X0, y0,y1,y2: %10.3f %10.3f %12.5e %12.5e %12.5e\n",dx,X0,y0,y1,y2);
	 

    
  // printf("dydx = %10.5f\n",dydx);
  
    fFunc->SetParameter(0,y0);
    fFunc->SetParameter(1,1);
    fFunc->SetParameter(2,dydx);
    fFunc->SetParameter(3,dydx2);
    fFunc->SetParameter(4,dydx3);
    fFunc->SetParameter(5,dydx3/10);
    //    fFunc->SetParameter(4,dydx3);
  }
  else if (fgMode == 2) {

    X0 = XMax;

    TSpline3* s = new TSpline3(Hist);
    
    y0   = s->Eval(X0);  
    dx = XMax-XMin;
    // double xm = (XMax+XMin)/2;

    ddx = 0.5;
    y1  = s->Eval(XMin);  
    y12 = s->Eval(XMin+ddx);  
    y2  = s->Eval(XMax);  
    y22 = s->Eval(XMax-ddx);  

    dydx  = (y2-y1)/dx; 
    dydx2 = ((y12-y1)/ddx - (y2-y22)/ddx)/dx;
    dydx3 = dydx2/dx/10;

    fFunc->SetParameter(0,y0);
    fFunc->SetParameter(1,1);
    fFunc->SetParameter(2,dydx);
    fFunc->SetParameter(3,dydx2);
    fFunc->SetParameter(4,dydx3);
    fFunc->SetParameter(5,dydx3/10);
    //    fFunc->SetParameter(4,dydx3);

  }
  else    {
    printf("ERROR: undefined Mode = %i\n",Mode);
    return;
  }

  Hist->GetXaxis()->SetRangeUser(XMin-1,XMax+1);
  Hist->Fit(fFunc,"","",XMin,XMax);
}


//-----------------------------------------------------------------------------
// get DIO spectrum with radiative corrections
// normalize the histogram to 10
// assume bin of 0.1, such that the integral, accounting for the bin size,
// would be equal to 1.
//-----------------------------------------------------------------------------
int TFitDioRcSpectrum::init_hist() {

  TFile* f   = TFile::Open("~/hist/photos/muon_decay/photos_4.hist");
  TH1F* hist = (TH1F*) f->Get("Ana/PhotosAna/Hist/genp_0/e_1");

  Hist = (TH1F*) hist->Clone("Hist");

  double qn = Hist->Integral();

  int nb = Hist->GetNbinsX();

  for (int i=1; i<=nb; i++) {
    double y = Hist->GetBinContent(i)*(10/qn);
    double err = y/200.;
    Hist->SetBinContent(i,y  );
    Hist->SetBinError  (i,err);
  }
  return 0;
}



//-----------------------------------------------------------------------------
// main function: returns parameterization of the DIO spectrum
//-----------------------------------------------------------------------------
double TFitDioRcSpectrum::dio_energy_spectrum(double E) {
// find range

  double f(0);

  Range_t* r;

  //  printf(" dio_energy_spectrum: E = %12.4f\n",E);

  int ir = -1;
  for (int i=0; i<fNRanges; i++) {
    //    printf("i, xmin, xmax = %3i %12.5f %12.5f\n",i,range[3*i+1],range[3*i+2]);

    r = fRange+i;
    
    if (E >= r->fXMin+0.1) {
      if (E < r->fXMax-0.1) {
	ir = i;
      }
      else {
	if (i == fNRanges-1) {
	  if (E < 104.8) ir = fNRanges-1;
	  else                                         goto END;
	}
      }
    }
    else {
      if (i == 0) {
//-----------------------------------------------------------------------------
// first range
//-----------------------------------------------------------------------------
	if   (E > 0.5) ir = 0;
	else                                           goto END;
      }
    }
  }
//-----------------------------------------------------------------------------
// range found, determine the function value
//-----------------------------------------------------------------------------

  fgMode = fRange[ir].fMode;

  if      (fgMode == 1) X0 = fRange[ir].fXMin;
  else if (fgMode == 2) X0 = fRange[ir].fXMax;
  
  f      = f_pol_dio(&E, fP[ir]);

  if (fgDebugMode > 0) {
    printf(" dio_energy_spectrum: E = %10.3f range = %i, f = %12.5g\n",E,ir,f);
  }

 END:;
  return f;
}

//-----------------------------------------------------------------------------
// print spectrum in a format required by Rob's .tbl file
// not sure why the existing printout starts from the high end
//-----------------------------------------------------------------------------
void TFitDioRcSpectrum::PrintSpectrum(double EMin, double EMax, double Step) {
  double e, f;

  e = EMax;
  
  while (e >= EMin) {
    f = dio_energy_spectrum(e);
    printf("%5.2f  %12.5e\n",e,f);
    e -= Step;
  }
}


//-----------------------------------------------------------------------------
// at this point the parameters are fixed
//-----------------------------------------------------------------------------
double TFitDioRcSpectrum::f_dio_energy_spectrum(double* X, double* P) {
  return dio_energy_spectrum(X[0]);
}

//-----------------------------------------------------------------------------
int TFitDioRcSpectrum::main_pol_dio() {

  Range_t*  r;
  int       mode;
  double    xmin, xmax;

  init_hist();
  
  for (int i=0; i<fNRanges; i++) {
    // for (int i=4; i<5; i++) {
    r = fRange+i;
    
    mode = r->fMode;
    xmin = r->fXMin;
    xmax = r->fXMax;

    printf(" i, mode, xmin, xmax: %3i %3i %10.3f %10.3f\n",i,fNRanges,xmin,xmax);

    fit_pol_dio(mode,xmin,xmax);
//------------------------------------------------------------------------------
// extract the parameter values
//------------------------------------------------------------------------------
    for (int ip=0; ip<fNPar; ip++) {
      fP[i][ip] = fFunc->GetParameter(ip);
      printf(" %12.5e\n",fP[i][ip]);
    }
  }

  printf("done with fitting\n");
//-----------------------------------------------------------------------------
// at this point parameters for all ranges are saved
//-----------------------------------------------------------------------------

  fFull = new TF1("full",f_dio_energy_spectrum,0,105,0);

  Hist->Draw();
  fFull->Draw("same");

  return 0;
}
