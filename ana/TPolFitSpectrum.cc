///////////////////////////////////////////////////////////////////////////////
// for an arbitrary energy spectrum, define a number of energy ranges and
// fit the spectrum within each range with a polynomial or a similar function
// ranges overlap by 0.2 MeV.
// store fit results for each range
// when using the parameterization,  first , determine the range number,
// then use the stored parameterization within this range to calculate dN/dE
// range=i if (e > emin[i]+0.1) && (e < emax[i]-0.1)
// 
// example:
/*
   TPolFitSpectrum* x = new TPolFitSpectrum("DIO_RC");
   x->main_pol_dio();
   x->PrintSpectrum(0.05,105.05,0.1); > conversion_electrons_rc.tbl
*/

#include "ana/TPolFitSpectrum.hh"
#include "TMath.h"
#include "TFile.h"

// ClassImp(TPolFitSpectrum)

namespace{
//-----------------------------------------------------------------------------
// local for the file
// Mode = 0: dx is calculated starting from XMin
// Mode = 1: dx is calculated starting from XMax
//-----------------------------------------------------------------------------
					// DIO ranges
  double dio_rc_range[] = { 
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

  double cnv_rc_range[] = { 
   1,   0.1,   10.0,
   1,  10.0,   50.0,
   1,  50.0,   86.0,
   1,  86.0,   97.0,
   1,  97.0,  103.0,
   1, 103.0,  103.9,
   1, 103.9,  104.75,
   1, 104.75, 104.85,
   1, 104.85, 104.95,
   1, 104.95, 105.05,
   -1
  };
};


TPolFitSpectrum*              TPolFitSpectrum::fgFit (0);

// TPolFitSpectrum::Range_t*     TPolFitSpectrum::fRange(NULL);
double                        TPolFitSpectrum::fP   [100][10];
int                           TPolFitSpectrum::fgMode;
int                           TPolFitSpectrum::fNRanges;
int                           TPolFitSpectrum::fgDebugMode;
double                        TPolFitSpectrum::X0;
//-----------------------------------------------------------------------------
// Name: 'dio_rc'
//       'cnv_rc'
//-----------------------------------------------------------------------------
TPolFitSpectrum::TPolFitSpectrum(const char* Name): TNamed(Name,Name) {

  TString name(Name);
  name.ToUpper();

  if (fgFit != 0) {
    printf("attempt to reinitialize TPolFitSpectrum. BAIL OUT\n");
  }

  fgFit = this;

  if (name == "DIO_RC") {
    for (int i=0; dio_rc_range[3*i] > 0; i++) fNRanges++;
  }
  else if (name == "CNV_RC") {
    for (int i=0; cnv_rc_range[3*i] > 0; i++) fNRanges++;
  }

  printf(" nranges = %3i\n",fNRanges);

  fRange = new Range_t[fNRanges];
  
  for (int i=0; i<fNRanges; i++) {
    if (name == "DIO_RC") {
      fRange[i].fMode = dio_rc_range[3*i  ];
      fRange[i].fXMin = dio_rc_range[3*i+1];
      fRange[i].fXMax = dio_rc_range[3*i+2];
    }
    else if (name == "CNV_RC") {
      fRange[i].fMode = cnv_rc_range[3*i  ];
      fRange[i].fXMin = cnv_rc_range[3*i+1];
      fRange[i].fXMax = cnv_rc_range[3*i+2];
    }
  }

  init_hist();

  fgDebugMode = 0;

  fStep = 0.1;  // 100 KeV

}


//-----------------------------------------------------------------------------
TPolFitSpectrum::TPolFitSpectrum(const char* Name, const TH1F* Hist, const Range_t* R): TNamed(Name,Name) {

  TString name(Name);
  name.ToUpper();

  if (fgFit != 0) {
    printf("attempt to reinitialize TPolFitSpectrum. BAIL OUT\n");
  }

  fgFit = this;

//   if (name == "DIO_RC") {
//     for (int i=0; dio_rc_range[3*i] > 0; i++) fNRanges++;
//   }
//   else if (name == "CNV_RC") {
//     for (int i=0; cnv_rc_range[3*i] > 0; i++) fNRanges++;
//   }

  fRange   = 0;

  fNRanges = 0;
  while (R[fNRanges].fMode >= 0) fNRanges++;

  printf(" nranges = %3i\n",fNRanges);
  if (fNRanges > 0)  fRange = new Range_t[fNRanges];

  for (int i=0; i<fNRanges; i++) {
    fRange[i].fMode = R[i].fMode;
    fRange[i].fXMin = R[i].fXMin;
    fRange[i].fXMax = R[i].fXMax;
  }
//-----------------------------------------------------------------------------
// initialize the histogram
//-----------------------------------------------------------------------------
  fHist = (TH1F*) Hist->Clone("Hist");

  //  double qn = fHist->Integral();
					// so far, assume constant bin
  fStep = fHist->GetBinWidth(1);


  fgDebugMode = 0;
}


//-----------------------------------------------------------------------------
TPolFitSpectrum::~TPolFitSpectrum() {
  delete fRange;
}

//-----------------------------------------------------------------------------
double TPolFitSpectrum::f_pol_dio(double* X, double* P) {

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
double TPolFitSpectrum::f_pol(double* X, double* P) {

  double dx;

  if (fgMode == 1) dx = X[0]-X0;
  else             dx = X0-X[0];

  double f     = P[0] + P[1]*dx+P[2]*dx*dx+P[3]*dx*dx*dx;

  // printf("X[0],P[0],P[1],P[2],P[3],f = %10.3f %12.5e %12.5e %12.5e %12.5e %12.5e\n",
  // 	 X[0],P[0],P[1],P[2],P[3],f);

  return f;
}



//-----------------------------------------------------------------------------
// Mode = 0: dx is calculated starting from XMin
// Mode = 1: dx is calculated starting from XMax
//-----------------------------------------------------------------------------
void TPolFitSpectrum::fit_pol(int Mode, double XMin, double XMax) {

  TSpline3  *s;
  
  fNPar = 4;
  fFunc = new TF1("f",f_pol,XMin,XMax,fNPar);
  //  fFunc->SetNpx(1);

  int    ib;
  double y0, y1, y2, y22, y12, dydx, dydx2, dx, ddx, xmin;
  
  if (fNRanges <= 0) {
    printf(">>> TPolFitSpectrum::fit_pol : fNRanges <=0 , RETURN\n");
    return;
  }

  fgMode = Mode;
  dx     = XMax-XMin;
  s      = new TSpline3(fHist);

//-----------------------------------------------------------------------------
// if the range is just one bin (1+1+1), take the histogram value in that bin
//-----------------------------------------------------------------------------
  if (dx < 4*fStep) {
    X0   = (XMax+XMin)/2.;
    xmin = fHist->GetXaxis()->GetXmin();
    ib   = (X0-xmin)/fStep+1;
    y0   = fHist->GetBinContent(ib);
    fFunc->SetParameter(0,y0);
    fFunc->SetParameter(1,0);
    fFunc->SetParameter(2,0);
    fFunc->SetParameter(3,0);
    return;
  }
  else if      (fgMode == 1) {
    X0  = XMin;
    y0  = s->Eval(X0);

    ddx = 0.5;
    y1  = s->Eval(XMin);  
    y12 = s->Eval(XMin+ddx);  
    y2  = s->Eval(XMax);  
    y22 = s->Eval(XMax-ddx);  

    dydx  = (y2-y1)/dx; 
    dydx2 = ((y12-y1)/ddx - (y2-y22)/ddx)/dx;
    //    dydx3 = dydx2/dx/10;

    
    // printf("dx,X0, y0,y1,y2: %10.3f %10.3f %12.5e %12.5e %12.5e\n",dx,X0,y0,y1,y2);
	 

    
  // printf("dydx = %10.5f\n",dydx);
  
    fFunc->SetParameter(0,y0);
    fFunc->SetParameter(1,1);
    fFunc->SetParameter(2,dydx);
    fFunc->SetParameter(3,dydx2);
  }
  else if (fgMode == 2) { 



    X0  = XMax;
    y0  = s->Eval(X0);  

    ddx = 0.5;
    y1  = s->Eval(XMin);  
    y12 = s->Eval(XMin+ddx);  
    y2  = s->Eval(XMax);  
    y22 = s->Eval(XMax-ddx);  

    dydx  = (y2-y1)/dx; 
    dydx2 = ((y12-y1)/ddx - (y2-y22)/ddx)/dx;
    //    dydx3 = dydx2/dx/10;

    fFunc->SetParameter(0,y0);
    fFunc->SetParameter(1,1);
    fFunc->SetParameter(2,dydx);
    fFunc->SetParameter(3,dydx2);
  }
  else    {
    printf("ERROR: undefined Mode = %i\n",Mode);
    return;
  }

  fHist->GetXaxis()->SetRangeUser(XMin-1,XMax+1);
  fHist->Fit(fFunc,"","",XMin,XMax);
}

//-----------------------------------------------------------------------------
// Mode = 0: dx is calculated starting from XMin
// Mode = 1: dx is calculated starting from XMax
//-----------------------------------------------------------------------------
void TPolFitSpectrum::fit_pol_dio(int Mode, double XMin, double XMax) {

  TSpline3  *s;
  
  fNPar = 6;
  fFunc = new TF1("f",f_pol_dio,XMin,XMax,fNPar);
  //  fFunc->SetNpx(1);

  int    ib;
  double y0, y1, y2, y22, y12, dydx, dydx2, dydx3, dx, ddx, xmin;
  
  fgMode = Mode;
  dx     = XMax-XMin;
  s      = new TSpline3(fHist);

//-----------------------------------------------------------------------------
// if the range is just one bin (1+1+1), take the histogram value in that bin
//-----------------------------------------------------------------------------
  if (dx < 4*fStep) {
    X0   = (XMax+XMin)/2.;
    xmin = fHist->GetXaxis()->GetXmin();
    ib   = (X0-xmin)/fStep+1;
    y0   = fHist->GetBinContent(ib);
    fFunc->SetParameter(0,y0);


    fFunc->SetParameter(1,0);
    fFunc->SetParameter(2,0);
    fFunc->SetParameter(3,0);
    fFunc->SetParameter(4,0);
    fFunc->SetParameter(5,0);
    return;
  }
  else if      (fgMode == 1) {
    X0  = XMin;
    y0  = s->Eval(X0);  

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

    X0  = XMax;
    y0  = s->Eval(X0);  

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

  fHist->GetXaxis()->SetRangeUser(XMin-1,XMax+1);
  fHist->Fit(fFunc,"","",XMin,XMax);
}


//-----------------------------------------------------------------------------
// photos_4.hist: DIO spectrum with radiative corrections
// photos_7.hist: CE spectrum with radiative corrections
// normalize the histogram to 10
// assume bin of 0.1, such that the integral, accounting for the bin size,
// would be equal to 1.
//-----------------------------------------------------------------------------
int TPolFitSpectrum::init_hist() {
  TFile*  f;
  TString name;

  name = GetName();
  name.ToUpper();
  
  if      (name == "DIO_RC") f = TFile::Open("~/hist/photos/muon_decay/photos_4.hist");
  else if (name == "CNV_RC") f = TFile::Open("~/hist/photos/muon_decay/photos_2.hist");
  else {
    printf(" NAME: %s is not defined. BAIL OUT\n",name.Data());
    return -1;
  }
  
  TH1F* hist = (TH1F*) f->Get("Ana/PhotosAna/Hist/genp_0/e_1");

  fHist = (TH1F*) hist->Clone("Hist");

  double qn = fHist->Integral();
					// so far, assume constant bin
  fStep = fHist->GetBinWidth(1);

  int nb = fHist->GetNbinsX();

  for (int i=1; i<=nb; i++) {
    double y = fHist->GetBinContent(i)*(10/qn);
    double err = y/200.;
    fHist->SetBinContent(i,y  );
    fHist->SetBinError  (i,err);
  }
  
  return 0;
}



//-----------------------------------------------------------------------------
// main function: returns parameterization of the DIO spectrum
//-----------------------------------------------------------------------------
double TPolFitSpectrum::dio_energy_spectrum(double E) {
// find range  double x, fwhm; 



  double f(0);

  Range_t* r;

  //  printf(" dio_energy_spectrum: E = %12.4f\n",E);

  int ir = -1;

  //  double step = fgFit->fStep;
  for (int i=0; i<fNRanges; i++) {
    //    printf("i, xmin, xmax = %3i %12.5f %12.5f\n",i,range[3*i+1],range[3*i+2]);

    r = fgFit->fRange+i;
    
    if (E >= r->fXMin /*+step/2*/) {
      if (E < r->fXMax/*-step/2*/) {
	ir = i;
      }
      else {
	if (i == fNRanges-1) {
					// this 'if' is a DIO thing - why is it needed ?
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

  fgMode = fgFit->fRange[ir].fMode;

  if      (fgMode == 1) X0 = fgFit->fRange[ir].fXMin;
  else if (fgMode == 2) X0 = fgFit->fRange[ir].fXMax;
  
  f      = f_pol_dio(&E, fP[ir]);

  if (fgDebugMode > 0) {
    printf(" dio_energy_spectrum: E = %10.3f range = %i, f = %12.5g\n",E,ir,f);
  }

 END:;
  return f;
}

//-----------------------------------------------------------------------------
// main function: returns parameterization of the DIO spectrum
//-----------------------------------------------------------------------------
double TPolFitSpectrum::spectrum(double E) {
// find range

  double f(0);

  Range_t* r;

  int ir = -1;

  for (int i=0; i<fNRanges; i++) {
    //    printf("i, xmin, xmax = %3i %12.5f %12.5f\n",i,range[3*i+1],range[3*i+2]);

    r = fgFit->fRange+i;
    
    if (E >= r->fXMin ) {
      if (E < r->fXMax) {
	ir = i;
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

  fgMode = fgFit->fRange[ir].fMode;

  if      (fgMode == 1) X0 = fgFit->fRange[ir].fXMin;
  else if (fgMode == 2) X0 = fgFit->fRange[ir].fXMax;
  
  f      = f_pol(&E, fP[ir]);

  if (fgDebugMode > 0) {
    printf(" spectrum: E = %10.3f range = %i, f = %12.5g\n",E,ir,f);
  }

 END:;
  return f;
}



//-----------------------------------------------------------------------------
int TPolFitSpectrum::GetFWHM(double* XMax, double* Fwhm) {

//-----------------------------------------------------------------------------
// step 1: find maximum xx - x , corresponding to the maximum, fmax - the function value 
//-----------------------------------------------------------------------------
  int      found(0); 
  double   fm, fmax(-1), xx, xext[2], p[3], q[2], val;
  double   xtot;

  for (int ir=0; ir<fNRanges; ir++) {
    // 3rd order polynomial - find derivative

    double xmin = fRange[ir].fXMin;
    double xmax = fRange[ir].fXMax;

    p[0] =   fP[ir][1];
    p[1] = 2*fP[ir][2];
    p[2] = 3*fP[ir][3];
//-----------------------------------------------------------------------------
// solve quadratic equation to find points where the derivative is equal to zero
//-----------------------------------------------------------------------------
    int nsol = 0;
    if (p[2] == 0) {
      nsol    = 1;
      xext[0] = -p[0]/p[1];
    }
    else {
      double a = p[0]/p[2];
      double b = p[1]/p[2];

      double det = b*b/4-a;
      if (det > 0) {
	nsol = 2;
	xext[0] = -b/2+sqrt(det);
	xext[1] = -b/2-sqrt(det);
      }
      else if (det == 0) {
	nsol    = 1;
	xext[0] = -b/2;
      }
    }

    double xt;
    for (int i=0; i<nsol; i++) {
      if      (fRange[ir].fMode == 1) xt = fRange[ir].fXMin+xext[i];
      else {
	printf(" TROUBLE #2\n");
      }

      if ((xt >= xmin) && (xt < xmax)) {
					// extremum is within the range, look at the 2nd derivative
	q[0] = p[1];
	q[1] = 2*p[2];

	val = q[0]+q[1]*xext[i];
	if (val < 0) {
					// maximum, assume the only one
	  xx    = xext[i];
	  fm    = fP[ir][0]+fP[ir][1]*xx+fP[ir][2]*xx*xx+fP[ir][3]*xx*xx*xx;

	  if (fm > fmax) {
	    found = 1;
	    xtot  = xt;
	    fmax  = fm;
	    break;
	  }
	}
      }
    }
  }
  
  if (found == 0) {
					// in trouble
    printf(" GetFWHM : IN TROUBLE, EXIT\n");
  }

  printf(" fmax = %10.3f\n",fmax);
//-----------------------------------------------------------------------------
// maximum found, look for the full width at half max
//-----------------------------------------------------------------------------
  double xw1(-1), xw2(0), fwhm(-1);

  int nbinx = fHist->GetNbinsX();
  for (int i=1; i<nbinx; i++) {
					// find range for a given bin
    double x1 = fHist->GetBinCenter(i);
    double x2 = fHist->GetBinCenter(i+1);
    
    double y1 = spectrum(x1);
    double y2 = spectrum(x2);
//-----------------------------------------------------------------------------
// in principle, should use Cardano formula
// for now - a kludge
//-----------------------------------------------------------------------------
    if ((y1 < fmax/2) && (y2 > fmax/2)) {
					// find the first point 
      xw1 = x1+(fmax/2-y1)/(y2-y1)*(x2-x1);
    }
    else if ((y1 > fmax/2) && (y2 < fmax/2)) {
      xw2 = x1+(fmax/2-y1)/(y2-y1)*(x2-x1);
    }
  }

  if (xw1 < xw2) fwhm = xw2-xw1;

  *XMax = xtot;
  *Fwhm = fwhm; 
  
  return 0;
}






//-----------------------------------------------------------------------------
// print spectrum in a format required by Rob's .tbl file
// not sure why the existing printout starts from the high end
// so far, only makes sense for the interval [0,105.0] - normalize
//-----------------------------------------------------------------------------
void TPolFitSpectrum::PrintSpectrum(double EMin, double EMax, double Step) {
  double e, f;

  e = EMax;
					// original histogram was normalized to 10, step - 0.1,
					// so the integral should be normalized to 1
  
  double scale = 1./fFull->Integral(EMin,EMax);
  
  while (e >= EMin) {
    f = dio_energy_spectrum(e)*scale;
    printf("%6.2f  %12.5e\n",e,f);
    e -= Step;
  }
}


//-----------------------------------------------------------------------------
// at this point the parameters are fixed
//-----------------------------------------------------------------------------
double TPolFitSpectrum::f_dio_energy_spectrum(double* X, double* P) {
  return dio_energy_spectrum(X[0]);
}

//-----------------------------------------------------------------------------
// at this point the parameters are fixed
//-----------------------------------------------------------------------------
double TPolFitSpectrum::f_spectrum(double* X, double* P) {
  return spectrum(X[0]);
}

//-----------------------------------------------------------------------------
int TPolFitSpectrum::main_pol_dio() {

  Range_t*  r;
  int       mode;
  double    xmin, xmax;
  int       rc(0);

  rc  = init_hist();

  if (rc < 0) return rc;
  
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

  fHist->Draw();
  fFull->Draw("same");

  return 0;
}



//-----------------------------------------------------------------------------
int TPolFitSpectrum::main_pol() {

  Range_t*  r;
  int       mode;
  double    xmin, xmax;
  int       rc(0);

  for (int i=0; i<fNRanges; i++) {
    // for (int i=4; i<5; i++) {
    r = fRange+i;
    
    mode = r->fMode;
    xmin = r->fXMin;
    xmax = r->fXMax;

    printf(" i, mode, xmin, xmax: %3i %3i %10.3f %10.3f\n",i,fNRanges,xmin,xmax);

    fit_pol(mode,xmin,xmax);
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

  fFull = new TF1("full",f_spectrum,fRange[0].fXMin,fRange[fNRanges-1].fXMax,0);

  fHist->Draw();
  fFull->Draw("same");

  return rc;
}


