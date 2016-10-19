///////////////////////////////////////////////////////////////////////////////
// 2016-10-19: 
// look at the scale miscalibration due to the use of the wrong momentum spectrum
// radiative corrections change the shape - check what happens to the scale if 
// LO spectrum is used to fit the LL one
//-----------------------------------------------------------------------------


#include "murat/plot/smooth.hh"

namespace {

  TF1* func;
}

//-----------------------------------------------------------------------------
// P[0] : vertical scale
// P[1] : horizontal scale
// f =P[0] * f ((1+P[1])*X[0])
//-----------------------------------------------------------------------------
double fit_f(double* X, double* P) {

  double x = X[0]*(1+P[1]);

  double f = P[0]*func->Eval(x);

  return f;
}


//-----------------------------------------------------------------------------
void fit_dio_rc_spectrum(double EMin = 40., double EMax = 55.) {

  TH1F* h1 = (TH1F*) gh1("hist/d50s6000.track_comp_use_mva.hist","TrackComp","trk_3/p0")->Clone("hd");

  TH1F* hg = (TH1F*) gh2("hist/g50s6000.track_comp_use_mva.hist","TrackComp","trk_3/p0")->Clone("hg");

  h1->Rebin(2);
  smooth* s = new smooth(h1,30,60);
  
  func = s->GetFunc();

  TF1* ff = new TF1("ff",fit_f,30,60,2);

  ff->SetParameters(0.5,0.);

  hg->Fit(ff,"","",EMin,EMax);

  
}
