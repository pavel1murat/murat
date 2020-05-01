//-----------------------------------------------------------------------------
// test lognormal function - used for parameterization of 
// the beam fluctuations
//-----------------------------------------------------------------------------
#include "TF1.h"
#include "math.h"

//-----------------------------------------------------------------------------
double f_lgn(double* X, double* P) {

  double mu    = P[0];
  double sigma = P[1];
  double x     = X[0];
  double c     = 1./sqrt(2*M_PI)/sigma;

  double dx    = (log(x)-mu)/sigma;
  double f     = c*exp(-dx*dx/2);

  return f;
}

//-----------------------------------------------------------------------------
// default values of Mu and Sigma  - from Dave's mu2e-17192:
// https://mu2e-docdb.fnal.gov/cgi-bin/private/ShowDocument?docid=17192
// assume XMin=0
//-----------------------------------------------------------------------------
int plot_lognormal(double Mu=0.9047,double Sigma=0.3814,double XMax=20) {

  TF1* fun = new TF1("fun",f_lgn,0,100,2);

  fun->SetParameter(0,Mu);
  fun->SetParameter(1,Sigma);

  double xmin = 0;
  double xmax = XMax;
  int nb      = 1000;

  double bw = (xmax-xmin)/nb;

  TH1F* h = new TH1F("h","lognormal",nb,xmin,xmax);

  for(int i=1; i<=1000; i++) {
    double x = (i+0.5)*bw;
    double y = fun->Eval(x);
    h->SetBinContent(i,y);
  }

  h->Draw();

  return 0;
}
