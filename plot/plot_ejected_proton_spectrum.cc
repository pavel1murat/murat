//-----------------------------------------------------------------------------
// an example function referencing Mu2e Offline code which can be compiled
// on the fly
//-----------------------------------------------------------------------------
#include "TCanvas.h"
#include "TAxis.h"
#include "TF1.h"
#include "Offline/Mu2eUtilities/inc/EjectedProtonSpectrum.hh"
#include "murat/plot/murat_plot_functions.hh"

TF1* ep_fun;

//-----------------------------------------------------------------------------
// P[0]=0: energy, 
//     =1: momentum
// P[1]  : normalization of the integral
//-----------------------------------------------------------------------------
double GetEjectedProtontWeight(double* X, double* P) {
  static double MP = 938. ; 

  double e, p, w1(1.);

  if (P[0] == 0) {
    e = X[0];
  }
  else {
    p  = X[0];
    e  = p*p/(2*MP);
    w1 = p/MP;
  }

  return w1*P[1]*mu2e::EjectedProtonSpectrum::getWeight(e);
}


//-----------------------------------------------------------------------------
// ProtonsPerCapture : default = 0.05, see the include file
//-----------------------------------------------------------------------------
void plot_ejected_proton_spectrum(const char* Variable, double ProtonsPerCapture, const char* Opt) {

  TString var = Variable;

  var.ToLower();

  if (var == "e") {

    ep_fun = new TF1("ep_fun",GetEjectedProtontWeight,0,100,2);
    ep_fun->SetParameter(0,0);
    ep_fun->SetParameter(1,1);
    ep_fun->SetParameter(1,ProtonsPerCapture/ep_fun->Integral(0,100));
    ep_fun->SetNpx(1000);

    ep_fun->SetTitle("Ejected proton energy");
    ep_fun->GetXaxis()->SetTitle("E, MeV");

  }
  else if (var == "p") {

    ep_fun = new TF1("ep_fun",GetEjectedProtontWeight,0,1000,2);
    ep_fun->SetParameter(0,1);
    ep_fun->SetParameter(1,1);
    ep_fun->SetParameter(1,ProtonsPerCapture/ep_fun->Integral(0,1000));
    ep_fun->SetNpx(1000);

    ep_fun->SetTitle("Ejected proton momentum");
    ep_fun->GetXaxis()->SetTitle("P, MeV/c");
  }
  
  ep_fun->Draw(Opt);
}
