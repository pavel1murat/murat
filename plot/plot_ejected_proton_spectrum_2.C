///////////////////////////////////////////////////////////////////////////////
// the code from  Mu2eUtilities/src/EjectedProtonSpectrum.cc
// 'e' is the kinetic energy in MeV
//
// call:  plot_ejected_proton_spectrum("e") 
//   or 
//        plot_ejected_proton_spectrum("p")
///////////////////////////////////////////////////////////////////////////////
#include "TCanvas.h"
#include "TF1.h"
#include "Mu2eUtilities/inc/EjectedProtonSpectrum.hh"

double GetEjectedProtontWeight(double* X, double* Par) {
  static double MP = 938. ; 
  double e, p, w1(1.);

  if (Par[0] == 0) {
    e = X[0];
  }
  else {
    p = X[0];
    e = p*p/(2*MP);
    w1 = p/MP;
  }
  
  double w = w1*Par[1]*mu2e::EjectedProtonSpectrum::getWeight(e);
  
  return w;
}


//-----------------------------------------------------------------------------
// plot energy and momentum distributions
//-----------------------------------------------------------------------------
TF1* ep_fun;

void plot_ejected_proton_spectrum_2(const char* Variable = "e") {

  TString var = Variable;

  var.ToLower();

  if (var == "e") {

    ep_fun = new TF1("ep_fun",GetEjectedProtontWeight,0,100,2);
    ep_fun->SetParameter(0,0);
    ep_fun->SetParameter(1,1);
    ep_fun->SetParameter(1,0.05/ep_fun->Integral(0,100));
    ep_fun->SetNpx(1000);

    ep_fun->SetTitle("Ejected proton energy");
    ep_fun->GetXaxis()->SetTitle("E, MeV");

  }
  else if (var == "p") {

    ep_fun = new TF1("ep_fun",GetEjectedProtontWeight,0,1000,2);
    ep_fun->SetParameter(0,1);
    ep_fun->SetParameter(1,1);
    ep_fun->SetParameter(1,0.05/ep_fun->Integral(0,1000));
    ep_fun->SetNpx(1000);

    ep_fun->SetTitle("Ejected proton momentum");
    ep_fun->GetXaxis()->SetTitle("P, MeV/c");
  }
  
  ep_fun->Draw();
}
