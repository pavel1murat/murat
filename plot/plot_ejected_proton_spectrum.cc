//-----------------------------------------------------------------------------
// an example function referencing Mu2e Offline code which can be compiled
// on the fly
//-----------------------------------------------------------------------------
#include "TCanvas.h"
#include "TF1.h"
#include "Mu2eUtilities/inc/EjectedProtonSpectrum.hh"


//-----------------------------------------------------------------------------
double fun(double* X, double* P) {
  return mu2e::EjectedProtonSpectrum::getWeight(X[0]);
}


//-----------------------------------------------------------------------------
void plot_ejected_proton_spectrum() {

  TCanvas* c = new TCanvas("c","c",1000,800);

  TF1* f = new TF1("f1",fun,0,100,0);

  f->Draw();
}
