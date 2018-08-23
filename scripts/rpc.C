//

#include "TH1F.h"

// Ivano's parameterization from mu2e-665 ? from PiCaptureEffects

double rpc(double E) {

  constexpr double emax  = 138.2;
  constexpr double alpha =   2.691;
  constexpr double gamma =   1.051;
  constexpr double tau   =   8.043;
  constexpr double c0    =   2.741;
  constexpr double c1    =  -0.005;

  double de = emax-E;

  if (de < 0) return 0;

  return pow(emax-E,alpha) * exp(-(emax-gamma*E)/tau) * (c0 + c1*E);
}



//-----------------------------------------------------------------------------
void plot_rpc() {
  

  TH1F* h = new TH1F("h1","RPC mu2e-665 ?",200,0,200);

  for (int i=0; i<200; i++) {
    double e = 0.5+i;

    double w = rpc(e);

    h->Fill(e,w);
  }

  h->Draw("h");
}
