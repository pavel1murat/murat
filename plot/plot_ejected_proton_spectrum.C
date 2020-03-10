///////////////////////////////////////////////////////////////////////////////
// the code is cloned from  Mu2eUtilities/src/EjectedProtonSpectrum.cc
// 'e' is the kinetic energy
//
// call:  plot_ejected_proton_spectrum("e") 
//   or 
//        plot_ejected_proton_spectrum("p")
///////////////////////////////////////////////////////////////////////////////
#include "Mu2eUtilities/inc/EjectedProtonSpectrum.hh"

double GetEjectedProtontWeight(double* X, double* Par) {

    //taken from GMC
    //
    //   Ed Hungerford  Houston University May 17 1999
    //   Rashid Djilkibaev New York University (modified) May 18 1999
    //
    //   e - proton kinetic energy (MeV)
    //   p - proton Momentum (MeV/c)
    //
    //   Generates a proton spectrum similar to that observed in
    //   mu capture in Si.  JEPT 33(1971)11 and PRL 20(1967)569

    //these numbers are in MeV!!!!
    static const double emn  = 1.4; // replacing par1 from GMC
    static const double par2 = 1.3279;
    static const double par3 = 17844.0;
    static const double par4 = .32218;
    static const double par5 = 100.;
    static const double par6 = 10.014;
    static const double par7 = 1050.;
    static const double par8 = 5.103;

    static double MP = 938. ; 
    
    double spectrumWeight;

    double e, p, w1(1.);

    if (Par[0] == 0) e = X[0];
    else {
      p  = X[0];
      e  = p*p/(2*MP);
      w1 = p/MP;
    }

    if (e >= 20) {
        spectrumWeight=par5*exp(-(e-20.)/par6);
    }
    else if (e >= 8.0 && e <= 20.0) {
        spectrumWeight=par7*exp(-(e-8.)/par8);
    }
    else if (e > emn) {
        double xw=(1.-emn/e);
        double xu=std::pow(xw,par2);
        double xv=par3*exp(-par4*e);
        spectrumWeight=xv*xu;
    }
    else {
        spectrumWeight = 0.;
    }

    return w1*Par[1]*spectrumWeight;
  }


//-----------------------------------------------------------------------------
// plot energy / momentum distributions
//-----------------------------------------------------------------------------
TF1* ep_fun;

void plot_ejected_proton_spectrum(const char* Variable = "e") {

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
