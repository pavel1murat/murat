///////////////////////////////////////////////////////////////////////////////
// integrate radiatively corrected electron spectrum from mu2e-7615 page 12
// from emin=0.7 MeV up to emax=mmu-de,
// where de is the call parameter
// use GNU scientific library
//
// as the highest bin is undefined, print the value of 1-integral(emin,emax)
//
//    de     |  0.1        1.5      2.0
// ----------+------------------------------
// mu2e-7615 | 0.861     0.923     0.930
// -----------------------------------------
// this int  | 0.85966   0.92252   0.92910
//
// to run from ROOT prompt:
// -------
// .L murat/scripts/test_rad_corr_gsl_integration.C+
// test_rad_corr_gsl_integration()
///////////////////////////////////////////////////////////////////////////////
#include "stdio.h"
#include "gsl/gsl_integration.h"

gsl_function F;

struct my_f_params {
  double a;
  double b;
  double c;
};

struct my_f_params params = { 3.0, 2.0, 1.0 };

namespace {
  double me    = 0.511; // MeV
  double mmu   = 105.6; // 
  double alpha = 1./137.;
};

//-----------------------------------------------------------------------------
double my_f (double E, void * p) {
  struct my_f_params * params = (struct my_f_params *)p;
  double a = (params->a);
  double b = (params->b);
  double c = (params->c);

  double f = (1./mmu)*(alpha / (2*M_PI))*(log(4*E*E/me/me)-2.)*((E*E+mmu*mmu)/mmu/(mmu-E));
  
  return  f;
}

//-----------------------------------------------------------------------------
// need log(4*E*E/me/me) > 2: Emin = 0.7 MeV
// integrate up to Emax = mmu-de
//-----------------------------------------------------------------------------
void test_rad_corr_gsl_integration(double de = 0.1) {
  
  F.function = &my_f;
  F.params   = &params;

  size_t limit  = 1000;
  double epsabs = 0.001;
  double epsrel = 0.001;
  
  gsl_integration_workspace * ws = gsl_integration_workspace_alloc(10000);

  double result, abserr;

  double emin = 0.7;
  double emax = mmu-de;

  gsl_integration_qags(&F,
		       emin, emax,
		       epsabs,
		       epsrel,
		       limit,
		       ws,
		       &result,
		       &abserr);

  printf("enemoe result=%12.5f\n",1-result);
}
