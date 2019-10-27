//

#include "TRandom3.h"
#include "TMath.h"
#include "TF1.h"

// #include "TPoisson.h"

#include "murat/ana/poisson_fluct.hh"

ClassImp(poisson_fluct)




//-----------------------------------------------------------------------------
double poisson_fluct::fun(double* x, double* p) {
  return TMath::Poisson(x[0],p[0]);
}


//-----------------------------------------------------------------------------
poisson_fluct::poisson_fluct() {
  fNTries = 100000000; // 1.e8
}
//-----------------------------------------------------------------------------
// integrate probability of the poisson distribution with the given 'Mean'
// from 0 to N
//-----------------------------------------------------------------------------
int poisson_fluct::poi1(double N, double Mean) {

  TF1 f("f",poisson_fluct::fun,0,20.,1);

  f.SetParameter(0,Mean);

  double s = f.Integral(0,N);

  printf ("prob = %10.6f\n",s);

  return 0;
}

//-----------------------------------------------------------------------------
// sum poisson probabilities for the distribution with given 'Mean' from 0 to N
//-----------------------------------------------------------------------------
int poisson_fluct::poi(double N, double Mean) {

  double sum(0);

  for (int i=0; i<=N; i++) {
    sum += TMath::Poisson(i,Mean);
  }

  printf ("prob = %10.6f\n",sum);

  return 0;
  
}


//-----------------------------------------------------------------------------
int poisson_fluct::l95(double Signal) {
  
  // starting value

  double mean = Signal+3.;

  double prob;
  int nj;

  // calculate probability to observe less than signal using MC

  TRandom3* r = new TRandom3();

  double prr = 0;
  
  double y;

  int ntries = 100;

  for (int i=0;i<ntries; i++) {
				// calculate probability for 
    prob = 0;
    nj   = 10000;
    for (int j=0; j<nj; j++) {
      y = r->Poisson(mean);
      if (y < Signal) prob++;
    }

    prob = prob/nj;

    prr  = prr+prob;
  }

  prr = prr/ntries;

  printf(" mean = %10.5f probability=%10.6f \n",mean, prr);

  return 0;

}


//-----------------------------------------------------------------------------
// assume that we have an expected background which 'Expected' mean we have
// estimated
// What 'Observed' signal is higher than 'Expected' at 95% confidence level?
// p_expected = 0: p_observed > 0
// p_expected = 1: p_observed > 1
// ... etc ... 
//-----------------------------------------------------------------------------
int poisson_fluct::limit95(double Expected, double Observed) {
  
  // starting value

  //  double mean = Expected;

  double prob;
  int    nj, nk;

  // calculate probability to observe less than signal using MC

  TRandom3* r = new TRandom3();

  double prr = 0;
  
  double y0, y1, pp;

  int ntries = 10;
  nj         = 1000;
  nk         = 1000;

  for (int i=0;i<ntries; i++) {
				// calculate probability for 
    prob = 0;
    for (int j=0; j<nj; j++) {
      y1 = r->Poisson(Expected);

      pp = 0;

      for(int k=0; k<nk; k++) {
	y0 = r->Poisson(Observed);
	if (y0 > y1) pp++;
      }
      pp   = pp/nk;
      prob = prob+pp;
    }

    prob = prob/nj;

    prr  = prr+prob;
  }

  prr = prr/ntries;

  printf(" expected = %10.5f observed=%10.5f prob=%10.6f \n",Expected, Observed, prr);

  return 0;

}


/*
 
root [41] limit95(1.,5.42)
 expected =    1.00000 observed=   5.42000 prob=  0.9502

root [46] limit95(1.6,6.6)
 expected =    1.60000 observed=   6.60000 prob=  0.9507

 6.6/5.4 = 1.22

*/


//-----------------------------------------------------------------------------
// assume that we have a background with the 'Expected' mean 
// we observe 'Observed; events, what is the probability for the background
// to fluctuate that high?
//-----------------------------------------------------------------------------
int poisson_fluct::fluct_prob(double Expected, double Observed) {
  
  // starting value

  //  double mean = Expected;

  //  double prob;
  //  int    nj, nk;

  // calculate probability to observe less than signal using MC

  TRandom3* r = new TRandom3();

  double prr = 0;
  
  double y1, pp;

  pp = 0;
  for (int i=0; i<fNTries; i++) {
				// calculate probability for 
    y1   = r->Poisson(Expected);
    if (y1 > Observed) pp++;
  }

  prr = pp/fNTries;

  printf(" expected = %10.5g observed=%10.5f prob=%12.5e \n",Expected, Observed, prr);

  return 0;

}

