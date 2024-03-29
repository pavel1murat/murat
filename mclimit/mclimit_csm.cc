/* This is a set of routines that works like mclimit.f but uses csm.c to compute
   the chisquared of the data histogram compared to a sum of models, each of which
   may (or may not) have sensitivities to nuisance parameters in their shapes,
   normalizations, and even statistical uncertainty in each bin */

//    25 May 2007 fix MC statistical error problem in Bayesian limit calc.
//    30 May 2007 add a version printout method to the mclimit_csm class
//                add a histogram of ln(1+s/b) in the spirit of plotwithdata, but which
//                combines all channels together.
//    30 July 2007 Work around a problem in Joel's genlimit involving double-precision
//                underflow computing small likelihood functions
//    26 Sep 2007  Add an access method for the MINUIT covariance matrix
//                 put in protection against zero background in a systematic variation in
//                 the Bayesian limit calculators
//     2 Oct 2007 Update the coverage checker to allow for an arbitrary "true" signal rate for
//                 which we'd like to check the coverage.
//     3 Oct 2007 Found another place where the max minuit calls was set to 500 -- raise to 5000
//     5 Oct 2007 Put in access methods for setting max minuit calls, printout switches for pseudoexperiments
//                and whether to call MINOS or not.
//    23 Oct 2007 Put in access methods to control MINUIT's initial step size
//     8 Nov 2007 Fix memory leak -- fitcov not deleted after csm is deleted.  Other fixes from valgrind output --
//                (output should not change). thanks to Kevin Lannon!
//    28 Nov 2007 Scale horizontally interpolated histograms like the vertically interpolated ones.
//     8 Dec 2007 Add a new MC statistical error mode which pays attention to the error bars supplied with
//                each template.  Treatment is approximate -- meant to ease the cases where MC contributions
//                come from a sample of differently-weighted events.
//     9 Dec 2007 Minor change to speed up the new histogram interpolation with errors.
//                Avoid cloning TH1's as that's a very slow operation.
//    11 Dec 2007 Change to formatting of pseudoexperiment printout in the CLs calculation.  MINOS
//                has some unavoidable printf statements which can corrupt the printout, so label
//                the pseudoexperiment data and make sure they both go on the same line.
//    12 Dec 2007 Small changes to the error in each bin handling -- bugfix to option=2
//    18 Dec 2007 Another bugfix to option=2 -- swapping contents and errors was buggy if two of the
//                histogram pointers were identical
//    24 Dec 2007 Speedup -- make a histogram interpolator that does not interpolate
//                errors for use with options 0 and 1 for poissflag
//     7 Feb 2008 Have a separate parameter controlling MINOS's maxcalls
//    11 Feb 2008 Always call IMPROVE after MINIMIZE in the fits
//    12 Feb 2008 Protect against divide by zero in MC stat fitter if a channel's scale factor is zero
//    14 Feb 2008 Always clear bayes_posterior_points when clearing the bayes_posterior vector
//    14 Feb 2008 Add on Marginalized posterior methods for measuring cross sections bh_xsfit and
//                bh_xsfit_withexpect
//    17 Feb 2008 Minor bugfix in bh_xsfit and bh_xsfit_withexpect -- make sure not to drop the
//                largest point in the marginalized posterior when integrating it!
//    26 Feb 2008 Update the interpolator style option to allow or disallow shape extrapolations
//    27 Feb 2008 Add a feature to put bounds on nuisance parameters
//    27 Feb 2008 Fix two bugs in the 2D interpolation with error bars -- one cleared the histogram contents
//                if all errors were zero, and the other set the wrong contents.
//    28 Feb 2008 Add in gclsexpb* routines -- Tom Wright's analysis with unconstrianed fits makes for
//                a case where CLs is not a monotonic function of -2lnQ, so need to scan over all
//                possible outcomes.
//     7 Mar 2008 Add in a quicker Bayesian limit integrator
//                with a cutoff in finite steps.  Sometimes bayes_heinrich and bayes_heinrich_withexpect
//                take too long, particularly for channels with lots of candidates.
//    26 Mar 2008 Fix a minor bug in the gclsexp* routines -- it didn't deal with discrete outcomes
//                as well as it could.
//    21 Apr 2008 Make gclsexp* the default in s95aux.
//    21 Apr 2008 Add a 2D cspdf2d and a function to print it out.  No 2D limits yet -- just print them out
//                for now for further plotting and analysis.  Also wrap a pseudoexperiment loop
//                around it.  New functions:  bh_2d_scan and bh_2d_scan_expect
//    22 Apr 2008 Add to bh_2d_scan_expect the same cross section fits and integrals as in bh_2d_scan.  Also
//                add to the argument list the input cross section scale factors so that coverages can be
//                computed.
//    23 Apr 2008 Add a printout of nuisance parameters for extreme null-hypothesis pseudoexperiments
//     8 May 2008 Fix bugs in the nuisance parameter constraint equations.
//    15 May 2008 2D cross section fit bugfix -- wrong low bound chosen.
//    18 Jun 2008 Add a printout of all bins' s, b, d
//    18 Sep 2008 Add some accessors to get a hold of normalization factors applied to the Bayesian posteriors.
//                Useful when parallelizing a big computation in order to add up many posteriors
//    14 Jan 2009 Protect against a segfault in gclsaux

#define MCLIMIT_CSM_VERSION_NUMBER 3.47
#define MCLIMIT_CSM_VERSION_DATE "Jan 14, 2009"

// Author:  Tom Junk, Fermilab.  trj@fnal.gov

#include "mclimit/mclimit_csm.h"

#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <assert.h>
#include <stddef.h>
#include <algorithm>
#include "TRandom.h"
#include "TMinuit.h"
#include "THStack.h"
#include "TLegend.h"
#include "TList.h"
#include "TMath.h"
#include "TString.h"

using namespace std;

#define max(max_a,max_b) (((max_a)>(max_b)) ? (max_a) : (max_b))
#define min(min_a,min_b) (((min_a)>(min_b)) ? (min_b) : (min_a))

// Minuit ugliness -- data communication with the function to fit is either
// via member functions of a new class inherited from TObject (not done here),
// or in global data (at least global to this source file) storage, which is
// the method chosen here because it's easier.

double csint0(double xlo,double logscale,
		     int nchan,const int nobs[],const EB chan[],
		     int ngl,const double xgl[],const double lwgl[],
		     PRIOR prior);

void csint02cut(double xlo1,double xlo2,double xhi,double logscale,
	          int nchan,const int nobs[],const EB chan[],
		  int ngl,const double xgl[],const double lwgl[],PRIOR prior,
	          double* int1,double* int2);

void gameansigma(double *mean,double *sigma,
			int nchan,int nens,const int nobs[],const EB* ens);
double arcfreq(double y);
#define freq(x) (0.5*erfc(-0.707106781186547524*(x)))

void csint02(double xlo1,double xlo2,double logscale,
		    int nchan,const int nobs[],const EB chan[],
		    int ngl,const double xgl[],const double lwgl[],PRIOR prior,
		    double* int1,double* int2);

// some globals which really should be put into classes, but the Bayesian routines are
// written in C and not C++

double dlcsn=0;
double dlcsn2d=0;
double bhnorm=0;

/*----------------------------------------------------------------------------*/

// constructor

mclimit_csm::mclimit_csm()
{
  nmc = 0;
  nmc_req = 10000;
  recalctsflag = 1;
  // set null pointers to our cumulative histograms -- if we
  // don't get any from the user, don't bother filling them.
  nullnullchisquare = 0;
  nulltestchisquare = 0;
  testnullchisquare = 0;
  testtestchisquare = 0;
  // null pointers to the test statistic arrays -- need to allocate memory when
  // we know how many to do
  tss = 0;
  tsb = 0;
  wtss = 0;
  wtsb = 0;
  itss = 0;
  itsb = 0;

  bayes_interval_begin = 0;
  bayes_interval_end = 0;
  bayes_interval_step = 0;
  bayes_posterior.clear();
  bayes_posterior_points.clear();
  bayes_pseudoexperiment_limits = 0;

  minuitmaxcalls = 500;
  minosmaxcalls = 500;
  minuitstepsize = 0.1;
  minuitprintflag = 0;
  minosflag = 0;
  pxprintflag = 0;
  bayesintegralmethod = CSM_BAYESINTEGRAL_JOEL;

  extremenpprintflag = 0;
  extremenpprintvalue = 0;

  fBayesRandom = new TRandom3();
}

/*----------------------------------------------------------------------------*/

// destructor

mclimit_csm::~mclimit_csm()
{
  Int_t i;

  // deallocate cloned input data histograms and their names.

  for (i=0; i<(Int_t) datahist.size(); i++)
    {
      delete datahist[i];
      delete[] dhname[i]; 
    }

  delete fBayesRandom;
}

/*----------------------------------------------------------------------------*/

void mclimit_csm::print_version()
{
  cout << "Version information for mclimit_csm.C: " << endl;
  cout << "Version number: " << MCLIMIT_CSM_VERSION_NUMBER << endl;
  cout << "Version date:   " << MCLIMIT_CSM_VERSION_DATE << endl;
}

void mclimit_csm::setminuitmaxcalls(Int_t maxcalls)
{
  minuitmaxcalls = maxcalls;
}
Int_t mclimit_csm::getminuitmaxcalls()
{
  return(minuitmaxcalls);
}

void mclimit_csm::setminosmaxcalls(Int_t maxcalls)
{
  minosmaxcalls = maxcalls;
}
Int_t mclimit_csm::getminosmaxcalls()
{
  return(minosmaxcalls);
}

void mclimit_csm::setminuitstepsize(Double_t stepsize)
{
  minuitstepsize = stepsize;
}
Double_t mclimit_csm::getminuitstepsize()
{
  return(minuitstepsize);
}

void mclimit_csm::setprintflag(bool pf)
{
  minuitprintflag = pf;
}
bool mclimit_csm::getprintflag()
{
  return(minuitprintflag);
}

void mclimit_csm::setminosflag(bool mf)
{
  minosflag = mf;
}
bool mclimit_csm::getminosflag()
{
  return(minosflag);
}

void mclimit_csm::setpxprintflag(bool pf)
{
  pxprintflag = pf;
}
bool mclimit_csm::getpxprintflag()
{
  return(pxprintflag);
}

void mclimit_csm::set_bayes_integration_method(int imethod)
{
  bayesintegralmethod = imethod;
}

int mclimit_csm::get_bayes_integration_method()
{
  return(bayesintegralmethod);
}

/*----------------------------------------------------------------------------*/
/* Build the list of channel data histograms and channel names that is sorted by channel */
/* name during the building process.  Use the same sorting procedure as used in */
/* csm_model::lookup_add_channame so that our data histograms are stored in */
/* the same order as our models */

void mclimit_csm::set_datahist(TH1 *h, const char *cname)
{
  Int_t i,ifound,j,jfound;
  char *s;
  vector<char*>::iterator nhi;
  vector<TH1*>::iterator dhi;

  recalctsflag = 1;

  ifound = -1;
  jfound = -1;
  for (i=0; i < (Int_t) dhname.size(); i++)
    {
      j = (Int_t) strcmp(cname,dhname[i]);
      if (j == 0)
	{
	  ifound = i;
	}
      if (j>0 && jfound == -1)
	{
	  jfound = i;
	}
    }
  /* if the name isn't already in the list, add it to the vector of names and
     make a blank model for it too.  Put the new name in it sorted place, sorted
     by increasing sort order of the name strings.  If the name is on the 
     list, replace the existing data histogram with a clone of the one supplied. */

  if (ifound == -1)
    {
      s = new char[strlen(cname)+1];
      strcpy(s,cname);
      if (jfound == -1)
	{
          dhname.push_back(s);
          datahist.push_back((TH1*) h->Clone());
	}
      else
	{
	  nhi = dhname.begin() + jfound;
	  dhname.insert(nhi,s);
	  dhi = datahist.begin() + jfound;
	  datahist.insert(dhi,(TH1*) h->Clone());
	}
    }
  else
    {
      delete datahist[ifound];
      datahist[ifound] = (TH1*) h->Clone();
    }
}

/*----------------------------------------------------------------------------*/

void mclimit_csm::set_npe(Int_t nperequest)
{
  if (nperequest < 0)
    {
      cout << "mclimit_csm::set_npe: Invalid pseudoexperiment request: " << nperequest << endl;
      exit(0);
    }
  nmc_req = nperequest;
}

/*----------------------------------------------------------------------------*/

Int_t mclimit_csm::get_npe()
{
  return(nmc_req);
}

/*----------------------------------------------------------------------------*/

void mclimit_csm::set_null_hypothesis(csm_model *model)
{
  null_hypothesis = model;
  recalctsflag = 1;
  nmc = 0;
}

/*----------------------------------------------------------------------------*/

void mclimit_csm::set_test_hypothesis(csm_model *model)
{
  test_hypothesis = model;
  recalctsflag = 1;
  nmc = 0;
}

/*----------------------------------------------------------------------------*/

void mclimit_csm::set_null_hypothesis_pe(csm_model *model)
{
  null_hypothesis_pe = model;
  nmc = 0;
}

/*----------------------------------------------------------------------------*/

void mclimit_csm::set_test_hypothesis_pe(csm_model *model)
{
  test_hypothesis_pe = model;
  nmc = 0;
}

//-----------------------------------------------------------------------------
// test statistic of observed data : ts = chi2(test,data)-chi2(null,data)
//-----------------------------------------------------------------------------
Double_t mclimit_csm::ts()
{
  Int_t i;

  if (recalctsflag) {
      // copy the data histogram pointers into a flat array
      // for the chisquare calculator.
    const TH1** darray = new const TH1*[datahist.size()];
    for (i=0;i<(Int_t)datahist.size();i++) {
      darray[i] = datahist[i];
    }
    tsd          = calc_chi2(test_hypothesis,darray)-calc_chi2(null_hypothesis,darray);
    recalctsflag = 0;
    delete[] darray;
  }
  return tsd;
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::tsbm2()
{
  Int_t i;
  if (nmc==0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  i = (Int_t) nearbyint(nmc*MCLIMIT_CSM_MCLM2S);
  if (i<0) 
    { i=0; }
  if (i>nmc-1) 
    {i=nmc-1;}
 
  return(tsb[itsb[i]]);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::tsbm1()
{
  Int_t i;
  if (nmc==0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  i = (Int_t) nearbyint(nmc*MCLIMIT_CSM_MCLM1S);
  if (i<0) 
    { i=0; }
  if (i>nmc-1) 
    {i=nmc-1;}
 
  return(tsb[itsb[i]]);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::tsbmed()
{
  Int_t i;
  if (nmc==0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  i = (Int_t) nearbyint(nmc*MCLIMIT_CSM_MCLMED);
  if (i<0) 
    { i=0; }
  if (i>nmc-1) 
    {i=nmc-1;}
 
  return(tsb[itsb[i]]);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::tsbp1()
{
  Int_t i;
  if (nmc==0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  i = (Int_t) nearbyint(nmc*MCLIMIT_CSM_MCLP1S);
  if (i<0) 
    { i=0; }
  if (i>nmc-1) 
    {i=nmc-1;}
 
  return(tsb[itsb[i]]);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::tsbp2()
{
  Int_t i;
  if (nmc==0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  i = (Int_t) nearbyint(nmc*MCLIMIT_CSM_MCLP2S);
  if (i<0) 
    { i=0; }
  if (i>nmc-1) 
    {i=nmc-1;}
 
  return(tsb[itsb[i]]);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::tssm2()
{
  Int_t i;
  if (nmc==0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  i = (Int_t) nearbyint(nmc*MCLIMIT_CSM_MCLM2S);
  if (i<0) 
    { i=0; }
  if (i>nmc-1) 
    {i=nmc-1;}
 
  return(tss[itss[i]]);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::tssm1()
{
  Int_t i;
  if (nmc==0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  i = (Int_t) nearbyint(nmc*MCLIMIT_CSM_MCLM1S);
  if (i<0) 
    { i=0; }
  if (i>nmc-1) 
    {i=nmc-1;}
 
  return(tss[itss[i]]);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::tssmed()
{
  Int_t i;
  if (nmc==0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  i = (Int_t) nearbyint(nmc*MCLIMIT_CSM_MCLMED);
  if (i<0) 
    { i=0; }
  if (i>nmc-1) 
    {i=nmc-1;}
 
  return(tss[itss[i]]);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::tssp1()
{
  Int_t i;
  if (nmc==0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  i = (Int_t) nearbyint(nmc*MCLIMIT_CSM_MCLP1S);
  if (i<0) 
    { i=0; }
  if (i>nmc-1) 
    {i=nmc-1;}
 
  return(tss[itss[i]]);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::tssp2()
{
  Int_t i;
  if (nmc==0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  i = (Int_t) nearbyint(nmc*MCLIMIT_CSM_MCLP2S);
  if (i<0) 
    { i=0; }
  if (i>nmc-1) 
    {i=nmc-1;}
 
  return(tss[itss[i]]);
}

/*----------------------------------------------------------------------------*/
// confidence levels for an arbitrary test statistic.  Used for computing
// actual and expected confidence levels.
// may need to work on this a bit because of finite MC statistics -- it's hard to
// compute the confidence levels on the tails unless the histograms are properly
// filled out there.  Particularly values of clb very close to 1 need to be computed
// with extreme care.
// This is addressed with the reweighting procedure suggested by Alex Read --
// the likelihood ratio can be used to reweight test hypothesis pseudoexperiments
// to model the null hypotheis background distribution, and vice versa.
//-----------------------------------------------------------------------------
Double_t mclimit_csm::clsbaux(Double_t tsaux) {
  Int_t i;
  Double_t clsbloc;
  if (nmc == 0) {
    cout << "Need to run pseudoexperiments after defining/changing models" << endl;
    cout << "and before calling results routines -- mclimit_csm" << endl;
    return(0);
  }
  clsbloc = 0;
  for (i=0;i<nmc;i++) {
    if (tss[itss[i]] < tsaux) { 
      clsbloc = (i+1.)/((Double_t) nmc);
    }
  }
  return(1-clsbloc);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsbauxw(Double_t tsaux)
{
  Int_t i;
  Double_t clsbloc;
  if (nmc == 0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  clsbloc = 0;
  for (i=0;i<nmc;i++)
    {
      if (tsb[itsb[i]] >= tsaux)
	{ 
          clsbloc += wtsb[itsb[i]];
	}
    }
  clsbloc /= ((Double_t) nmc);
  return(clsbloc);
}

//-----------------------------------------------------------------------------
Double_t mclimit_csm::clbaux(Double_t tsaux) {

  Int_t    i;
  Double_t clbloc;

  if (nmc==0) {
    cout << "Need to run pseudoexperiments after defining/changing models" << endl;
    cout << "and before calling results routines -- mclimit_csm" << endl;
    return 0;
  }

  clbloc = 0;
  for (i=0;i<nmc;i++) {
    if (tsb[itsb[i]] < tsaux) { 
      clbloc = (i+1.)/((Double_t) nmc);
    }
  }
  return (1-clbloc);
}

/*----------------------------------------------------------------------------*/

// compute 1-CLb (as a p-value, including the outcomes with exactly the -2lnQ
// observed) using null hypothesis px's.

Double_t mclimit_csm::omclbaux(Double_t tsaux)
{
  Int_t i;
  Double_t omclbloc;
  if (nmc==0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  omclbloc = 1;
  for (i=0;i<nmc;i++)
    {
      if (tsb[itsb[i]] <= tsaux)
	{ 
	  omclbloc = (i+1.)/((Double_t) nmc);
	}
    }
  return(omclbloc);
}

/*----------------------------------------------------------------------------*/

/*  Compute 1-CLb using reweighted test hypothesis pseudoexperiments.  Reweight
    using the inverese of the likelihood ratio, 1/(p(data|test)/p(data|null)) */

Double_t mclimit_csm::omclbauxw(Double_t tsaux)
{
  Int_t i;
  Double_t omclbloc = 0;

  if (nmc==0)
    {
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  omclbloc = 0;
  for (i=0;i<nmc;i++)
    {
      if (tss[itss[i]] <= tsaux)
	{ 
          omclbloc += wtss[itss[i]];
	}
    }
  omclbloc /= ((Double_t) nmc);
  return(omclbloc);
}

/*----------------------------------------------------------------------------*/
Double_t mclimit_csm::clsaux(Double_t tsaux) {
  Double_t clbloc, clsloc;
  clbloc = clbaux(tsaux);
  if (clbloc > 0) {
    clsloc = clsbaux(tsaux)/clbloc;
  }
  else {
    clsloc = 1;
  }
  return clsloc;
}

/*----------------------------------------------------------------------------*/
Double_t mclimit_csm::clsauxw(Double_t tsaux) {

  Double_t clbloc,clsloc;
  clbloc = clbaux(tsaux);
  if (clbloc > 0)
    {
      clsloc = clsbauxw(tsaux)/clbloc;
    }
  else
    {
      clsloc = 1;
    }
  return(clsloc);
}

//-----------------------------------------------------------------------------
// confidence levels using the data test statistic.  Recompute the data test
// statistic if need be.
//-----------------------------------------------------------------------------
Double_t mclimit_csm::cls() {
  return(clsaux(ts()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsb()
{
  return(clsbaux(ts()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsw()
{
  return(clsauxw(ts()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsbw()
{
  return(clsbauxw(ts()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clb()
{
  return(clbaux(ts()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclb()
{
  return(omclbaux(ts()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbw()
{
  return(omclbauxw(ts()));
}

/*----------------------------------------------------------------------------*/
// keep the convention of the other clsexp routines

Double_t mclimit_csm::gclsexpbm2()  { return(gclsaux(MCLIMIT_CSM_MCLP2S)); }
Double_t mclimit_csm::gclsexpbm1()  { return(gclsaux(MCLIMIT_CSM_MCLP1S)); }
Double_t mclimit_csm::gclsexpbmed() { return(gclsaux(MCLIMIT_CSM_MCLMED)); }
Double_t mclimit_csm::gclsexpbp1()  { return(gclsaux(MCLIMIT_CSM_MCLM1S)); }
Double_t mclimit_csm::gclsexpbp2()  { return(gclsaux(MCLIMIT_CSM_MCLM2S)); }

//-----------------------------------------------------------------------------
Double_t mclimit_csm::gclsaux(Double_t thresh) {

  vector<Double_t> clslist;
  double           tstmp;
  if (nmc == 0) {
    cout << "Need to run pseudoexperiments after defining/changing models" << endl;
    cout << "and before calling results routines -- mclimit_csm" << endl;
    return 0;
  }
  // get CLs for each background outcome, sort them (they could be out of order as a function of -2lnQ)
  // tss and tsb are already sorted in ascending order, using the sort indices
  // fix up by counting duplicates properly

  int k=0;
  int j=0;
  for (int i=0; i<nmc; i++) {
    tstmp = tsb[itsb[i]];
    while ((tsb[itsb[k+1]] < tstmp) && (k < nmc-2)) k++;
    while ((tss[itss[j+1]] < tstmp) && (j < nmc-2)) j++;
    clslist.push_back(  ((Double_t) (nmc-j)) / ((Double_t) (nmc-k)) );
  }
  std::sort(clslist.begin(),clslist.end());
  int i =  (int) nearbyint( ((Double_t) nmc)*thresh);
  return(clslist[i]);
 }

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpbm2()
{
  return(clsaux(tsbm2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpbm1()
{
  return(clsaux(tsbm1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpbmed()
{
  return(clsaux(tsbmed()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpbp1()
{
  return(clsaux(tsbp1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpbp2()
{
  return(clsaux(tsbp2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpbm2w()
{
  return(clsauxw(tsbm2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpbm1w()
{
  return(clsauxw(tsbm1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpbmedw()
{
  return(clsauxw(tsbmed()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpbp1w()
{
  return(clsauxw(tsbp1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpbp2w()
{
  return(clsauxw(tsbp2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsbexpbm2()
{
  return(clsbaux(tsbm2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsbexpbm1()
{
  return(clsbaux(tsbm1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsbexpbmed()
{
  return(clsbaux(tsbmed()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsbexpbp1()
{
  return(clsbaux(tsbp1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsbexpbp2()
{
  return(clsbaux(tsbp2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpsm2()
{
  return(clsaux(tssm2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpsm1()
{
  return(clsaux(tssm1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpsmed()
{
  return(clsaux(tssmed()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpsp1()
{
  return(clsaux(tssp1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clsexpsp2()
{
  return(clsaux(tssp2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clbexpbm2()
{
  return(clbaux(tsbm2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clbexpbm1()
{
  return(clbaux(tsbm1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clbexpbmed()
{
  return(clbaux(tsbmed()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clbexpbp1()
{
  return(clbaux(tsbp1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clbexpbp2()
{
  return(clbaux(tsbp2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clbexpsm2()
{
  return(clbaux(tssm2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clbexpsm1()
{
  return(clbaux(tssm1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clbexpsmed()
{
  return(clbaux(tssmed()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clbexpsp1()
{
  return(clbaux(tssp1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::clbexpsp2()
{
  return(clbaux(tssp2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpbm2()
{
  return(omclbaux(tsbm2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpbm1()
{
  return(omclbaux(tsbm1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpbmed()
{
  return(omclbaux(tsbmed()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpbp1()
{
  return(omclbaux(tsbp1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpbp2()
{
  return(omclbaux(tsbp2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpsm2()
{
  return(omclbaux(tssm2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpsm1()
{
  return(omclbaux(tssm1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpsmed()
{
  return(omclbaux(tssmed()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpsp1()
{
  return(omclbaux(tssp1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpsp2()
{
  return(omclbaux(tssp2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpbm2w()
{
  return(omclbauxw(tsbm2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpbm1w()
{
  return(omclbauxw(tsbm1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpbmedw()
{
  return(omclbauxw(tsbmed()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpbp1w()
{
  return(omclbauxw(tsbp1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpbp2w()
{
  return(omclbauxw(tsbp2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpsm2w()
{
  return(omclbauxw(tssm2()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpsm1w()
{
  return(omclbauxw(tssm1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpsmedw()
{
  return(omclbauxw(tssmed()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpsp1w()
{
  return(omclbauxw(tssp1()));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::omclbexpsp2w()
{
  return(omclbauxw(tssp2()));
}



/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::p2sigmat()
{
  Int_t i;
  Double_t p2s;
  p2s = 0;

  if (nmc == 0)
    {
      cout << "mclimit_csm::p2sigmat: " << endl;
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  for (i=0;i<nmc;i++)
    {
      if (omclbauxw(tss[itss[i]]) <= MCLIMIT_CSM_2S)
	{
	  p2s = ((Double_t) i)/((Double_t) nmc);
	}
    }
  return(p2s);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::p3sigmat()
{
  Int_t i;
  Double_t p3s;
  p3s = 0;

  if (nmc == 0)
    {
      cout << "mclimit_csm::p3sigmat: " << endl;
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  for (i=0;i<nmc;i++)
    {
      if (omclbauxw(tss[itss[i]]) <= MCLIMIT_CSM_3S)
	{
	  p3s = ((Double_t) i)/((Double_t) nmc);
	}
    }
  return(p3s);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::p5sigmat()
{
  Int_t i;
  Double_t p5s;
  p5s = 0;

  if (nmc == 0)
    {
      cout << "mclimit_csm::p5sigmat: " << endl;
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  for (i=0;i<nmc;i++)
    {
      if (omclbauxw(tss[itss[i]]) <= MCLIMIT_CSM_5S)
	{
	  p5s = ((Double_t) i)/((Double_t) nmc);
	}
    }
  return(p5s);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::p2sigman()
{
  Int_t i;
  Double_t p2s;
  p2s = 0;

  if (nmc == 0)
    {
      cout << "mclimit_csm::p2sigman: " << endl;
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  for (i=0;i<nmc;i++)
    {
      if (omclbaux(tsb[itsb[i]]) <= MCLIMIT_CSM_2S)
	{
	  p2s = ((Double_t) i)/((Double_t) nmc);
	}
    }
  return(p2s);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::p3sigman()
{
  Int_t i;
  Double_t p3s;
  p3s = 0;

  if (nmc == 0)
    {
      cout << "mclimit_csm::p3sigman: " << endl;
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  for (i=0;i<nmc;i++)
    {
      if (omclbaux(tsb[itsb[i]]) <= MCLIMIT_CSM_3S)
	{
	  p3s = ((Double_t) i)/((Double_t) nmc);
	}
    }
  return(p3s);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::p5sigman()
{
  Int_t i;
  Double_t p5s;
  p5s = 0;

  if (nmc == 0)
    {
      cout << "mclimit_csm::p5sigman: " << endl;
      cout << "Need to run pseudoexperiments after defining/changing models" << endl;
      cout << "and before calling results routines -- mclimit_csm" << endl;
      return(0);
    }
  for (i=0;i<nmc;i++)
    {
      if (omclbaux(tsb[itsb[i]]) <= MCLIMIT_CSM_5S)
	{
	  p5s = ((Double_t) i)/((Double_t) nmc);
	}
    }
  return(p5s);
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::s95()
{
  return(s95aux(MCLIMIT_CSM_CLS));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::s95m2()
{
  return(s95aux(MCLIMIT_CSM_CLSM2));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::s95m1()
{
  return(s95aux(MCLIMIT_CSM_CLSM1));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::s95med()
{
  return(s95aux(MCLIMIT_CSM_CLSMED));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::s95p1()
{
  return(s95aux(MCLIMIT_CSM_CLSP1));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::s95p2()
{
  return(s95aux(MCLIMIT_CSM_CLSP2));
}


//-----------------------------------------------------------------------------
Double_t mclimit_csm::s95aux(Int_t itype) {

  Double_t   sf,sfl,sfh,cltest,cll,cla,clh;
  Int_t      j,foundit;
  csm_model  *testhypsave , *testhyppesave;
  csm_model  *scaledsignal, *scaledsignalpe;

  sfl    = 0;
  sfh    = 0;
  cltest = 0;

  // this hypothesis never gets destroyed, but we lose the pointer to it

  testhypsave   = test_hypothesis;
  testhyppesave = test_hypothesis_pe;

  sf      = 1.0;
  cll     = 0.0;
  cla     = 0.0;
  clh     = 0.0; 
  foundit = 0;

  for (j=0;j<32;j++) {
    //cout << "in s95aux seek " << j << " scale: " << sf << endl;
    if (foundit == 0) {
      scaledsignal       = testhypsave->scalesignal  (sf);
      scaledsignalpe     = testhyppesave->scalesignal(sf);
      test_hypothesis    = scaledsignal;
      test_hypothesis_pe = scaledsignalpe;
      run_pseudoexperiments();
      recalctsflag = 1;
      if (itype == MCLIMIT_CSM_CLS) {
					// s95 gets here
	cltest = cls(); 
      }
      else if (itype == MCLIMIT_CSM_CLSM2) {
	cltest = gclsexpbm2();
      }
      else if (itype == MCLIMIT_CSM_CLSM1) {
	cltest = gclsexpbm1();
      }
      else if (itype == MCLIMIT_CSM_CLSMED)  {
	cltest = gclsexpbmed();
      }
      else if (itype == MCLIMIT_CSM_CLSP1)  {
	cltest = gclsexpbp1();
      }
      else if (itype == MCLIMIT_CSM_CLSP2) {
	cltest = gclsexpbp2();
      }

      delete scaledsignal;
      delete scaledsignalpe;

      if (j==0){
	cla = cltest;
      }
      if (cltest<0.05) {
	if (cla>0.05) {
	  sfh = sf;
	  clh = cltest;
	  sfl = sf/2.0;
	  cll = cla;
	  foundit = 1;
	}
	sf /= 2.0;
      }
      else if (cltest>0.05){
	if (cla<0.05) {
	  sfl = sf;
	  cll = cltest;
	  sfh = sf*2.0;
	  clh = cla;
	  foundit = 1;
	}
	sf *= 2.0;
      }
      else {
	test_hypothesis    = testhypsave;
	test_hypothesis_pe = testhyppesave;
	return(sf);
      }
      cla = cltest;
    }
  } // end of loop over 32 powers of 2 in search of a signal scale factor which brackets
  // 95% CL exclusion

  //cout << "done with seek loop " << sf << endl;
  sf = sfh;
  if (foundit == 0) { 
    cout << "mclimit_csm::s95** could not find s95 within 2**32 of original guess" << endl;
    sf = 0;
  }
  else {
    // From Tom Wright -- speed up by doing a deterministic five more
    // calcs of CL and a linear fit of log(CL) vs. sf.

    // find error on 0.05 CL for number of PEs
    double dcl=sqrt(0.05*0.95/nmc);

    // put in some protection against logarithms of negative numbers
    // makes sure -5*dcl + 0.05 is not negative.

    dcl = min(dcl,0.0099);
//-----------------------------------------------------------------------------
// try +6sigma, +3sigma, 0sigma, -3sigma, -6sigma regions
// increment stuff used for linear fit of ln(CL) vs sf
//-----------------------------------------------------------------------------
    double lf_a=0, lf_b=0, lf_c=0, lf_d=0, lf_e=0, lf_f=0;
    for( int j=-5; j<6; j+=2 ) {
      sf = sfl + (log(0.05+j*dcl) - log(cll))*
	(sfl-sfh)/(log(cll)-log(clh));

      double clsbtest = 0;
      double clbtest  = 0;

      // calculate CL for this sf

      csm_model* scaledsignal   = testhypsave->scalesignal(sf);
      csm_model* scaledsignalpe = testhyppesave->scalesignal(sf);

      test_hypothesis    = scaledsignal;
      test_hypothesis_pe = scaledsignalpe;

      run_pseudoexperiments();
      recalctsflag = 1;

      if (itype == MCLIMIT_CSM_CLS) {
	cltest   = cls ();
	clbtest  = clb ();
	clsbtest = clsb();
      }
      else if (itype == MCLIMIT_CSM_CLSM2) {
	cltest = clsexpbm2();
	clbtest = clbexpbm2();
	clsbtest = clsbexpbm2();
      }
      else if (itype == MCLIMIT_CSM_CLSM1) {
	cltest = clsexpbm1();
	clbtest = clbexpbm1();
	clsbtest = clsbexpbm1();
      }
      else if (itype == MCLIMIT_CSM_CLSMED) {
	cltest = clsexpbmed();
	clbtest = clbexpbmed();
	clsbtest = clsbexpbmed();
      }
      else if (itype == MCLIMIT_CSM_CLSP1) {
	cltest = clsexpbp1();
	clbtest = clbexpbp1();
	clsbtest = clsbexpbp1();
      }
      else if (itype == MCLIMIT_CSM_CLSP2) {
	cltest = clsexpbp2();
	clbtest = clbexpbp2();
	clsbtest = clsbexpbp2();
      }
      
      delete scaledsignal;
      delete scaledsignalpe;
      
      // double dcltest=sqrt(cltest*(1-cltest)/nmc);
      // 	  double dcltest = cltest*sqrt((1-clbtest)/clbtest/nmc +
      // 				       (1-clsbtest)/clsbtest/nmc);
      double dcltest=sqrt(clsbtest*(1-clsbtest)/nmc)/clbtest;

      //  printf("%f %f %f %f %f\n",sf,clbtest,clsbtest,cltest,dcltest);

      double lcl = log(cltest);
      double dlcl = dcltest/cltest;

      lf_a += sf/dlcl/dlcl;
      lf_b += 1/dlcl/dlcl;
      lf_c += lcl/dlcl/dlcl;
      lf_d += sf*sf/dlcl/dlcl;
      lf_e += sf*lcl/dlcl/dlcl;
      lf_f += lcl*lcl/dlcl/dlcl;
    }

    // Find fit parameters for log(CL)=p1+p2*sf
    double lf_p1 = (lf_d*lf_c-lf_e*lf_a)/(lf_d*lf_b-lf_a*lf_a);
    double lf_p2 = (lf_e*lf_b-lf_c*lf_a)/(lf_d*lf_b-lf_a*lf_a);

    //double lf_dp1 = sqrt(lf_d/(lf_b*lf_d-lf_a*lf_a));
    //double lf_dp2 = sqrt(lf_b/(lf_b*lf_d-lf_a*lf_a));
    //double lf_rho = -lf_a/(lf_b*lf_d-lf_a*lf_a)/lf_dp1/lf_dp2;
    
    //printf("fit results %f %f %f %f %f\n",lf_p1,lf_dp1,lf_p2,lf_dp2,lf_rho);

    //double lf_x2 = lf_f-2*lf_p2*lf_e-2*lf_p1*lf_c+lf_p2*lf_p2*lf_d+
    //	2*lf_p1*lf_p2*lf_a+lf_p1*lf_p1*lf_b;
    //printf("chisuare/dof: %f\n",lf_x2/4);
    
    // invert to get sf at 0.05 and its error
    // assuming 100% anticorrelation
    double lf_sf = (log(0.05)-lf_p1)/lf_p2;
    //printf("CL variation at %f: %f %f %f\n",lf_sf,
    //     exp(lf_p1-lf_dp1+(lf_p2+lf_dp2)*lf_sf),
    //     exp(lf_p1+lf_p2*lf_sf),
    //     exp(lf_p1+lf_dp1+(lf_p2-lf_dp2)*lf_sf));
    
    //double lf_dsf1 = (log(0.05)-lf_p1-lf_dp1)/(lf_p2-lf_dp2);
    //double lf_dsf2 = (log(0.05)-lf_p1+lf_dp1)/(lf_p2+lf_dp2);
    
    //printf("SF variation: %f %f %f\n",lf_dsf1,lf_sf,lf_dsf2);
    sf = lf_sf;
  }

  test_hypothesis    = testhypsave;
  test_hypothesis_pe = testhyppesave;
  recalctsflag       = 1;
  return sf;
}


/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::lumi95()
{
  return(lumipaux(MCLIMIT_CSM_LUMI95));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::lumi3s()
{
  return(lumipaux(MCLIMIT_CSM_LUMI3S));
}

/*----------------------------------------------------------------------------*/

Double_t mclimit_csm::lumi5s()
{
  return(lumipaux(MCLIMIT_CSM_LUMI5S));
}

/*----------------------------------------------------------------------------*/
// compute median amounts of luminosity needed for 95% CL exclusion, 3 sigma
// evidence, or 5 sigma discovery -- scale the systematic errors with 1/sqrt(lumi/lumi_0)

Double_t mclimit_csm::lumipaux(Int_t itype)
{
  Double_t sf,sfl,sfh,cltest,cll,cla,clh;
  Int_t foundit,j;
  csm_model *testhypsave,*nullhypsave;
  csm_model *testhyppesave,*nullhyppesave;
  Double_t resdes;

  resdes = 0.5; // do median luminosity thresholds

  sfl = 0;
  sfh = 0;
  cltest = 0;

  sf = 1.0;
  cll = 0.0;
  cla = 0.0;
  clh = 0.0; 
  foundit = 0;

  testhypsave = test_hypothesis;
  nullhypsave = null_hypothesis;
  testhyppesave = test_hypothesis_pe;
  nullhyppesave = null_hypothesis_pe;
  
  for (j=0;j<32;j++)
    {
      if (foundit == 0)
	{
	  csm_model* scaledtest = testhypsave->scale_err(sf);
	  csm_model* scalednull = nullhypsave->scale_err(sf);
	  csm_model* scaledtestpe = testhyppesave->scale_err(sf);
	  csm_model* scalednullpe = nullhyppesave->scale_err(sf);
	  test_hypothesis = scaledtest;
	  null_hypothesis = scalednull;
	  test_hypothesis_pe = scaledtestpe;
	  null_hypothesis_pe = scalednullpe;
          run_pseudoexperiments();
          recalctsflag = 1;

          if (itype == MCLIMIT_CSM_LUMI95)
	    {
	      cltest = p2sigmat();
	    }
          else if (itype == MCLIMIT_CSM_LUMI3S)
	    {
	      cltest = p3sigmat();
	    }
          else if (itype == MCLIMIT_CSM_LUMI5S)
	    {
	      cltest = p5sigmat();
	    }
	  delete scaledtest;
	  delete scalednull;
	  delete scaledtestpe;
	  delete scalednullpe;

	  if (j==0)
	    {
	      cla = cltest;
	    }
	  if (cltest < resdes)
	    {
	      if (cla > resdes)
	        {
		  sfl = sf;
	  	  cll = cltest;
		  sfh = sf*2.0;
		  clh = cla;
		  foundit = 1;
	        }
	      sf *= 2.0;
	    }
  	  else if (cltest > resdes)
	    {
	      if (cla < resdes)
	        {
 		  sfh = sf;
		  clh = cltest;
		  sfl = sf/2.0;
		  cll = cla;
		  foundit = 1;
	        }
	      sf /= 2.0;
	    }
	  else
	    {
	      test_hypothesis = testhypsave;
	      null_hypothesis = nullhypsave;
	      test_hypothesis_pe = testhyppesave;
	      null_hypothesis_pe = nullhyppesave;
	      return(sf);
	    }
	  cla = cltest;
        }
    } // end of loop over 32 powers of 2 in search of a luminosity scale factor which
      // brackets the desired sensitvity

  sf = sfl;
  if (foundit == 0)
    { 
      cout << "mclimit_csm::lumipaux** could not find s95 within 2**32 of original guess" << endl;
      sf = 0;
    }
  else
    {

      // From Tom Wright -- speed up by doing a deterministic five more
      // calcs of CL and a linear fit of log(CL) vs. sf.

      // find error on resdes CL for number of PEs
      double dcl=sqrt(resdes*(1.0-resdes)/nmc);

      // put in some protection against logarithms of negative numbers
      // makes sure -5*dcl + resdes is not negative.

      dcl = min(dcl,resdes/5 - 0.0001);

      // try +6sigma, +3sigma, 0sigma, -3sigma, -6sigma regions
      // increment stuff used for linear fit of ln(CL) vs sf
      double lf_a=0, lf_b=0, lf_c=0, lf_d=0, lf_e=0, lf_f=0;
      for( int j=-5; j<6; j+=2 )
	{
	  sf = sfl + (log(resdes+j*dcl) - log(cll))*
	    (sfl-sfh)/(log(cll)-log(clh));

	  //double clsbtest, clbtest;

	  // calculate CL for this sf
	  csm_model* scaledtest = testhypsave->scale_err(sf);
	  csm_model* scalednull = nullhypsave->scale_err(sf);
	  csm_model* scaledtestpe = testhyppesave->scale_err(sf);
	  csm_model* scalednullpe = nullhyppesave->scale_err(sf);
	  test_hypothesis = scaledtest;
	  null_hypothesis = scalednull;
	  test_hypothesis_pe = scaledtestpe;
	  null_hypothesis_pe = scalednullpe;
          run_pseudoexperiments();
          recalctsflag = 1;

          if (itype == MCLIMIT_CSM_LUMI95)
	    {
	      cltest = p2sigmat();
	    }
          else if (itype == MCLIMIT_CSM_LUMI3S)
	    {
	      cltest = p3sigmat();
	    }
          else if (itype == MCLIMIT_CSM_LUMI5S)
	    {
	      cltest = p5sigmat();
	    }

	  delete scaledtest;
	  delete scalednull;
	  delete scaledtestpe;
	  delete scalednullpe;

	  // double dcltest=sqrt(cltest*(1-cltest)/nmc);
// 	  double dcltest = cltest*sqrt((1-clbtest)/clbtest/nmc +
// 				       (1-clsbtest)/clsbtest/nmc);
	  double dcltest=sqrt(cltest*(1-cltest)/nmc)/cltest;

          //  printf("%f %f %f %f %f\n",sf,clbtest,clsbtest,cltest,dcltest);

	  double lcl = log(cltest);
	  double dlcl = dcltest/cltest;

	  lf_a += sf/dlcl/dlcl;
	  lf_b += 1/dlcl/dlcl;
	  lf_c += lcl/dlcl/dlcl;
	  lf_d += sf*sf/dlcl/dlcl;
	  lf_e += sf*lcl/dlcl/dlcl;
	  lf_f += lcl*lcl/dlcl/dlcl;
	}

      // Find fit parameters for log(CL)=p1+p2*sf
      double lf_p1 = (lf_d*lf_c-lf_e*lf_a)/(lf_d*lf_b-lf_a*lf_a);
      double lf_p2 = (lf_e*lf_b-lf_c*lf_a)/(lf_d*lf_b-lf_a*lf_a);

      //double lf_dp1 = sqrt(lf_d/(lf_b*lf_d-lf_a*lf_a));
      //double lf_dp2 = sqrt(lf_b/(lf_b*lf_d-lf_a*lf_a));
      //double lf_rho = -lf_a/(lf_b*lf_d-lf_a*lf_a)/lf_dp1/lf_dp2;

      //printf("fit results %f %f %f %f %f\n",lf_p1,lf_dp1,lf_p2,lf_dp2,lf_rho);

      //double lf_x2 = lf_f-2*lf_p2*lf_e-2*lf_p1*lf_c+lf_p2*lf_p2*lf_d+
      //	2*lf_p1*lf_p2*lf_a+lf_p1*lf_p1*lf_b;
      //printf("chisuare/dof: %f\n",lf_x2/4);

      // invert to get sf at resdes and its error
      // assuming 100% correlation
      double lf_sf = (log(resdes)-lf_p1)/lf_p2;
      //printf("CL variation at %f: %f %f %f\n",lf_sf,
      //     exp(lf_p1-lf_dp1+(lf_p2+lf_dp2)*lf_sf),
      //     exp(lf_p1+lf_p2*lf_sf),
      //     exp(lf_p1+lf_dp1+(lf_p2-lf_dp2)*lf_sf));

      //double lf_dsf1 = (log(resdes)-lf_p1-lf_dp1)/(lf_p2-lf_dp2);
      //double lf_dsf2 = (log(resdes)-lf_p1+lf_dp1)/(lf_p2+lf_dp2);

      //printf("SF variation: %f %f %f\n",lf_dsf1,lf_sf,lf_dsf2);
      sf = lf_sf;

    }
   
  test_hypothesis = testhypsave;
  null_hypothesis = nullhypsave;
  test_hypothesis_pe = testhyppesave;
  null_hypothesis_pe = nullhyppesave;
  return(sf);
}

/*----------------------------------------------------------------------------*/

void mclimit_csm::tshists(TH1* testhypts, TH1* nullhypts)
{
  int i;
  testhypts->Reset();
  nullhypts->Reset();
  for (i=0;i<nmc;i++)
    { 
      testhypts->Fill(tss[i]);
      nullhypts->Fill(tsb[i]);
    }
}

/*----------------------------------------------------------------------------*/
// makes and overlaid plot of ln(1+s/b) in the user's binning
// assumes a canvas and pad are already set up and plot options are set up

void mclimit_csm::plotlnsb(TH1 *mcb_hist, TH1 *mcs_hist, TH1 *data_hist)
{
  Int_t i,ibinx,ibiny,nbinsx,nbinsy,ic,nc;
  Double_t s,b,dtb,gbc,sbln;
  csm_channel_model *cm;
  TH1 *dhp;

  mcb_hist->Reset();
  mcs_hist->Reset();
  data_hist->Reset();

  for (i=0;i<(Int_t)(test_hypothesis_pe->chanmodel.size());i++)
    {
      // dhp is the data histogram we're using to fill in the ln(1+s/b) histo
      // cm is the channel model for this data histogram
      dhp = datahist[i];
      cm = test_hypothesis_pe->chanmodel[i];
      nc = (Int_t) cm->histotemplate.size();
      nbinsx = dhp->GetNbinsX();
      nbinsy = dhp->GetNbinsY();
      for (ibinx=0;ibinx<nbinsx;ibinx++)
	{
	  for (ibiny=0;ibiny<nbinsy;ibiny++)
	    {
	      s = 0;
	      b = 0;
	      if (nbinsy == 1)
		{ dtb = dhp->GetBinContent(ibinx+1); }
	      else
		{ dtb = dhp->GetBinContent(ibinx+1,ibiny+1); }
	      for (ic=0;ic<nc;ic++)
		{
		  if (nbinsy == 1)
		    { gbc = cm->histotemplate_varied[ic]->GetBinContent(ibinx+1); }
		  else
		    { gbc = cm->histotemplate_varied[ic]->GetBinContent(ibinx+1,ibiny+1); }
		  gbc *= cm->sft_varied[ic];
		  if (cm->scaleflag[ic])
		    { s += gbc; }
		  else
		    { b += gbc; }
		}
	      if (b>0 && s>= 0) 
		{ 
		  sbln = log(1.0+s/b);
		  mcs_hist->Fill(sbln,s);
		  mcb_hist->Fill(sbln,b);
		  data_hist->Fill(sbln,dtb);
		}
	    }
	}
    }

  THStack *hs = new THStack("lnsbstack",data_hist->GetTitle());
  mcs_hist->SetFillColor(2);
  mcb_hist->SetFillColor(3);
  hs->Add(mcb_hist);
  hs->Add(mcs_hist);
  data_hist->GetXaxis()->SetTitle("ln(1+s/b)");
  TLegend *slegend = new TLegend(0.7,0.6,0.89,0.89);
  slegend->AddEntry(mcs_hist,"Signal","F");
  slegend->AddEntry(mcb_hist,"Background","F");
  slegend->AddEntry(data_hist,data_hist->GetName(),"P");
  //slegend->SetHeader(data_hist->GetTitle());

  Double_t mcmax,datamax,plotmax;
  mcmax = hs->GetMaximum();
  datamax = data_hist->GetMaximum();
  datamax += sqrt(datamax);
  plotmax = max(datamax,mcmax);
  hs->SetMaximum(plotmax);
  data_hist->SetMaximum(plotmax);
  hs->Draw("HIST");
  data_hist->SetMarkerStyle(20);
  data_hist->SetMarkerColor(kBlack);
  data_hist->DrawCopy("E0SAME");
  slegend->Draw();
}


/*----------------------------------------------------------------------------*/

void mclimit_csm::set_chisquarehistos(TH1 *nn, TH1 *nt, TH1 *tn, TH1 *tt)
{
  nullnullchisquare = nn;
  nulltestchisquare = nt;
  testnullchisquare = tn;
  testtestchisquare = tt;

}

// steer printing of nuisance parameters for extreme pseudoexperiments

void mclimit_csm::setprintnpextremeflag(bool prf)
{
  extremenpprintflag = prf;
}

void mclimit_csm::setprintextremevalue(Double_t pexval)
{
  extremenpprintvalue = pexval;
}

//-----------------------------------------------------------------------------
// run pseudoexperiments, allocating and filling the tss and tsb arrays
// also fill in wtss and wtsb, for use in reweighting pseudoexperiments, using
// the varied nuisance parameters.  Sort tss and tsb, and keep the sort order
// arrays itss and itsb around so the corresponding wtss and wtsb arrays can
// be used with them
//-----------------------------------------------------------------------------
void mclimit_csm::run_pseudoexperiments() {
  Int_t i, nch;
  char  pdname[1000];
  // Double_t tmp;
  Double_t csnull,cstest;

  // make some histograms to store the pseudodata.

  nch = null_hypothesis_pe->channame.size();

  TH1** pdarray = new TH1*[nch];

  for (i=0; i<nch; i++) {
    sprintf(pdname,"%s pseudodata",null_hypothesis_pe->channame[i]);
    pdarray[i] = (TH1*) null_hypothesis_pe->chanmodel[i]->histotemplate[0]->Clone(pdname);
  }

  // allocate memory for test statistic and weight and sort order storage
  if (tss  != 0) delete[] tss;
  if (tsb  != 0) delete[] tsb;
  if (wtss != 0) delete[] wtss;
  if (wtsb != 0) delete[] wtsb;
  if (itss != 0) delete[] itss;
  if (itsb != 0) delete[] itsb;

  tss  = new Double_t[nmc_req];
  tsb  = new Double_t[nmc_req];
  wtss = new Double_t[nmc_req];
  wtsb = new Double_t[nmc_req];
  itss = new Int_t   [nmc_req];
  itsb = new Int_t   [nmc_req];

  for (i=0; i<nmc_req; i++) {
//-----------------------------------------------------------------------------
// generate a single pseudoexperiment for the NULL hypothesis, 
// 'pdarray' - just a data storage, an array of empty histograms filled in 
//             csm_model::single_pseudoexperiment
//-----------------------------------------------------------------------------
    null_hypothesis_pe->single_pseudoexperiment(pdarray);
    wtsb[i] = weightratio(test_hypothesis_pe,null_hypothesis_pe,pdarray);
    csnull  = calc_chi2(null_hypothesis,(const TH1**) pdarray);
    cstest  = calc_chi2(test_hypothesis,(const TH1**) pdarray);
    if (nullnullchisquare != 0) { nullnullchisquare->Fill(csnull); }
    if (nulltestchisquare != 0) { nulltestchisquare->Fill(cstest); }
    tsb[i] = cstest - csnull;
//     if (pxprintflag != 0) {
//       printf("null hyp chi^2: %14.7e test hyp chi^2: %14.7e\n",csnull,cstest);
//     }
    // diagnostic code below to print out all the nuisance parameter values in the case that the null hypothesis
    // has fluctuated into a very signal-like region of -2lnQ.
    if (extremenpprintflag && tsb[i]<extremenpprintvalue) {
      cout << "Extreme value of -2lnQ for a null hypothesis pseudoexperiment: " 
	   << tsb[i] << " < " << extremenpprintvalue << endl;
      null_hypothesis_pe->print_nuisance_params();
    }
//-----------------------------------------------------------------------------
// generate a single pseudoexperiment for the TEST hypothesis 
// note, that 'pdarray' doesn't include signal - its dimension is defined by 
// the NULL hypothesis
//-----------------------------------------------------------------------------
    test_hypothesis_pe->single_pseudoexperiment(pdarray);
    wtss[i] = weightratio(null_hypothesis_pe,test_hypothesis_pe,pdarray);
    csnull = calc_chi2(null_hypothesis,(const TH1**) pdarray);
    cstest = calc_chi2(test_hypothesis,(const TH1**) pdarray);
    if (testnullchisquare != 0) { testnullchisquare->Fill(csnull); }
    if (testtestchisquare != 0) { testtestchisquare->Fill(cstest); }
    tss[i] = cstest - csnull;
    // cout << "null hyp chisquared: " << csnull << " test hyp chisquared: " << cstest << endl;
    
    if (pxprintflag) {
      cout << " Null, Test hyp px: " << tsb[i] << " " << tss[i] << endl;
    }
  }

  TMath::Sort(nmc_req,tss,itss,0);
  TMath::Sort(nmc_req,tsb,itsb,0);

//   for (i=0; i<nch; i++) {
//     delete pdarray[i];
//   }

  delete[] pdarray;
  nmc = nmc_req;
}


/*----------------------------------------------------------------------------*/

// calculate the ratio p(data|nmodel)/p(data|dmodel)

Double_t mclimit_csm::weightratio(csm_model *nmodel, csm_model *dmodel, TH1 *hist[])
{
  int nchans = nmodel->channame.size();
  int ichan;
  csm_channel_model *ncm;
  csm_channel_model *dcm;
  int nbinsx,nbinsy,ibin,jbin,ic;
  double wr,pn,pd;
  int dtb,ncc,dcc;

  wr = 0;
  for (ichan=0;ichan<nchans;ichan++)
    {
      ncm = nmodel->chanmodel[ichan];
      ncc = ncm->histotemplate.size();
      dcm = dmodel->chanmodel[ichan];
      dcc = dcm->histotemplate.size();

      nbinsx = hist[ichan]->GetNbinsX();
      nbinsy = hist[ichan]->GetNbinsY();
      for (ibin=0;ibin<nbinsx;ibin++)
	{
	  for (jbin=0;jbin<nbinsy;jbin++)
	    {
	      dtb = 0;
	      if (nbinsy == 1)
		{
		  dtb = max(0,(int) nearbyint(hist[ichan]->GetBinContent(ibin+1)));
		}
	      else
		{
		  dtb = max(0,(int) nearbyint(hist[ichan]->GetBinContent(ibin+1,jbin+1)));
		}

	      pn = 0; // prediction for numerator model
	      for (ic=0;ic<ncc;ic++)
		{
		  if (nbinsy == 1)
		    {
		      pn += ncm->histotemplate_varied[ic]->GetBinContent(ibin+1)*
                              ncm->sft_varied[ic];
		    }
		  else
		    {
		      pn += ncm->histotemplate_varied[ic]->GetBinContent(ibin+1,jbin+1)*
                              ncm->sft_varied[ic];
		    }
		}
	      pd = 0; // prediction for denominator model
	      for (ic=0;ic<dcc;ic++)
		{
		  if (nbinsy == 1)
		    {
		      pd += dcm->histotemplate_varied[ic]->GetBinContent(ibin+1)*
                              dcm->sft_varied[ic];
		    }
		  else
		    {
		      pd += dcm->histotemplate_varied[ic]->GetBinContent(ibin+1,jbin+1)*
                              dcm->sft_varied[ic];
		    }
		}
	      if (pd > 0)
		{
	          wr += dtb*log(pn/pd) - pn + pd;
		}
	    }
	}
    }
  // about the limit of exponentials
  if (wr>680.) 
    { wr = 680.; }
  wr = exp(wr);
  return(wr);
}

/*----------------------------------------------------------------------------*/

// minimize the chisquared function over the nuisance parameters

Double_t mclimit_csm::calc_chi2(csm_model* Model, const TH1** Hist) {

  Double_t chisquared;
  csm* mycsm = new csm;
  mycsm->setminuitmaxcalls(minuitmaxcalls);
  mycsm->setminosmaxcalls(minosmaxcalls);
  mycsm->setminuitstepsize(minuitstepsize);
  mycsm->setprintflag(minuitprintflag);
  mycsm->setminosflag(minosflag);

  mycsm->set_modeltofit(Model,Hist);

  chisquared = mycsm->chisquared();
  //  cout << "total number of nuisance parameters: " << mycsm->getnparams() << endl;
  //for (int i=0;i<mycsm->getnparams();i++)
  //  {
  //    cout << mycsm->getpname(i) << endl;
  //  }
  delete mycsm;
  return(chisquared);
}

//-----------------------------------------------------------------------------
TH1* mclimit_csm::get_datahist(Int_t i) {
  return(datahist[i]);
}


//-----------------------------------------------------------------------------
// 2010-01-27 P.Murat
// return pointer to a data histogram corresponding to a 'ChannelName' or NULL
//-----------------------------------------------------------------------------
TH1* mclimit_csm::get_datahist(const char* ChannelName) {
  TH1* h1(0);
  int nch = datahist.size();
  for (int i=0; i<nch; i++) {
    if (strcmp(ChannelName,dhname[i]) == 0) {
      h1 = datahist[i];
      break;
    }
  }

  return h1;
}



//-----------------------------------------------------------------------------
// Print out signal, background and data for all bins in all channels.
// use test_hypothesis_pe which has all components.
//-----------------------------------------------------------------------------
void mclimit_csm::printsbd() {
  Int_t i,j,k;
  Int_t nbinsx,nbinsy;
  
  Int_t nchans = (Int_t) test_hypothesis_pe->channame.size();
  
  for (i=0;i<nchans;i++) {
    nbinsx = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsX();
    nbinsy = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsY();
    for (j=0;j<nbinsx;j++) {
      for (k=0;k<nbinsy;k++) {
	Int_t nobs;
	Double_t nsig = 0;
	Double_t nbkg = 0;
	if (nbinsy==1) { 
	  nobs = (Int_t) nearbyint(datahist[i]->GetBinContent(j+1)); 
	}
	else { 
	  nobs = (Int_t) nearbyint(datahist[i]->GetBinContent(j+1,k+1)); 
	}
	csm_channel_model *cm = test_hypothesis_pe->chanmodel[i];
	Int_t ntemplates = (Int_t) cm->histotemplate.size();
	for(Int_t itpl=0;itpl<ntemplates;itpl++) {
	  Double_t r;
	  if (nbinsy==1) { 
	    r = cm->histotemplate_varied[itpl]->GetBinContent(j+1); 
	  }
	  else { 
	    r = cm->histotemplate_varied[itpl]->GetBinContent(j+1,k+1); 
	  }
	  r *= cm->sft_varied[itpl];
	  if (cm->scaleflag[itpl] != 0) { 
	    nsig += r; 
	  }
	  else { 
	    nbkg += r; 
	  }
	}
	cout << "DumpSBD: " << nsig << " " << nbkg << " " << nobs << endl;
      }
    }
  }  
}

/*-------------------------------------------------------------------------*/
/*    Interface to Joel's genlimit.c program for a Bayesian calcualtion    */
/*    of an upper limit.  Also run pseudoexperiments if need be to compute */
/*    expected limits.                                                     */
/*    Bayesian limit calculation uses test_hypothesis_pe to compute the    */
/*    "Bayesian ensemble" because it should have signal and background     */
/*    components marked, and because it should have all systematic         */
/*    errors included.  It uses null_hypothesis_pe in order to generate    */
/*    pseudoexperiments to compute expected limits however.                */
/*-------------------------------------------------------------------------*/
/* arguments:  beta: credibility level:L  0.95 for 95% CL limits
               sflimit:  the observed limit
               unc:      MC statistical unc. on observed limit
               npx:      Number of pseudoexperiments to run to compute expected limits
               sm2:      -2 sigma expected limit      *put in null pointers for all
               sm1:      -1 sigma expected limit      *five of these to skip the
               smed:     median expected limit        *calculation and speed it up.
               sp1:      +1 sigma expected limit
               sp2:      +2 sigma expected limit
*/

void mclimit_csm::bayes_heinrich_withexpect(Double_t  beta,
                                            Double_t* sflimit,
                                            Double_t* unc,
					    Int_t     npx,
                                            Double_t* sm2,
                                            Double_t* sm1,
                                            Double_t* smed,
                                            Double_t* sp1,
                                            Double_t* sp2)
{
  Int_t nbinstot;
  Int_t i,j,k,ibin,nbinsx,nbinsy,ipx,nens,iens;
  Int_t nchans,ntemplates,itpl;
  csm_channel_model* cm;
  TH1* ht;
  Double_t r;
  int ngl;
  const PRIOR prior=corr;
  vector<Double_t> cslist;
  Int_t nobstot;
  int nglmax;
  double *xgl;
  double *lwgl;

  vector<Double_t> bpploc;
  vector<Double_t> bploc;

  unc = 0;
  nbinstot = 0;

  // figure out the total number of bins in all of our histograms
  nchans = (Int_t) test_hypothesis_pe->channame.size();
  for (i=0;i<nchans;i++)
    {
      nbinsx = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsX();
      nbinsy = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsY();
      nbinstot += nbinsx*nbinsy;
    }

  int* nobs = new int[nbinstot];
  EB* ens = new EB[nbinstot*nmc_req];

  // copy the observed candidates from histograms into nobs -- be sure to have
  // the same association of bins and the flat array as for the model histogram sums

  nobstot = 0;
  ibin = 0;
  for (i=0;i<nchans;i++)
    {
      nbinsx = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsX();
      nbinsy = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsY();
      for (j=0;j<nbinsx;j++)
	{
	  for (k=0;k<nbinsy;k++)
	    {
	      if (nbinsy==1)
		{ nobs[ibin] = (Int_t) nearbyint(datahist[i]->GetBinContent(j+1)); }
	      else
		{ nobs[ibin] = (Int_t) nearbyint(datahist[i]->GetBinContent(j+1,k+1)); }
	      //cout << ibin << "    " << nobs[ibin] << endl;
	      nobstot += nobs[ibin];
	      ibin++;
	    }
	}
    }  

  // The prior ensemble is constructed in the same way mclimit_csm does pseudoexperiments

  iens = 0;
  for (ipx=0;ipx<nmc_req;ipx++)
    {
      test_hypothesis_pe->varysyst();
      for (i=0;i<nchans;i++)
	{
	  cm = test_hypothesis_pe->chanmodel[i];
	  ntemplates = (Int_t) cm->histotemplate.size();
	  nbinsx = cm->histotemplate[0]->GetNbinsX();
	  nbinsy = cm->histotemplate[0]->GetNbinsY();
	  for (j=0;j<nbinsx;j++)
	    {
	      for (k=0;k<nbinsy;k++)
		{
                  ens[iens].e = 0;
                  ens[iens].b = 0;
		  for(itpl=0;itpl<ntemplates;itpl++)
		    {
		      ht = cm->histotemplate_varied[itpl];
		      if (nbinsy==1)
			{ r = ht->GetBinContent(j+1); }
		      else
			{r = ht->GetBinContent(j+1,k+1); }
		      if (cm->poissflag[itpl] == CSM_POISSON_BINERR)
			{ r = fBayesRandom->Poisson(r); }
		      else if (cm->poissflag[itpl] == CSM_GAUSSIAN_BINERR)
			{ 
			  double histerr,edraw;
			  if (nbinsy==1)
			    { histerr = ht->GetBinError(j+1);}
			  else
			    { histerr = ht->GetBinError(j+1,k+1);}
			  do
			    { edraw = fBayesRandom->Gaus(0,histerr); }
			  while (edraw+r<r*1E-6); // don't let it hit zero or go negative.
                                                  // if r is already zero this won't get stuck in a loop
			  r += edraw;
			}
		      r *= cm->sft_varied[itpl];
		      if (cm->scaleflag[itpl] != 0)
			{ ens[iens].e += r; }
		      else
			{ ens[iens].b += r; }
		    }
                  if (ens[iens].b<=0) {ens[iens].b = 1E-9;}
		  //cout << iens << " " << ens[iens].b << " " << ens[iens].e << endl;
		  iens++;
		}
	    }
	}
    }

  //be generous here -- we really just need nobstot/2 entries here,
  //but this memory is fairly inexpensive.  We will enlarge these arrays
  //later if the need arises.

  nglmax = nobstot;
  if (nglmax<10000) {nglmax = 10000;}
  xgl = new double[nglmax];
  lwgl = new double[nglmax];

  nens = nmc_req;
  ngl = 0;

  *sflimit = 0;
  if (bayesintegralmethod == CSM_BAYESINTEGRAL_JOEL)
    {
       *sflimit = (Double_t) cslimit(beta,nbinstot,nens,nobs,ens,&ngl,xgl,lwgl,prior,unc);
    }
  else
    {
      setdlcsn(nbinstot,nens,nobs,ens);
    }

  // make a vector of the posterior likelihood function
  if (bayes_interval_step > 0)
    {
      if ( (bayes_interval_end-bayes_interval_begin)>0 )
	{
	  bayes_posterior.clear();
	  bayes_posterior_points.clear();
	  Double_t b,p;
	  for (b=bayes_interval_begin;b<=bayes_interval_end;b += bayes_interval_step)
	    {
	      p = cspdf(b,1.0,nbinstot,nens,nobs,ens,prior);
	      bayes_posterior.push_back(p);
	      bploc.push_back(p);
	      bayes_posterior_points.push_back(b);
	      bpploc.push_back(b);
	    }
	  int jsiz = bayes_posterior.size();
	  int ic;
	  Double_t bptot = 0;
	  for (ic=0;ic<jsiz;ic++) bptot += bayes_posterior[ic];
	  if (bptot>0)
	    {
    	      Double_t scaleb = 1.0/(bptot*bayes_interval_step);
	      bhnorm = scaleb;
	      for (ic=0;ic<jsiz;ic++) bayes_posterior[ic] *= scaleb;
	    }
	  if (bayesintegralmethod == CSM_BAYESINTEGRAL_QUICK)
	    {
	       *sflimit = quickbint(beta);
	    }
	}
    }


  // compute expected limits

  cslist.clear();
  TH1** pdarray = new TH1*[nchans];
  char *pdname;

  int* nobslist = new int[nbinstot*npx];
  Int_t* nobstotlist = new Int_t[npx];
  Int_t* nobsindex = new Int_t[npx];

  for (i=0;i<(Int_t) null_hypothesis_pe->channame.size(); i++)
    {
      pdname = new char[strlen(test_hypothesis_pe->channame[i])+strlen(" pseudodata ")];
      strcpy(pdname,null_hypothesis_pe->channame[i]);
      strcat(pdname," pseudodata");
      pdarray[i] = (TH1*) null_hypothesis_pe->chanmodel[i]->histotemplate[0]->Clone(pdname);
      delete [] pdname;
    }
  for (ipx=0;ipx<npx;ipx++)
    {
      null_hypothesis_pe->single_pseudoexperiment(pdarray);

      nobstotlist[ipx] = 0;
      ibin = 0;
      for (i=0;i<nchans;i++)
        {
          nbinsx = pdarray[i]->GetNbinsX();
          nbinsy = pdarray[i]->GetNbinsY();
          for (j=0;j<nbinsx;j++)
            {
              for (k=0;k<nbinsy;k++)
                {
		  if (nbinsy==1)
		    { nobslist[ibin+ipx*nbinstot] = (Int_t) nearbyint(pdarray[i]->GetBinContent(j+1)); }
		  else
		    { nobslist[ibin+ipx*nbinstot] = (Int_t) nearbyint(pdarray[i]->GetBinContent(j+1,k+1)); }
		  nobstotlist[ipx] += nobslist[ibin+ipx*nbinstot];
		  ibin++;
                }
            }
        } 
    }
  TMath::Sort(npx,nobstotlist,nobsindex,kTRUE);

  if (nglmax < nobstotlist[nobsindex[0]]/2 + 1)
    {
      nglmax = nobstotlist[nobsindex[0]]/2 + 1; 
      delete[] xgl;
      delete[] lwgl;
      xgl = new double[nglmax];
      lwgl = new double[nglmax];
    }

  ngl = 0;
  for (ipx=0;ipx<npx;ipx++) {
    Double_t p=0;
    if (bayesintegralmethod == CSM_BAYESINTEGRAL_JOEL) {
      if (ipx>0) { 
	if (nobstotlist[nobsindex[ipx]] != nobstotlist[nobsindex[ipx-1]]) { 
	  ngl = 0;
	}
      }
      p = cslimit(beta,nbinstot,nens,&(nobslist[nbinstot*nobsindex[ipx]]),ens,&ngl,xgl,lwgl,prior,unc);
    }
    else {
      setdlcsn(nbinstot,nens,&(nobslist[nbinstot*nobsindex[ipx]]),ens);

	  // make a vector of the posterior likelihood function
      if (bayes_interval_step > 0) {
	if ( (bayes_interval_end-bayes_interval_begin)>0 ) {
	  bayes_posterior.clear();
	  bayes_posterior_points.clear();
	  Double_t b,p1;
	  for (b=bayes_interval_begin;b<=bayes_interval_end;b += bayes_interval_step) {
	    p1 = cspdf(b,1.0,nbinstot,nens,&(nobslist[nbinstot*nobsindex[ipx]]),ens,prior);
	    bayes_posterior.push_back(p1);
	    bayes_posterior_points.push_back(b);
	  }
	  int jsiz = bayes_posterior.size();
	  int ic;
	  Double_t bptot = 0;
	  for (ic=0;ic<jsiz;ic++) bptot += bayes_posterior[ic];
	  if (bptot>0) {
	    Double_t scaleb = 1.0/(bptot*bayes_interval_step);
	    bhnorm = scaleb;
	    for (ic=0;ic<jsiz;ic++) bayes_posterior[ic] *= scaleb;
	  }
	}
	p = quickbint(beta);
      }
    }


    if (pxprintflag) {
      cout << "bayespx: " << p << endl;
    }
    if (bayes_pseudoexperiment_limits != 0) bayes_pseudoexperiment_limits->Fill(p);
    cslist.push_back(p);
  }
  std::sort(cslist.begin(),cslist.end());

// // 2011-05-02 P.Murat: debug


//   for (i=0; i< cslist.size(); i++) {
//     printf ("[bayes_heinrich_withexpect] i=%5i  cslist[i]=%12.5f \n",i,cslist[i]);
//   }

  i =  (Int_t) nearbyint(npx*MCLIMIT_CSM_MCLM2S);
  if (i<0) 
    { i=0; }
  if (i>npx-1)
    { i=npx-1; }
  *sm2 = cslist[i];

  i =  (Int_t) nearbyint(npx*MCLIMIT_CSM_MCLM1S);
  if (i<0) 
    { i=0; }
  if (i>npx-1)
    { i=npx-1; }
  *sm1 = cslist[i];

  i =  (Int_t) nearbyint(npx*MCLIMIT_CSM_MCLMED);
  if (i<0) 
    { i=0; }
  if (i>npx-1)
    { i=npx-1; }
  *smed = cslist[i];

  i =  (Int_t) nearbyint(npx*MCLIMIT_CSM_MCLP1S);
  if (i<0) 
    { i=0; }
  if (i>npx-1)
    { i=npx-1; }
  *sp1 = cslist[i];

  i =  (Int_t) nearbyint(npx*MCLIMIT_CSM_MCLP2S);
  if (i<0) 
    { i=0; }
  if (i>npx-1)
    { i=npx-1; }
  *sp2 = cslist[i];

  for (i=0;i<(Int_t) null_hypothesis_pe->channame.size(); i++) {
    delete pdarray[i];
  }
  delete[] pdarray;
  delete[] nobslist;
  delete[] nobsindex;
  delete[] nobstotlist;
  delete[] nobs;
  delete[] ens;
  delete[] xgl;
  delete[] lwgl;

  // copy the observed posterior curve into the externally visible vectors

  bayes_posterior.clear();
  bayes_posterior_points.clear();
  for (int k=0;k< (int) bpploc.size();k++) {
    bayes_posterior.push_back(bploc[k]);
    bayes_posterior_points.push_back(bpploc[k]);
  }
}

/*-------------------------------------------------------------------------*/
//  Same thing as above, but only compute observed limit (much quicker)
/*-------------------------------------------------------------------------*/

void mclimit_csm::bayes_heinrich(Double_t beta,
                                 Double_t* sflimit,
                                 Double_t* unc)
{
  Int_t nbinstot;
  Int_t i,j,k,ibin,nbinsx,nbinsy,ipx,nens,iens;
  Int_t nchans,ntemplates,itpl;
  csm_channel_model* cm;
  TH1* ht;
  Double_t r;
  int ngl;
  const PRIOR prior=corr;
  Int_t nobstot;
  int nglmax;

  unc = 0;
  nbinstot = 0;

  // figure out the total number of bins in all of our histograms
  nchans = (Int_t) test_hypothesis_pe->channame.size();
  for (i=0;i<nchans;i++)
    {
      nbinsx = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsX();
      nbinsy = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsY();
      nbinstot += nbinsx*nbinsy;
    }

  int* nobs = new int[nbinstot];
  EB* ens = new EB[nbinstot*nmc_req];

  // copy the observed candidates from histograms into nobs -- be sure to have
  // the same association of bins and the flat array as for the model histogram sums

  nobstot = 0;
  ibin = 0;
  for (i=0;i<nchans;i++)
    {
      nbinsx = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsX();
      nbinsy = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsY();
      for (j=0;j<nbinsx;j++)
	{
	  for (k=0;k<nbinsy;k++)
	    {
	      if (nbinsy==1)
		{ nobs[ibin] = (Int_t) nearbyint(datahist[i]->GetBinContent(j+1)); }
	      else
		{ nobs[ibin] = (Int_t) nearbyint(datahist[i]->GetBinContent(j+1,k+1)); }
	      nobstot += nobs[ibin];
	      ibin++;
	    }
	}
    }  

  // The prior ensemble is constructed in the same way mclimit_csm does pseudoexperiments

  iens = 0;
  for (ipx=0;ipx<nmc_req;ipx++)
    {
      test_hypothesis_pe->varysyst();
      for (i=0;i<nchans;i++)
	{
	  cm = test_hypothesis_pe->chanmodel[i];
	  ntemplates = (Int_t) cm->histotemplate.size();
	  nbinsx = cm->histotemplate[0]->GetNbinsX();
	  nbinsy = cm->histotemplate[0]->GetNbinsY();
	  for (j=0;j<nbinsx;j++)
	    {
	      for (k=0;k<nbinsy;k++)
		{
                  ens[iens].e = 0;
                  ens[iens].b = 0;
		  for(itpl=0;itpl<ntemplates;itpl++)
		    {
		      ht = cm->histotemplate_varied[itpl];
		      if (nbinsy==1)
			{ r = ht->GetBinContent(j+1); }
		      else
			{r = ht->GetBinContent(j+1,k+1); }
		      if (cm->poissflag[itpl] == CSM_POISSON_BINERR)
			{ r = fBayesRandom->Poisson(r); }
		      else if (cm->poissflag[itpl] == CSM_GAUSSIAN_BINERR)
			{ 
			  double histerr,edraw;
			  if (nbinsy==1)
			    { histerr = ht->GetBinError(j+1);}
			  else
			    { histerr = ht->GetBinError(j+1,k+1);}
			  do
			    { edraw = fBayesRandom->Gaus(0,histerr); }
			  while (edraw+r<r*1E-6); // don't let it hit zero or go negative.
			  r += edraw;
			}
		      r *= cm->sft_varied[itpl];
		      if (cm->scaleflag[itpl] != 0)
			{ ens[iens].e += r; }
		      else
			{ ens[iens].b += r; }
		    }
                  if (ens[iens].b<=0) {ens[iens].b = 1E-9;}
		  iens++;
		}
	    }
	}
    }

  //be generous here -- we really just need nobstot/2 entries here,
  //but this memory is fairly inexpensive.

  nglmax = nobstot;
  if (nglmax<10000) {nglmax = 10000;}
  double* xgl = new double[nglmax];
  double* lwgl = new double[nglmax];

  nens = nmc_req;
  ngl = 0;

  *sflimit = 0;
  if (bayesintegralmethod == CSM_BAYESINTEGRAL_JOEL)
    {
       *sflimit = (Double_t) cslimit(beta,nbinstot,nens,nobs,ens,&ngl,xgl,lwgl,prior,unc);
    }
  else
    {
      setdlcsn(nbinstot,nens,nobs,ens);
    }

  //cout << "ngl: " << ngl << " " << nobstot << endl;

  // make a vector of the posterior likelihood function
  if (bayes_interval_step > 0)
    {
      if ( (bayes_interval_end-bayes_interval_begin)>0 )
	{
	  bayes_posterior.clear();
	  bayes_posterior_points.clear();
	  Double_t b,p;
	  for (b=bayes_interval_begin;b<=bayes_interval_end;b += bayes_interval_step)
	    {
	      p = cspdf(b,1.0,nbinstot,nens,nobs,ens,prior);
	      bayes_posterior.push_back(p);
	      bayes_posterior_points.push_back(b);
	    }
	  int jsiz = bayes_posterior.size();
	  int ic;
	  Double_t bptot = 0;
	  for (ic=0;ic<jsiz;ic++) bptot += bayes_posterior[ic];
	  if (bptot>0)
	    {
    	      Double_t scaleb = 1.0/(bptot*bayes_interval_step);
	      bhnorm = scaleb;
	      for (ic=0;ic<jsiz;ic++) bayes_posterior[ic] *= scaleb;
	    }
      if (bayesintegralmethod == CSM_BAYESINTEGRAL_QUICK)
	{
	  *sflimit = quickbint(beta);
	}
	}
    }

  delete[] nobs;
  delete[] ens;
  delete[] xgl;
  delete[] lwgl;
}


// Fit a cross section using Joel Heinrich's marginalized posterior
// use testhyp_pe for this fit (has all nuisance parameters defined)

void mclimit_csm::bh_xsfit(Double_t *xsfit, Double_t *downerr, Double_t *uperr)
{
  Int_t nbinstot;
  Int_t i,j,k,ibin,nbinsx,nbinsy,ipx,nens,iens;
  Int_t nchans,ntemplates,itpl;
  csm_channel_model* cm;
  TH1* ht;
  Double_t r;
  //  int ngl;
  const PRIOR prior=flat;
  Int_t nobstot;

  *xsfit = 0.0;
  *downerr = 0.0;
  *uperr = 0.0;
  
  nbinstot = 0;

  // figure out the total number of bins in all of our histograms
  nchans = (Int_t) test_hypothesis_pe->channame.size();
  for (i=0;i<nchans;i++)
    {
      nbinsx = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsX();
      nbinsy = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsY();
      nbinstot += nbinsx*nbinsy;
    }

  int* nobs = new int[nbinstot];
  EB* ens = new EB[nbinstot*nmc_req];

  // copy the observed candidates from histograms into nobs -- be sure to have
  // the same association of bins and the flat array as for the model histogram sums

  nobstot = 0;
  ibin = 0;
  for (i=0;i<nchans;i++)
    {
      nbinsx = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsX();
      nbinsy = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsY();
      for (j=0;j<nbinsx;j++)
	{
	  for (k=0;k<nbinsy;k++)
	    {
	      if (nbinsy==1)
		{ nobs[ibin] = (Int_t) nearbyint(datahist[i]->GetBinContent(j+1)); }
	      else
		{ nobs[ibin] = (Int_t) nearbyint(datahist[i]->GetBinContent(j+1,k+1)); }
	      nobstot += nobs[ibin];
	      ibin++;
	    }
	}
    }  

  // The prior ensemble is constructed in the same way mclimit_csm does pseudoexperiments

  iens = 0;
  for (ipx=0;ipx<nmc_req;ipx++)
    {
      test_hypothesis_pe->varysyst();
      for (i=0;i<nchans;i++)
	{
	  cm = test_hypothesis_pe->chanmodel[i];
	  ntemplates = (Int_t) cm->histotemplate.size();
	  nbinsx = cm->histotemplate[0]->GetNbinsX();
	  nbinsy = cm->histotemplate[0]->GetNbinsY();
	  for (j=0;j<nbinsx;j++)
	    {
	      for (k=0;k<nbinsy;k++)
		{
                  ens[iens].e = 0;
                  ens[iens].b = 0;
		  for(itpl=0;itpl<ntemplates;itpl++)
		    {
		      ht = cm->histotemplate_varied[itpl];
		      if (nbinsy==1)
			{ r = ht->GetBinContent(j+1); }
		      else
			{r = ht->GetBinContent(j+1,k+1); }
		      if (cm->poissflag[itpl] == CSM_POISSON_BINERR)
			{ r = fBayesRandom->Poisson(r); }
		      else if (cm->poissflag[itpl] == CSM_GAUSSIAN_BINERR)
			{ 
			  double histerr,edraw;
			  if (nbinsy==1)
			    { histerr = ht->GetBinError(j+1);}
			  else
			    { histerr = ht->GetBinError(j+1,k+1);}
			  do
			    { edraw = fBayesRandom->Gaus(0,histerr); }
			  while (edraw+r<r*1E-6); // don't let it hit zero or go negative.
			  r += edraw;
			}
		      r *= cm->sft_varied[itpl];
		      if (cm->scaleflag[itpl] != 0)
			{ ens[iens].e += r; }
		      else
			{ ens[iens].b += r; }
		    }
                  if (ens[iens].b<=0) {ens[iens].b = 1E-9;}
		  iens++;
		}
	    }
	}
    }

  nens = nmc_req;
  //  ngl = 0;

  //cout << "ngl: " << ngl << " " << nobstot << endl;


  // make a vector of the posterior likelihood function
  if (bayes_interval_step > 0)
    {
      if ( (bayes_interval_end-bayes_interval_begin)>0 )
	{
          double lmax = 0;
          int imax=0;
	  bayes_posterior.clear();
	  bayes_posterior_points.clear();
	  Double_t b,p;
	  int k=0;
	  for (b=bayes_interval_begin;b<=bayes_interval_end;b += bayes_interval_step)
	    {
	      p = cspdf(b,1.0,nbinstot,nens,nobs,ens,prior);
	      bayes_posterior.push_back(p);
	      bayes_posterior_points.push_back(b);
	      if (p>lmax)
		{ 
		  lmax=p; 
		  imax = k;
		  *xsfit = b;
		}
	      k++;
	    }

	  // normalize to unit area, using the trapezoidal rule

	  int jsiz = bayes_posterior.size();
	  int ic;
	  Double_t bptot = 0;
	  for (ic=0;ic<jsiz;ic++) bptot += bayes_posterior[ic];
	  bptot -= 0.5*bayes_posterior[0];
	  bptot -= 0.5*bayes_posterior.back();
	  if (bptot>0)
	    {
    	      Double_t scaleb = 1.0/(bptot);
	      bhnorm = scaleb;
	      for (ic=0;ic<jsiz;ic++) bayes_posterior[ic] *= scaleb;
	    }

	  // start at the maximum, and integrate left and right until we get 68%

	  int ilow=imax;
	  int ihigh=imax;
	  double psum = bayes_posterior[imax];
	  double psumtrap = 0.0;
	  do 
	    {
	      int ibest = ilow;
	      double bpbest = 0.0;
	      if (ilow>0)
		{
		  ibest = ilow-1;
		  bpbest = bayes_posterior[ibest];
		}
	      if (ihigh<jsiz-1)
		{
		  if (bayes_posterior[ihigh+1] > bpbest)
		    {
  	   	      ibest = ihigh+1;
		      bpbest = bayes_posterior[ibest];
		    }
		}
	      psum += bpbest;
	      if (ibest == ilow-1) 
		{ ilow --; }
	      else
		{ ihigh ++; }
      	      psumtrap = psum - 0.5*bayes_posterior[ilow] - 0.5*bayes_posterior[ihigh];
	    }
	  while (psumtrap < 0.68);
	  *downerr = *xsfit-bayes_posterior_points[ilow];
	  *uperr = bayes_posterior_points[ihigh] - *xsfit;
	}
    }

  delete[] nobs;
  delete[] ens;
}


// quick and dirty posterior integrator -- assumes bayes_posterior and bayes_posterior_points
// have been filled.

Double_t mclimit_csm::quickbint(Double_t beta)
{
  // normalize to unit area, using the trapezoidal rule

  int jsiz = bayes_posterior.size();
  int ic;
  Double_t bptot = 0;
  for (ic=0;ic<jsiz;ic++) bptot += bayes_posterior[ic];
  bptot -= 0.5*bayes_posterior[0];
  bptot -= 0.5*bayes_posterior.back();
  if (bptot>0)
    {
      Double_t scaleb = 1.0/(bptot);
      for (ic=0;ic<jsiz;ic++) bayes_posterior[ic] *= scaleb;
    }

  // count bacwards from the end until we get 1-beta of the integral
  // do a little trapezoidal trick to add back part of the interval if we overstep

  bptot = 0;
  for (ic=jsiz-1;ic>=0;ic--)
    {
      bptot += bayes_posterior[ic];
      if (bptot >= 1-beta) return(bayes_posterior_points[ic]+
				  bayes_interval_step*((bptot-(1-beta))/bayes_posterior[ic]));
    }
  return(0);
}

// expected fitted cross section distribuitons -- run pseudoexperiments from the test hypothesis
// and give distributions of fitted values.

void mclimit_csm::bh_xsfit_expect(Int_t npx, Double_t *xsfitavg, Double_t *m2s, 
				  Double_t *m1s, Double_t *med, Double_t *p1s, Double_t *p2s)
{
  Int_t nbinstot;
  Int_t i,j,k,ibin,nbinsx,nbinsy,ipx,nens,iens;
  Int_t nchans,ntemplates,itpl;
  csm_channel_model* cm;
  TH1* ht;
  Double_t r;
  const PRIOR prior=flat;
  vector<Double_t> xslist;
  double xsfit=0;
  double downerr=0;
  double uperr=0;

  nbinstot = 0;

  // figure out the total number of bins in all of our histograms
  nchans = (Int_t) test_hypothesis_pe->channame.size();
  for (i=0;i<nchans;i++)
    {
      nbinsx = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsX();
      nbinsy = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsY();
      nbinstot += nbinsx*nbinsy;
    }

  int* nobs = new int[nbinstot];
  EB* ens = new EB[nbinstot*nmc_req];

  // The prior ensemble is constructed in the same way mclimit_csm does pseudoexperiments

  iens = 0;
  for (ipx=0;ipx<nmc_req;ipx++)
    {
      test_hypothesis_pe->varysyst();
      for (i=0;i<nchans;i++)
	{
	  cm = test_hypothesis_pe->chanmodel[i];
	  ntemplates = (Int_t) cm->histotemplate.size();
	  nbinsx = cm->histotemplate[0]->GetNbinsX();
	  nbinsy = cm->histotemplate[0]->GetNbinsY();
	  for (j=0;j<nbinsx;j++)
	    {
	      for (k=0;k<nbinsy;k++)
		{
                  ens[iens].e = 0;
                  ens[iens].b = 0;
		  for(itpl=0;itpl<ntemplates;itpl++)
		    {
		      ht = cm->histotemplate_varied[itpl];
		      if (nbinsy==1)
			{ r = ht->GetBinContent(j+1); }
		      else
			{r = ht->GetBinContent(j+1,k+1); }
		      if (cm->poissflag[itpl] == CSM_POISSON_BINERR)
			{ r = fBayesRandom->Poisson(r); }
		      else if (cm->poissflag[itpl] == CSM_GAUSSIAN_BINERR)
			{ 
			  double histerr,edraw;
			  if (nbinsy==1)
			    { histerr = ht->GetBinError(j+1);}
			  else
			    { histerr = ht->GetBinError(j+1,k+1);}
			  do
			    { edraw = fBayesRandom->Gaus(0,histerr); }
			  while (edraw+r<r*1E-6); // don't let it hit zero or go negative.
                                                  // if r is already zero this won't get stuck in a loop
			  r += edraw;
			}
		      r *= cm->sft_varied[itpl];
		      if (cm->scaleflag[itpl] != 0)
			{ ens[iens].e += r; }
		      else
			{ ens[iens].b += r; }
		    }
                  if (ens[iens].b<=0) {ens[iens].b = 1E-9;}
		  //cout << iens << " " << ens[iens].b << " " << ens[iens].e << endl;
		  iens++;
		}
	    }
	}
    }

  nens = nmc_req;

  // run pseudoexperiments and fit the cross section for each one. -- Use test_hypothesis
  // to generate the pseudoexperiments.

  xslist.clear();
  TH1** pdarray = new TH1*[nchans];
  char *pdname;

  int* nobslist = new int[nbinstot*npx];
  Int_t* nobstotlist = new Int_t[npx];
  Int_t* nobsindex = new Int_t[npx];

  for (i=0;i<(Int_t) test_hypothesis->channame.size(); i++)
    {
      pdname = new char[strlen(test_hypothesis->channame[i])+strlen(" pseudodata ")];
      strcpy(pdname,test_hypothesis->channame[i]);
      strcat(pdname," pseudodata");
      pdarray[i] = (TH1*) test_hypothesis->chanmodel[i]->histotemplate[0]->Clone(pdname);
      delete [] pdname;
    }

  // don't really need to store and sort all of these, but code is take from bayes_heinrich_withexpect
  // which needed that for optimization

  for (ipx=0;ipx<npx;ipx++)
    {
      test_hypothesis->single_pseudoexperiment(pdarray);

      nobstotlist[ipx] = 0;
      ibin = 0;
      for (i=0;i<nchans;i++)
        {
          nbinsx = pdarray[i]->GetNbinsX();
          nbinsy = pdarray[i]->GetNbinsY();
          for (j=0;j<nbinsx;j++)
            {
              for (k=0;k<nbinsy;k++)
                {
		  if (nbinsy==1)
		    { nobslist[ibin+ipx*nbinstot] = (Int_t) nearbyint(pdarray[i]->GetBinContent(j+1)); }
		  else
		    { nobslist[ibin+ipx*nbinstot] = (Int_t) nearbyint(pdarray[i]->GetBinContent(j+1,k+1)); }
		  nobstotlist[ipx] += nobslist[ibin+ipx*nbinstot];
		  ibin++;
                }
            }
        } 
    }
  TMath::Sort(npx,nobstotlist,nobsindex,kTRUE);

  for (ipx=0;ipx<npx;ipx++)
    {
      // make a vector of the posterior likelihood function
      if (bayes_interval_step > 0)
	{
	  if ( (bayes_interval_end-bayes_interval_begin)>0 )
	    {
	      double lmax = 0;
	      int imax=0;
	      bayes_posterior.clear();
	      bayes_posterior_points.clear();
	      Double_t b,p;
	      int k=0;
	      for (b=bayes_interval_begin;b<=bayes_interval_end;b += bayes_interval_step)
		{
		  p = cspdf(b,1.0,nbinstot,nens,&(nobslist[nbinstot*nobsindex[ipx]]),ens,prior);
		  bayes_posterior.push_back(p);
		  bayes_posterior_points.push_back(b);
		  if (p>lmax)
		    { 
		      lmax=p; 
		      imax = k;
		      xsfit = b;
		    }
		  k++;
		}

	      // normalize to unit area, using the trapezoidal rule

	      int jsiz = bayes_posterior.size();
	      int ic;
	      Double_t bptot = 0;
	      for (ic=0;ic<jsiz;ic++) bptot += bayes_posterior[ic];
	      bptot -= 0.5*bayes_posterior[0];
	      bptot -= 0.5*bayes_posterior.back();
	      if (bptot>0)
		{
		  Double_t scaleb = 1.0/(bptot);
                  bhnorm = scaleb;
		  for (ic=0;ic<jsiz;ic++) bayes_posterior[ic] *= scaleb;
		}

	      // start at the maximum, and integrate left and right until we get 68%

	      int ilow=imax;
	      int ihigh=imax;
  	      double psum = bayes_posterior[imax];
	      double psumtrap = 0.0;
	      do 
		{
		  int ibest = ilow;
		  double bpbest = 0.0;
		  if (ilow>0)
		    {
		      ibest = ilow-1;
		      bpbest = bayes_posterior[ibest];
		    }
		  if (ihigh<jsiz-1)
		    {
		      if (bayes_posterior[ihigh+1] > bpbest)
			{
  	   	          ibest = ihigh+1;
		          bpbest = bayes_posterior[ibest];
			}
		    }
		  psum += bpbest;
		  if (ibest == ilow-1) 
		    { ilow --; }
		  else
		    { ihigh ++; }
		  psumtrap = psum - 0.5*bayes_posterior[ilow] - 0.5*bayes_posterior[ihigh];
		}
	      while (psumtrap < 0.68);
	      downerr = xsfit-bayes_posterior_points[ilow];
	      uperr = bayes_posterior_points[ihigh] - xsfit;
	    }
	}

      //Double_t p = cslimit(beta,nbinstot,nens,&(nobslist[nbinstot*nobsindex[ipx]]),ens,&ngl,xgl,lwgl,prior,unc);
      if (pxprintflag)
	{
	  cout << "Marginalized fit px: " << xsfit << " + " << uperr << " - " << downerr << endl;
	}
      xslist.push_back(xsfit);
      *xsfitavg += xsfit;
    }
  if (xslist.size()>0) 
    {
       *xsfitavg /= xslist.size();
    }

  std::sort(xslist.begin(),xslist.end());
  i =  (Int_t) nearbyint(npx*MCLIMIT_CSM_MCLM2S);
  if (i<0) 
    { i=0; }
  if (i>npx-1)
    { i=npx-1; }
  *m2s = xslist[i];

  i =  (Int_t) nearbyint(npx*MCLIMIT_CSM_MCLM1S);
  if (i<0) 
    { i=0; }
  if (i>npx-1)
    { i=npx-1; }
  *m1s = xslist[i];

  i =  (Int_t) nearbyint(npx*MCLIMIT_CSM_MCLMED);
  if (i<0) 
    { i=0; }
  if (i>npx-1)
    { i=npx-1; }
  *med = xslist[i];

  i =  (Int_t) nearbyint(npx*MCLIMIT_CSM_MCLP1S);
  if (i<0) 
    { i=0; }
  if (i>npx-1)
    { i=npx-1; }
  *p1s = xslist[i];

  i =  (Int_t) nearbyint(npx*MCLIMIT_CSM_MCLP2S);
  if (i<0) 
    { i=0; }
  if (i>npx-1)
    { i=npx-1; }
  *p2s = xslist[i];

  for (i=0;i<(Int_t) null_hypothesis_pe->channame.size(); i++)
    {
      delete pdarray[i];
    }
  delete[] pdarray;
  delete[] nobslist;
  delete[] nobsindex;
  delete[] nobstotlist;
  delete[] nobs;
  delete[] ens;
}

// scan over two signals -- always assume they are in the same order in all channels.  Assume
// a flat prior in the two signals and print out the marginalized posterior.
// If more than two signals are present, the first one in  each channel is called signal 1,
// and the sum of all others is called signal 2


void mclimit_csm::bh_2d_scan(Double_t s1low, Double_t s1high, Double_t ds1,
                   Double_t s2low, Double_t s2high, Double_t ds2)
{
  Int_t nbinstot;
  Int_t i,j,k,ibin,nbinsx,nbinsy,ipx,nens,iens;
  Int_t nchans,ntemplates,itpl;
  csm_channel_model* cm;
  TH1* ht;
  Double_t r;
  //  int ngl;
  Int_t nobstot;
  const PRIOR prior=flat;

  nbinstot = 0;

  // figure out the total number of bins in all of our histograms
  nchans = (Int_t) test_hypothesis_pe->channame.size();
  for (i=0;i<nchans;i++)
    {
      nbinsx = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsX();
      nbinsy = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsY();
      nbinstot += nbinsx*nbinsy;
    }

  int* nobs = new int[nbinstot];
  EB2D* ens = new EB2D[nbinstot*nmc_req];

  // copy the observed candidates from histograms into nobs -- be sure to have
  // the same association of bins and the flat array as for the model histogram sums

  nobstot = 0;
  ibin = 0;
  for (i=0;i<nchans;i++)
    {
      nbinsx = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsX();
      nbinsy = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsY();
      for (j=0;j<nbinsx;j++)
	{
	  for (k=0;k<nbinsy;k++)
	    {
	      if (nbinsy==1)
		{ nobs[ibin] = (Int_t) nearbyint(datahist[i]->GetBinContent(j+1)); }
	      else
		{ nobs[ibin] = (Int_t) nearbyint(datahist[i]->GetBinContent(j+1,k+1)); }
	      nobstot += nobs[ibin];
	      ibin++;
	    }
	}
    }  

  // The prior ensemble is constructed in the same way mclimit_csm does pseudoexperiments
  // but split up the first and subsequent signals.

  iens = 0;
  for (ipx=0;ipx<nmc_req;ipx++)
    {
      test_hypothesis_pe->varysyst();
      for (i=0;i<nchans;i++)
	{
	  cm = test_hypothesis_pe->chanmodel[i];
	  ntemplates = (Int_t) cm->histotemplate.size();
	  nbinsx = cm->histotemplate[0]->GetNbinsX();
	  nbinsy = cm->histotemplate[0]->GetNbinsY();
	  for (j=0;j<nbinsx;j++)
	    {
	      for (k=0;k<nbinsy;k++)
		{
                  ens[iens].e1 = 0;
                  ens[iens].e2 = 0;
                  ens[iens].b = 0;
		  int isc = 0;
		  for(itpl=0;itpl<ntemplates;itpl++)
		    {
		      ht = cm->histotemplate_varied[itpl];
		      if (nbinsy==1)
			{ r = ht->GetBinContent(j+1); }
		      else
			{r = ht->GetBinContent(j+1,k+1); }
		      if (cm->poissflag[itpl] == CSM_POISSON_BINERR)
			{ r = fBayesRandom->Poisson(r); }
		      else if (cm->poissflag[itpl] == CSM_GAUSSIAN_BINERR)
			{ 
			  double histerr,edraw;
			  if (nbinsy==1)
			    { histerr = ht->GetBinError(j+1);}
			  else
			    { histerr = ht->GetBinError(j+1,k+1);}
			  do
			    { edraw = fBayesRandom->Gaus(0,histerr); }
			  while (edraw+r<r*1E-6); // don't let it hit zero or go negative.
			  r += edraw;
			}
		      r *= cm->sft_varied[itpl];
		      if (cm->scaleflag[itpl] != 0)
			{ 
			  if (isc == 0)
			    { 
			      ens[iens].e1 += r;
			      isc++;
			    }
			  else
			    {
			      ens[iens].e2 += r;
			    }
			}
		      else
			{ ens[iens].b += r; }
		    }
                  if (ens[iens].b<=0) {ens[iens].b = 1E-9;}
		  iens++;
		}
	    }
	}
    }

  nens = nmc_req;
  //  ngl = 0;

  //cout << "ngl: " << ngl << " " << nobstot << endl;

  // for now just print out the posterior (normalized to unit sum over the interval chosen)

  // protect against crazy likelihoods underflowing or overflowing us.
  setdlcsn2d(nbinstot,nens,nobs,ens);

  double xsig1,xsig2;
  vector <double> xsig1v;
  vector <double> xsig2v;
  vector <double> postv;
  vector <double> postint;
  double psum = 0;
  int itot;

  for (xsig1=s1low;xsig1<=s1high;xsig1+=ds1)
    {
      for (xsig2=s2low;xsig2<=s2high;xsig2+=ds2)
         {
	   xsig1v.push_back(xsig1);
	   xsig2v.push_back(xsig2);
	   postv.push_back(cspdf2d(xsig1,xsig2,1.0,nbinstot,nens,nobs,ens,prior));
	   psum += postv.back();
	   postint.push_back(0);
	 }
    }
  itot=xsig1v.size();
  if (psum>0)
    {
      for (i=0;i<itot;i++)
        {
	  postv[i] = postv[i]/psum;
        }
    }
  else
    {
      cout << "bh_2d_scan -- normalization failed." << endl;
    }

  // find the best fit, and the 68% and 95% regions -- but only if it's normalized.

  if (psum>0)
    {
      int *idx = new int[itot];
      TMath::Sort(itot,&(postv[0]),idx,kTRUE);
      cout << "Two-Dimensional Maximum Posterior: " << xsig1v[idx[0]] << " " << xsig2v[idx[0]] << endl;
      double ps2 = 0;
      for (i=0;i<itot;i++)
	{
	  ps2 += postv[idx[i]];
          postint[idx[i]] = ps2;
	}
      delete[] idx;
    }

  for (i=0;i<itot;i++)
    {	   
      cout << "bh_2d_scan: " << xsig1v[i] << " " << xsig2v[i] << " " 
	   << postv[i] << " " << postint[i] << endl;
    }

  delete[] nobs;
  delete[] ens;
}

// scan over two signals -- always assume they are in the same order in all channels.  Assume
// a flat prior in the two signals and print out the marginalized posterior.
// If more than two signals are present, the first one in  each channel is called signal 1,
// and the sum of all others is called signal 2

// This routine below runs pseudoexperiments to get the expected distributions.  It does all the
// calculations that bh_2d_scan does for the real data, but for the pseudoexperiments instead, and prints
// out the 2D posterior distributions, cumulative integrals, and best-fit points for each pseudoexperiment.
// As an added feature, you put in the signal scale factors used in the pseudoexperiment generation
//  (testhyp is used for pseudoexperiments, and testhyp_pe -- yes, I know that seems backwards -- is the
//  model used to compute the posterior), relative to the model used to test as s1true and s2true, and it
// will compute 68% and 95% coverage fractions based on the nearest grid point.  Verbose printout is
// supplied.

void mclimit_csm::bh_2d_scan_expect(Int_t npx, Double_t s1true, Double_t s2true, 
                                    Double_t s1low, Double_t s1high, Double_t ds1,
				    Double_t s2low, Double_t s2high, Double_t ds2)
{
  Int_t nbinstot;
  Int_t i,j,k,ibin,nbinsx,nbinsy,ipx,nens,iens;
  Int_t nchans,ntemplates,itpl;
  csm_channel_model* cm;
  TH1* ht;
  Double_t r;
  const PRIOR prior=flat;

  nbinstot = 0;

  // figure out the total number of bins in all of our histograms
  nchans = (Int_t) test_hypothesis_pe->channame.size();
  for (i=0;i<nchans;i++)
    {
      nbinsx = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsX();
      nbinsy = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsY();
      nbinstot += nbinsx*nbinsy;
    }

  int* nobs = new int[nbinstot];
  EB2D* ens = new EB2D[nbinstot*nmc_req];

  // The prior ensemble is constructed in the same way mclimit_csm does pseudoexperiments
  // but split up the first and subsequent signals.

  iens = 0;
  for (ipx=0;ipx<nmc_req;ipx++)
    {
      test_hypothesis_pe->varysyst();
      for (i=0;i<nchans;i++)
	{
	  cm = test_hypothesis_pe->chanmodel[i];
	  ntemplates = (Int_t) cm->histotemplate.size();
	  nbinsx = cm->histotemplate[0]->GetNbinsX();
	  nbinsy = cm->histotemplate[0]->GetNbinsY();
	  for (j=0;j<nbinsx;j++)
	    {
	      for (k=0;k<nbinsy;k++)
		{
                  ens[iens].e1 = 0;
                  ens[iens].e2 = 0;
                  ens[iens].b = 0;
		  int isc = 0;
		  for(itpl=0;itpl<ntemplates;itpl++)
		    {
		      ht = cm->histotemplate_varied[itpl];
		      if (nbinsy==1)
			{ r = ht->GetBinContent(j+1); }
		      else
			{r = ht->GetBinContent(j+1,k+1); }
		      if (cm->poissflag[itpl] == CSM_POISSON_BINERR)
			{ r = fBayesRandom->Poisson(r); }
		      else if (cm->poissflag[itpl] == CSM_GAUSSIAN_BINERR)
			{ 
			  double histerr,edraw;
			  if (nbinsy==1)
			    { histerr = ht->GetBinError(j+1);}
			  else
			    { histerr = ht->GetBinError(j+1,k+1);}
			  do
			    { edraw = fBayesRandom->Gaus(0,histerr); }
			  while (edraw+r<r*1E-6); // don't let it hit zero or go negative.
			  r += edraw;
			}
		      r *= cm->sft_varied[itpl];
		      if (cm->scaleflag[itpl] != 0)
			{ 
			  if (isc == 0)
			    { 
			      ens[iens].e1 += r;
			      isc++;
			    }
			  else
			    {
			      ens[iens].e2 += r;
			    }
			}
		      else
			{ ens[iens].b += r; }
		    }
                  if (ens[iens].b<=0) {ens[iens].b = 1E-9;}
		  iens++;
		}
	    }
	}
    }

  nens = nmc_req;

  // run pseudoexperiments and fit the cross section for each one. -- Use test_hypothesis
  // to generate the pseudoexperiments.

  TH1** pdarray = new TH1*[nchans];
  char *pdname;

  int* nobslist = new int[nbinstot*npx];
  Int_t* nobstotlist = new Int_t[npx];
  Int_t* nobsindex = new Int_t[npx];

  for (i=0;i<(Int_t) test_hypothesis->channame.size(); i++)
    {
      pdname = new char[strlen(test_hypothesis->channame[i])+strlen(" pseudodata ")];
      strcpy(pdname,test_hypothesis->channame[i]);
      strcat(pdname," pseudodata");
      pdarray[i] = (TH1*) test_hypothesis->chanmodel[i]->histotemplate[0]->Clone(pdname);
      delete [] pdname;
    }

  // don't really need to store and sort all of these, but code is take from bayes_heinrich_withexpect
  // which needed that for optimization

  for (ipx=0;ipx<npx;ipx++)
    {
      test_hypothesis->single_pseudoexperiment(pdarray);

      nobstotlist[ipx] = 0;
      ibin = 0;
      for (i=0;i<nchans;i++)
        {
          nbinsx = pdarray[i]->GetNbinsX();
          nbinsy = pdarray[i]->GetNbinsY();
          for (j=0;j<nbinsx;j++)
            {
              for (k=0;k<nbinsy;k++)
                {
		  if (nbinsy==1)
		    { nobslist[ibin+ipx*nbinstot] = (Int_t) nearbyint(pdarray[i]->GetBinContent(j+1)); }
		  else
		    { nobslist[ibin+ipx*nbinstot] = (Int_t) nearbyint(pdarray[i]->GetBinContent(j+1,k+1)); }
		  nobstotlist[ipx] += nobslist[ibin+ipx*nbinstot];
		  ibin++;
                }
            }
        } 
    }
  TMath::Sort(npx,nobstotlist,nobsindex,kTRUE);

  int n68 = 0;
  int n95 = 0;
  int ncov = 0;

  for (ipx=0;ipx<npx;ipx++)
    {

      // protect against crazy likelihoods underflowing or overflowing us.

      setdlcsn2d(nbinstot,nens,&(nobslist[nbinstot*nobsindex[ipx]]),ens);

      double xsig1,xsig2;
      vector <double> xsig1v;
      vector <double> xsig2v;
      vector <double> postv;
      vector <double> postint;
      double psum = 0;
      int itot;
      int ibest = -1;
      double rtest;
      double rbest = 0;

      for (xsig1=s1low;xsig1<=s1high;xsig1+=ds1)
	{
	  for (xsig2=s2low;xsig2<=s2high;xsig2+=ds2)
	    {
	      xsig1v.push_back(xsig1);
	      xsig2v.push_back(xsig2);
	      postv.push_back(cspdf2d(xsig1,xsig2,1.0,nbinstot,nens,&(nobslist[nbinstot*nobsindex[ipx]]),ens,prior));
	      psum += postv.back();
	      postint.push_back(0);

	      rtest = sqrt((xsig1-s1true)*(xsig1-s1true) + (xsig2-s2true)*(xsig2-s2true));
	      if (ibest == -1 || (rtest < rbest) )
		{
		  ibest = xsig1v.size() -1;
		  rbest = rtest;
		}
		  
	    }
	}
      itot=xsig1v.size();
      if (psum>0)
	{
	  for (i=0;i<itot;i++)
	    {
	      postv[i] = postv[i]/psum;
	    }
	}
      else
	{
	  cout << "bh_2d_scan_expect -- normalization failed." << endl;
	}

      if (psum>0)
	{
	  int *idx = new int[itot];
	  TMath::Sort(itot,&(postv[0]),idx,kTRUE);
	  cout << "Two-Dimensional Maximum Posterior px: " << xsig1v[idx[0]] << " " << xsig2v[idx[0]] << endl;
	  double ps2 = 0;
	  for (i=0;i<itot;i++)
	    {
	      ps2 += postv[idx[i]];
	      postint[idx[i]] = ps2;
	    }

	  ncov++;
	  if (postint[ibest]<0.68)
	    {
	      cout << "bh_2d_scan_expect true signal falls within 68 percent region" << endl;
	      n68++;
	    }
	  else
	    {
	      cout << "bh_2d_scan_expect true signal falls outside 68 percent region" << endl;
	    }
	  if (postint[ibest]<0.95)
	    {
	      cout << "bh_2d_scan_expect true signal falls within 95 percent region" << endl;
	      n95++;
	    }
	  else
	    {
	      cout << "bh_2d_scan_expect true signal falls outside 95 percent region" << endl;
	    }
	  
	  delete[] idx;
	}

      for (i=0;i<itot;i++)
	{	   
	  //cout << "bh_2d_scan_expect px: " << xsig1v[i] << " " << xsig2v[i] << " " 
	  //     << postv[i] << " " << postint[i] << endl;
	}
    }

  cout << "bh_2d_scan_expect: n68, ncov, Fraction of time true signal lies within 68 percent region: " <<
    n68 << " " << ncov << " " << ((double) n68)/((double) ncov) << endl;
  cout << "bh_2d_scan_expect: n95, ncov, Fraction of time true signal lies within 95 percent region: " <<
    n95 << " " << ncov << " " << ((double) n95)/((double) ncov) << endl;

  delete[] pdarray;
  delete[] nobslist;
  delete[] nobsindex;
  delete[] nobstotlist;
  delete[] nobs;
  delete[] ens;
}

/*-------------------------------------------------------------------------*/
/* Coverage checker for Joel's Bayesian limit calc.  Based on              */
/* bayes_heinrich_withexpect, but now uses test_hypothesis_pe scaled       */
/* so the signal is at a user-setabble desired rate (a useful test is to   */
/* test the coverage for the rate excluded by bayes_heinrich, but one can  */
/* also test the coverage for any other signal rate, scaling the           */
/* signal in testhyp_pe by sflimit. The coverage should be 100% for zero   */
/* signal, for example.    The px's are done                               */
/* assuming the signal+background is present, and the false exclusion rate */
/* is computed.                                                            */
/*-------------------------------------------------------------------------*/
/* arguments:  beta: credibility level:L  0.95 for 95% CL limits
               sflimit:  INPUT -- desired multiplier on the signal in the testhyp_pe hypothesis
                         for which we'd like to check the coverage.
               npx:      Number of pseudoexperiments to run to compute fales exclusion rate
	       falsex:   false exclusion rate:  Should be no more than 1-beta.

*/

void mclimit_csm::bayes_heinrich_coverage_check(Double_t beta,
                                                Double_t sflimit,
					        Int_t npx,
                                                Double_t* falsex)
{
  Int_t nbinstot;
  Int_t i,j,k,ibin,nbinsx,nbinsy,ipx,nens,iens;
  Int_t nchans,ntemplates,itpl;
  csm_channel_model* cm;
  TH1* ht;
  Double_t r;
  int ngl;
  const PRIOR prior=corr;
  vector<Double_t> cslist;
  Int_t nobstot;
  int nglmax;
  double *xgl;
  double *lwgl;

  nbinstot = 0;

  // figure out the total number of bins in all of our histograms
  nchans = (Int_t) test_hypothesis_pe->channame.size();
  for (i=0;i<nchans;i++)
    {
      nbinsx = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsX();
      nbinsy = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsY();
      nbinstot += nbinsx*nbinsy;
    }

  int* nobs = new int[nbinstot];
  EB* ens = new EB[nbinstot*nmc_req];

  // copy the observed candidates from histograms into nobs -- be sure to have
  // the same association of bins and the flat array as for the model histogram sums

  nobstot = 0;
  ibin = 0;
  for (i=0;i<nchans;i++)
    {
      nbinsx = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsX();
      nbinsy = test_hypothesis_pe->chanmodel[i]->histotemplate[0]->GetNbinsY();
      for (j=0;j<nbinsx;j++)
	{
	  for (k=0;k<nbinsy;k++)
	    {
	      if (nbinsy==1)
		{ nobs[ibin] = (Int_t) nearbyint(datahist[i]->GetBinContent(j+1)); }
	      else
		{ nobs[ibin] = (Int_t) nearbyint(datahist[i]->GetBinContent(j+1,k+1)); }
	      nobstot += nobs[ibin];
	      ibin++;
	    }
	}
    }  

  // The prior ensemble is constructed in the same way mclimit_csm does pseudoexperiments

  iens = 0;
  for (ipx=0;ipx<nmc_req;ipx++)
    {
      test_hypothesis_pe->varysyst();
      for (i=0;i<nchans;i++)
	{
	  cm = test_hypothesis_pe->chanmodel[i];
	  ntemplates = (Int_t) cm->histotemplate.size();
	  nbinsx = cm->histotemplate[0]->GetNbinsX();
	  nbinsy = cm->histotemplate[0]->GetNbinsY();
	  for (j=0;j<nbinsx;j++)
	    {
	      for (k=0;k<nbinsy;k++)
		{
                  ens[iens].e = 0;
                  ens[iens].b = 0;
		  for(itpl=0;itpl<ntemplates;itpl++)
		    {
		      ht = cm->histotemplate_varied[itpl];
		      if (nbinsy==1)
			{ r = ht->GetBinContent(j+1); }
		      else
			{r = ht->GetBinContent(j+1,k+1); }
		      if (cm->poissflag[itpl] == CSM_POISSON_BINERR)
			{ r = fBayesRandom->Poisson(r); }
		      else if (cm->poissflag[itpl] == CSM_GAUSSIAN_BINERR)
			{ 
			  double histerr,edraw;
			  if (nbinsy==1)
			    { histerr = ht->GetBinError(j+1);}
			  else
			    { histerr = ht->GetBinError(j+1,k+1);}
			  do
			    { edraw = fBayesRandom->Gaus(0,histerr); }
			  while (edraw+r<r*1E-6); // don't let it hit zero or go negative.
			  r += edraw;
			}
		      r *= cm->sft_varied[itpl];
		      if (cm->scaleflag[itpl] != 0)
			{ ens[iens].e += r; }
		      else
			{ ens[iens].b += r; }
		    }
                  if (ens[iens].b<=0) {ens[iens].b = 1E-9;}
		  iens++;
		}
	    }
	}
    }

  //be generous here -- we really just need nobstot/2 entries here,
  //but this memory is fairly inexpensive.  We will enlarge these arrays
  //later if the need arises.

  nglmax = nobstot;
  if (nglmax<10000) {nglmax = 10000;}
  xgl = new double[nglmax];
  lwgl = new double[nglmax];

  nens = nmc_req;
  ngl = 0;
  // do not compute sflimit here -- input it instead as an adjustable parameter

  //  *sflimit = (Double_t) cslimit(beta,nbinstot,nens,nobs,ens,&ngl,xgl,lwgl,prior,unc);
  
  // Run signal+background pseudoexperiments at the observed limit, and see what
  // the distribution of limits we get out is.  Limits are computed using the
  // same Bayesian ensemble with the unscaled test_hypothesis_pe and so the
  // limits that are more restrictive than *sflimit are false exclusions.

  csm_model* testhyppescale = test_hypothesis_pe->scalesignal(sflimit);

  cslist.clear();
  TH1** pdarray = new TH1*[nchans];
  char *pdname;

  int* nobslist = new int[nbinstot*npx];
  Int_t* nobstotlist = new Int_t[npx];
  Int_t* nobsindex = new Int_t[npx];

  for (i=0;i<(Int_t) testhyppescale->channame.size(); i++)
    {
      pdname = new char[strlen(testhyppescale->channame[i])+strlen(" pseudodata ")];
      strcpy(pdname,testhyppescale->channame[i]);
      strcat(pdname," pseudodata");
      pdarray[i] = (TH1*) testhyppescale->chanmodel[i]->histotemplate[0]->Clone(pdname);
      delete [] pdname;
    }
  for (ipx=0;ipx<npx;ipx++)
    {
      testhyppescale->single_pseudoexperiment(pdarray);

      nobstotlist[ipx] = 0;
      ibin = 0;
      for (i=0;i<nchans;i++)
        {
          nbinsx = pdarray[i]->GetNbinsX();
          nbinsy = pdarray[i]->GetNbinsY();
          for (j=0;j<nbinsx;j++)
            {
              for (k=0;k<nbinsy;k++)
                {
                   if (nbinsy==1)
	             { nobslist[ibin+ipx*nbinstot] = (Int_t) nearbyint(pdarray[i]->GetBinContent(j+1)); }
	           else
	             { nobslist[ibin+ipx*nbinstot] = (Int_t) nearbyint(pdarray[i]->GetBinContent(j+1,k+1)); }
		   nobstotlist[ipx] += nobslist[ibin+ipx*nbinstot];
	           ibin++;
                }
            }
        } 
    }
  TMath::Sort(npx,nobstotlist,nobsindex,kTRUE);

  if (nglmax < nobstotlist[nobsindex[0]]/2 + 1)
    {
       nglmax = nobstotlist[nobsindex[0]]/2 + 1; 
       delete[] xgl;
       delete[] lwgl;
       xgl = new double[nglmax];
       lwgl = new double[nglmax];
    }

  ngl = 0;
  for (ipx=0;ipx<npx;ipx++)
    {
      if (ipx>0)
	{ if (nobstotlist[nobsindex[ipx]] != nobstotlist[nobsindex[ipx-1]])
	  { 
	    ngl = 0;
	  }
	}
      Double_t unc;
      Double_t p = (Double_t) cslimit(beta,nbinstot,nens,&(nobslist[nbinstot*nobsindex[ipx]]),ens,&ngl,xgl,lwgl,prior,&unc);
      if (pxprintflag)
	{
          cout << "Bayes Coverage Check px: " << p << endl;
	}
      cslist.push_back(p);
    }

  Int_t nfalse = 0;
  for (i = 0; i<npx; i++)
    {
      if (cslist[i] <= sflimit)
	{
	  nfalse++;
	}
    }

  *falsex = ((Double_t) nfalse)/((Double_t) npx);

  // clean up

  for (i=0;i<(Int_t) testhyppescale->channame.size(); i++)
    {
      delete pdarray[i];
    }
  delete[] pdarray;
  delete[] nobslist;
  delete[] nobsindex;
  delete[] nobstotlist;
  delete[] nobs;
  delete[] ens;
  delete[] xgl;
  delete[] lwgl;
  delete testhyppescale;
}

/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/
/*    genlimit Bayesian code from Joel                                     */
/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------*/

/*  Joel Heinrich  8 April 2005

Returns cross section posterior p.d.f. evaluated at s.

See CDF note 7587 for details:
http://www-cdf.fnal.gov/publications/cdf7587_genlimit.pdf

*/

void setdlcsn(int nchan, int nens, int nobs[], const EB* ens)
{
  dlcsn = 0;
  int set=0;
  double tmax=0;
  double tmin=0;
  double s=0;

  int i,k;
  double lgp = 0;
  const EB* p = ens;

  for(k=0;k<nchan;++k)
    lgp -= lgamma(nobs[k]+1);

  for(i=0;i<nens;++i) {
    double t = lgp, esum=0;
    for(k=0;k<nchan;++k) {
      const double mu = s*p->e + p->b;
      const int n = nobs[k];
      esum += p->e;
      //cout << "iens: " << i << " bin " << k << " nobs: " << n << " mc: " << mu << endl;
      t += ( (n>0) ? n*log(mu) : 0 ) - mu;
      ++p;
    }
   if (!set)
      { tmax = t; tmin = t; set = 1;}
    else
      { tmax=max(tmax,t); tmin = min(tmin,t); }
    dlcsn = -tmax;
 }
}


void setdlcsn2d(int nchan, int nens, int nobs[], const EB2D* ens)
{
  dlcsn2d = 0;
  int set=0;
  double tmax=0;
  double tmin=0;
  double s=0;

  int i,k;
  double lgp = 0;
  const EB2D* p = ens;

  for(k=0;k<nchan;++k)
    lgp -= lgamma(nobs[k]+1);

  for(i=0;i<nens;++i) {
    double t = lgp, esum=0;
    for(k=0;k<nchan;++k) {
      const double mu = s*(p->e1 + p->e2) + p->b;
      const int n = nobs[k];
      esum += (p->e1 + p->e2);
      //cout << "iens: " << i << " bin " << k << " nobs: " << n << " mc: " << mu << endl;
      t += ( (n>0) ? n*log(mu) : 0 ) - mu;
      ++p;
    }
   if (!set)
      { tmax = t; tmin = t; set = 1;}
    else
      { tmax=max(tmax,t); tmin = min(tmin,t); }
    dlcsn2d = -tmax;
 }
}

double cspdf(double s,double norm,
	     int nchan,int nens,const int nobs[],const EB* ens,PRIOR prior) {
  int i,k;
  double sum = 0, lgp = 0;
  const EB* p = ens;

  assert(nens>0);
  assert(prior==flat || prior==corr);

  for(k=0;k<nchan;++k)
    lgp -= lgamma(nobs[k]+1);

  for(i=0;i<nens;++i) {
    double t = lgp, esum=0;
    for(k=0;k<nchan;++k) {
      const double mu = s*p->e + p->b;
      const int n = nobs[k];
      esum += p->e;
      //cout << "iens: " << i << " bin " << k << " nobs: " << n << " mc: " << mu << endl;
      t += ( (n>0) ? n*log(mu) : 0 ) - mu;
      ++p;
    }
    sum += (prior==flat) ? exp(t+dlcsn) : esum*exp(t+dlcsn);
  }
  //cout << "pdf: " << s << " " << sum/nens << endl;
  return sum/(norm*nens);
}

// Tom Junk 21 April 2008 -- generalize to two signals

double cspdf2d(double s1, double s2, double norm,
	     int nchan,int nens,const int nobs[],const EB2D* ens,PRIOR prior) {
  int i,k;
  double sum = 0, lgp = 0;
  const EB2D* p = ens;

  assert(nens>0);
  assert(prior==flat || prior==corr);

  for(k=0;k<nchan;++k)
    lgp -= lgamma(nobs[k]+1);

  for(i=0;i<nens;++i) {
    double t = lgp, esum=0;
    for(k=0;k<nchan;++k) {
      const double mu = s1*p->e1+ s2*p->e2 + p->b;
      const int n = nobs[k];
      esum += (p->e1+p->e2);
      //cout << "iens: " << i << " bin " << k << " nobs: " << n << " mc: " << mu << endl;
      t += ( (n>0) ? n*log(mu) : 0 ) - mu;
      ++p;
    }
    sum += (prior==flat) ? exp(t+dlcsn2d) : esum*exp(t+dlcsn2d);
  }
  //cout << "pdf: " << s << " " << sum/nens << endl;
  return sum/(norm*nens);
}

/*  Joel Heinrich  8 April 2005

Function returns integral from xlo to infinity.

*uncertainty (if not null pointer) returned with uncertainty of
integral due to Monte Carlo statistical fluctuations of the prior
ensemble.

See CDF note 7587 for details:
http://www-cdf.fnal.gov/publications/cdf7587_genlimit.pdf


*/


double csint(double xlo,int nchan,int nens,const int nobs[],const EB* ens,
	     int* ngl,double xgl[],double lwgl[],
	     PRIOR prior,double* uncertainty) {
  int i,ntot=0;
  double sum=0, sum2=0, logscale=0;
  const EB* p=ens;

  assert(nens>0);
  assert(prior==flat || prior==corr);
  for(i=0;i<nchan;++i) {
    ntot += nobs[i];
    logscale -= lgamma(nobs[i]+1);
  }
  if(*ngl<=0)
    gausslaguerre(xgl,lwgl,*ngl=1+ntot/2,0.0);
  
  for(i=0;i<nens;++i) {
    const double t = csint0(xlo,logscale,nchan,nobs,p,*ngl,xgl,lwgl,prior);
    sum += t;
    sum2 += t*t;
    p+=nchan;
  }
  sum /= nens;
  sum2 /= nens;
  if(uncertainty)
    *uncertainty = (nens>1) ? sqrt((sum2-sum*sum)/(nens-1)) : 1.0 ;
  return sum;
}

double csint0(double xlo,double logscale,
		     int nchan,const int nobs[],const EB chan[],
		     int ngl,const double xgl[],const double lwgl[],
		     PRIOR prior) {
  int i,k;
  double sum=0, esum=0, bsum=0, resum;

  for(i=0;i<nchan;++i) {
    esum += chan[i].e;
    bsum += chan[i].b + xlo*chan[i].e;
  }
  assert(esum>0);
  resum=1/esum;

  for(k=0;k<ngl;++k) {
    const double xr = xgl[k]*resum + xlo;
    double t = logscale-bsum, v;
    for(i=0;i<nchan;++i)
      if(nobs[i]>0)
	t += nobs[i] * log( xr*chan[i].e + chan[i].b );
    if (dlcsn==0)
      {
	dlcsn = -t;
	if (dlcsn==0)
	  {
	    dlcsn = 1.0;
	  }
      }
    sum += v = exp(lwgl[k]+t+dlcsn);
    if(v<DBL_EPSILON*sum) break;
  }
  if(prior==flat)
    sum *= resum;
  return sum;
}


/*      Joel Heinrich  8 April 2005

    returns:
      *int1 = integral from xlo1 to xhi
      *int2 = integral from xlo2 to xhi
      *v11  = variance of *int1
      *v12  = covariance between *int1 and *int2
      *v22  = variance of *int2

See CDF note 7587 for details:
http://www-cdf.fnal.gov/publications/cdf7587_genlimit.pdf

*/



void csint2cut(double xlo1,double xlo2,double xhi,
	       int nchan,int nens,const int nobs[],const EB* ens,
	       int* ngl,double xgl[],double lwgl[],PRIOR prior,
	       double* int1,double* int2,
	       double* v11,double* v12,double* v22) {
  int i,ntot=0;
  double sum1=0, sum21=0, sum2=0, sum22=0, sump=0, logscale=0;
  const EB* p=ens;

  assert(nens>0);
  assert(prior==flat || prior==corr);
  for(i=0;i<nchan;++i) {
    ntot += nobs[i];
    logscale -= lgamma(nobs[i]+1);
  }
  if(*ngl<=0)
    gausslaguerre(xgl,lwgl,*ngl=1+ntot/2,0.0);
  
  for(i=0;i<nens;++i) {
    double t1=0, t2=0;
    csint02cut(xlo1,xlo2,xhi,logscale,nchan,nobs,p,*ngl,xgl,lwgl,
	       prior,&t1,&t2);
    sum1 += t1;
    sum21 += t1*t1;
    sum2 += t2;
    sum22 += t2*t2;
    sump += t1*t2;
    p+=nchan;
  }

  {
    const double rnens = 1.0/nens;
    sum1 *= rnens;
    sum21 *= rnens;
    sum2 *= rnens;
    sum22 *= rnens;
    sump *= rnens; 
    if(nens>1) {
      const double rn1 = 1.0/(nens-1);
      *v11 = (sum21-sum1*sum1)*rn1;
      *v22 = (sum22-sum2*sum2)*rn1;
      *v12 = (sump-sum1*sum2)*rn1;
    } else {
      *v11 = 1;
      *v22 = 1;
      *v12 = 0;
    }
  }

  *int1 = sum1;
  *int2 = sum2;
  return;
}

void csint02cut(double xlo1,double xlo2,double xhi,double logscale,
		    int nchan,const int nobs[],const EB chan[],
		    int ngl,const double xgl[],const double lwgl[],PRIOR prior,
		    double* int1,double* int2) {
  int i,k;
  double sum1=0, sum2=0, sum3=0, esum=0, bsum1=0, bsum2=0, bsum3=0, resum;

  for(i=0;i<nchan;++i) {
    const double ee=chan[i].e, bb=chan[i].b;
    esum += ee;
    bsum1 += bb + xlo1*ee;
    bsum2 += bb + xlo2*ee;
    bsum3 += bb + xhi*ee;
  }

  if(esum==0 && prior==corr) {
    *int1 = *int2 = 0;
    return;
  }

  assert(esum>0);
  resum=1/esum;

  for(k=0;k<ngl;++k) {
    const double xr = xgl[k]*resum + xlo1;
    double t1 = logscale-bsum1, v1;
    for(i=0;i<nchan;++i)
      if(nobs[i]>0)
	t1 += nobs[i] * log( xr*chan[i].e + chan[i].b );
    sum1 += v1 = exp(lwgl[k]+t1+dlcsn);
    if(v1<DBL_EPSILON*sum1) break;
  }

  for(k=0;k<ngl;++k) {
    const double xr = xgl[k]*resum + xlo2;
    double t2 = logscale-bsum2, v2;
    for(i=0;i<nchan;++i)
      if(nobs[i]>0)
	t2 += nobs[i] * log( xr*chan[i].e + chan[i].b );
    sum2 += v2 = exp(lwgl[k]+t2+dlcsn);
    if(v2<DBL_EPSILON*sum2) break;
  }

  for(k=0;k<ngl;++k) {
    const double xr = xgl[k]*resum + xhi;
    double t3 = logscale-bsum3, v3;
    for(i=0;i<nchan;++i)
      if(nobs[i]>0)
	t3 += nobs[i] * log( xr*chan[i].e + chan[i].b );
    sum3 += v3 = exp(lwgl[k]+t3+dlcsn);
    if(v3<DBL_EPSILON*sum3) break;
  }

  if(prior==flat) {
    *int1 = (sum1-sum3)*resum;
    *int2 = (sum2-sum3)*resum;
  } else {
    *int1 = sum1-sum3;
    *int2 = sum2-sum3;
  }
  return;
}


/*  Joel Heinrich  8 April 2005

Returns cross section upper limit.

*uncertainty (if not null pointer) returned with uncertainty of
limit due to Monte Carlo statistical fluctuations of the prior
ensemble.

See CDF note 7587 for details:
http://www-cdf.fnal.gov/publications/cdf7587_genlimit.pdf

*/


double cslimit(double beta,int nchan,int nens,const int nobs[],const EB* ens,
	       int* ngl,double xgl[],double lwgl[],
	       PRIOR prior,double* uncertainty) {
  const double eps=1.0e-6;

  
  /*
  cout << "nchan: " << nchan << endl;
  cout << " nens: " << nens << endl;
  cout << " prior: " << prior << endl;
  int i;
  for (i=0;i<nchan;i++)
    {
      if (nobs[i]>0)
	{  cout << "i, n, b, s: " << i << " " << nobs[i] << " " << ens[i].b << " " << ens[i].e << endl;
	//cout << "nobs(" << i << ") = " << nobs[i] << endl;
	}
    }
  */

  dlcsn = 0;

  double norm = csint(0,nchan,nens,nobs,ens,ngl,xgl,lwgl,prior,NULL);
  double limit = galim(beta,nchan,nens,nobs,ens);
  double dl=limit, rpdf=0;
  double lo=0, hi=1e200;

  while(fabs(dl)>1.0e-10*limit) {
    const double pbeta =
      1-csint(limit,nchan,nens,nobs,ens,ngl,xgl,lwgl,prior,NULL)/norm;
    rpdf = 1/cspdf(limit,norm,nchan,nens,nobs,ens,prior);
    if (pbeta > beta) {
      hi=limit*(1+eps);
    } 
    else {
      lo=limit*(1-eps);
    }
    dl = (pbeta-beta)*rpdf;
    if(limit-dl>=lo && limit-dl<=hi) {
      limit -= dl;
    } else {
      dl = limit - 0.5*(hi+lo);
      limit = 0.5*(hi+lo);
    }
  }

  if (uncertainty) {
    double i1=0, i2=0, v11=0, v22=0, v12=0;
    const double c = 1-beta;
    csint2(0,limit,nchan,nens,nobs,ens,ngl,xgl,lwgl,prior,
	   &i1,&i2,&v11,&v12,&v22);
    *uncertainty = rpdf*sqrt(v22 + v11*c*c - 2*v12*c)/norm   ;
  }

  return limit;
}


/*  Joel Heinrich  8 April 2005

Returns cross section upper limit.

*uncertainty (if not null pointer) returned with uncertainty of
limit due to Monte Carlo statistical fluctuations of the prior
ensemble.

See CDF note 7587 for details:
http://www-cdf.fnal.gov/publications/cdf7587_genlimit.pdf

*/




double cscutlimit(double beta,double smax,
		  int nchan,int nens,const int nobs[],const EB* ens,
		  int* ngl,double xgl[],double lwgl[],
		  PRIOR prior,double* uncertainty) {
  const double norm = csint(0,nchan,nens,nobs,ens,ngl,xgl,lwgl,prior,NULL);
  const double tail = csint(smax,nchan,nens,nobs,ens,ngl,xgl,lwgl,prior,NULL);
  const double eps=1.0e-6;
  double limit = galim(beta,nchan,nens,nobs,ens), rpdf=0;
  double dl=limit, lo=0, hi=smax;

  if(beta<=0) return 0;
  if(beta>=1) return smax;

  if(limit>smax || limit<0) dl = limit = 0.5*smax;

  while(fabs(dl)>1.0e-10*limit) {
    double pbeta =
      1-(csint(limit,nchan,nens,nobs,ens,ngl,xgl,lwgl,prior,NULL)-tail)/
      (norm-tail);
    rpdf = 1/cspdf(limit,norm-tail,nchan,nens,nobs,ens,prior);
    if(pbeta>beta) {
      hi=limit*(1+eps);
    } else {
      lo=limit*(1-eps);
    }
    dl = (pbeta-beta)*rpdf;
    if(limit-dl>=lo && limit-dl<=hi) {
      limit -= dl;
    } else {
      dl = limit - 0.5*(hi+lo);
      limit = 0.5*(hi+lo);
    }

  }


  if (uncertainty) {
    double i1=0, i2=0, v11=0, v22=0, v12=0;
    const double c = 1-beta;

    csint2cut(0,limit,smax,nchan,nens,nobs,ens,ngl,xgl,lwgl,prior,
	      &i1,&i2,&v11,&v12,&v22);
    *uncertainty = rpdf*sqrt(v22 + v11*c*c - 2*v12*c)/(norm-tail)   ;
  }


  return limit;
}


/*         Joel Heinrich  24 March 2005

   returns crude Gaussian approximation to upper limit.
   For use as a starting point.
*/

double galim(double beta,int nchan,int nens,const int nobs[],const EB* ens) {
  double mean=0,sigma=0;
  gameansigma(&mean,&sigma,nchan,nens,nobs,ens);
  return mean-sigma*arcfreq( (1-freq(-mean/sigma))*(1-beta) );
}



void gameansigma(double *mean,double *sigma,
		 int nchan,int nens,const int nobs[],const EB* ens) {

  double sum=0,sum2=0,vsum=0;
  const EB* p = ens;
  int i,j;
  for(i=0;i<nens;++i) {
    double s=0,s2=0;
    for(j=0;j<nchan;++j) {
      const int n = nobs[j];
      const double eps = p->e;
      s += (n-p->b)*eps/(n+1);
      s2 += eps*eps/(n+1);
      ++p;
    }
    s /= s2;
    vsum += 1/s2;
    sum += s;
    sum2 += s*s;
  }
  
  *mean = sum/nens;
  *sigma = sqrt(vsum/nens + sum2/nens - (*mean)*(*mean));
  return;
}

#define rdfreq(x) (exp(0.5*(x)*(x))*2.50662827463100050242)

#define C0 2.515517
#define C1 0.802853
#define C2 0.010328
#define D0 1.0
#define D1 1.432788
#define D2 0.189269
#define D3 0.001308

double arcfreq(double y) {
  const double yy = (y>0.5) ? 1-y : y, t = sqrt(-2*log(yy));
  double x = (C0+t*(C1+t*C2))/(D0+t*(D1+t*(D2+t*D3))) - t;
  x -= (freq(x) - yy)*rdfreq(x);
  x -= (freq(x) - yy)*rdfreq(x);
  return (y>0.5) ? -x : x;
}


/*

   Joel Heinrich
   February 10 2005

Returns Gauss-Laguerre quadrature abscissas and log(weights) which can
be used to approximate

      integral u=0 to infinity pow(u,alpha)*exp(-u)*f(u) du
as
      sum k=0 to n-1  exp(lw[k])*f(x[k])

or equivalently

      sum k=0 to n-1  exp(lw[k]+log(f(x[k])))

The quadrature is exact for polynomial f of degree 2n-1 or less.

*/


void gausslaguerre(double x[],double lw[],int n,double alpha){
  const int nshift = 20;
  const double shift = 1<<nshift, rshift=1/shift;
  int i;
  double z=0;
  
  for(i=0;i<n;++i) {
    int j=0, k=2, nscale=0;
    double dz=0.0, p1=0, p2=0;
    if(i==0) {
      z=(1.0+alpha)*(3.0+0.92*alpha)/(1.0+2.4*n+1.8*alpha);
    } else if(i==1) {
      z += (15.0+6.25*alpha)/(1.0+2.5*n+0.9*alpha);
    } else if(i==2) {
      const double ai=i-1;
      z += ( (1.0+2.55*ai)/(1.9*ai) + 1.26*ai*alpha/(1.0+3.5*ai) )*
	(z-x[i-2])/(1.0+0.3*alpha);
    } else if(i==3) {
      z = 3.0*(x[2]-x[1])+x[0];
    } else if(i==4) {
      z = 4.0*x[3] - 6.0*x[2] + 4.0*x[1] - x[0];
    } else if(i==5) {
      z = 5.0*x[4] - 10.0*x[3] + 10.0*x[2] - 5.0*x[1] + x[0];
    } else {
      z = 6.0*x[i-1] - 15.0*x[i-2] + 20.0*x[i-3] -
	15.0*x[i-4] + 6.0*x[i-5] - x[i-6];
    }
    while(k>0) {
      p1=1;
      p2=0;
      nscale=0;
      z -= dz;
      for(j=1;j<=n;++j){
	const double p3=p2;
	p2=p1;
	p1=((2*j-1+alpha-z)*p2 - (j-1+alpha)*p3)/j;
	if(fabs(p2)>shift) {
	  ++nscale;
	  p1 *= rshift;
	  p2 *= rshift;
	}
      }
      dz = p1*z/(n*p1-(n+alpha)*p2);
      if(fabs(dz)<1.0e-10*z)--k;
    }
    x[i]=z;
    lw[i] = log(z/(p2*p2)) - 2*nshift*nscale*M_LN2 ;
  }
  
  {
    double t = 0.0;
    for(i=n-1;i>=0;--i)
      t += exp(lw[i]);
    t = lgamma(alpha+1)-log(t);
    for(i=0;i<n;++i)
      lw[i] += t;
  }

  return;
}

/*    Joel Heinrich  8 April 2005

   returns:
      *int1 = integral from xlo1 to infinity
      *int2 = integral from xlo2 to infinity
      *v11  = variance of *int1
      *v12  = covariance between *int1 and *int2
      *v22  = variance of *int2


See CDF note 7587 for details:
http://www-cdf.fnal.gov/publications/cdf7587_genlimit.pdf


*/



void csint2(double xlo1,double xlo2,
	    int nchan,int nens,const int nobs[],const EB* ens,
	    int* ngl,double xgl[],double lwgl[],PRIOR prior,
	    double* int1,double* int2,
	    double* v11,double* v12,double* v22) {
  int i,ntot=0;
  double sum1=0, sum21=0, sum2=0, sum22=0, sump=0, logscale=0;
  const EB* p=ens;

  assert(nens>0);
  assert(prior==flat || prior==corr);
  for(i=0;i<nchan;++i) {
    ntot += nobs[i];
    logscale -= lgamma(nobs[i]+1);
  }
  if(*ngl<=0)
    gausslaguerre(xgl,lwgl,*ngl=1+ntot/2,0.0);
  
  for(i=0;i<nens;++i) {
    double t1=0, t2=0;
    csint02(xlo1,xlo2,logscale,nchan,nobs,p,*ngl,xgl,lwgl,prior,&t1,&t2);
    sum1 += t1;
    sum21 += t1*t1;
    sum2 += t2;
    sum22 += t2*t2;
    sump += t1*t2;
    p+=nchan;
  }

  {
    const double rnens = 1.0/nens;
    sum1 *= rnens;
    sum21 *= rnens;
    sum2 *= rnens;
    sum22 *= rnens;
    sump *= rnens; 
    if(nens>1) {
      const double rn1 = 1.0/(nens-1);
      *v11 = (sum21-sum1*sum1)*rn1;
      *v22 = (sum22-sum2*sum2)*rn1;
      *v12 = (sump-sum1*sum2)*rn1;
    } else {
      *v11 = 1;
      *v22 = 1;
      *v12 = 0;
    }
  }

  *int1 = sum1;
  *int2 = sum2;
  return;
}

void csint02(double xlo1,double xlo2,double logscale,
		    int nchan,const int nobs[],const EB chan[],
		    int ngl,const double xgl[],const double lwgl[],PRIOR prior,
		    double* int1,double* int2) {
  int i,k;
  double sum1=0, sum2=0, esum=0, bsum1=0, bsum2=0, resum;

  for(i=0;i<nchan;++i) {
    esum += chan[i].e;
    bsum1 += chan[i].b + xlo1*chan[i].e;
    bsum2 += chan[i].b + xlo2*chan[i].e;
  }
  assert(esum>0);
  resum=1/esum;

  for(k=0;k<ngl;++k) {
    const double xr = xgl[k]*resum + xlo1;
    double t1 = logscale-bsum1, v1;
    for(i=0;i<nchan;++i)
      if(nobs[i]>0)
	t1 += nobs[i] * log( xr*chan[i].e + chan[i].b );
    sum1 += v1 = exp(lwgl[k]+t1+dlcsn);
    if(v1<DBL_EPSILON*sum1) break;
  }

  for(k=0;k<ngl;++k) {
    const double xr = xgl[k]*resum + xlo2;
    double t2 = logscale-bsum2, v2;
    for(i=0;i<nchan;++i)
      if(nobs[i]>0)
	t2 += nobs[i] * log( xr*chan[i].e + chan[i].b );
    sum2 += v2 = exp(lwgl[k]+t2+dlcsn);
    if(v2<DBL_EPSILON*sum2) break;
  }

  if(prior==flat) {
    sum1 *=resum;
    sum2 *=resum;
  }

  *int1 = sum1;
  *int2 = sum2;
  return;
}

// access the globals in here with accessor methods.

double mclimit_csm::getdlcsn() { return dlcsn; }
double mclimit_csm::getdlcsn2d() { return dlcsn2d; }
double mclimit_csm::getbhnorm() { return bhnorm; }
