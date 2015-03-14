///////////////////////////////////////////////////////////////////////////////
// mclimit_csm -- use a chisquared minimized over nuisance parameters
// to compute cls.

// Class to run the TMinuit minimization of T. Devlin's chisquared
// defined in CDF 3126, minimized over the nuisance parameters.

// version dated January 14, 2009
// Author:  Tom Junk, Fermilab.  trj@fnal.gov

#ifndef CSM_H
#define CSM_H

#include "TH1.h"
#include "TRandom3.h"

#include <vector>

#include "csm_channel_model.h"
#include "csm.h"
using std::vector;

// one-sided or two-sided 3-sigma or 5-sigma

#define MCLIMIT_CSM_TWOSIDED

// what to do about 1-sided or 2-sided 2-sigmas?

#ifdef MCLIMIT_CSM_TWOSIDED
#define MCLIMIT_CSM_2S 0.02275
#define MCLIMIT_CSM_3S 1.349898E-3
#define MCLIMIT_CSM_5S 2.866516E-7
#else
#define MCLIMIT_CSM_2S 0.0455
#define MCLIMIT_CSM_3S 2.6998E-3
#define MCLIMIT_CSM_5S 5.7330E-7
#endif

// cumulative probabilities for defining bands on test statistic
// and CL plots

#define MCLIMIT_CSM_MCLM2S 0.02275
#define MCLIMIT_CSM_MCLM1S 0.16
#define MCLIMIT_CSM_MCLMED 0.5
#define MCLIMIT_CSM_MCLP1S 0.84
#define MCLIMIT_CSM_MCLP2S 0.97725

// some messages to pass around inside for the s95 calculator

#define MCLIMIT_CSM_CLS 1
#define MCLIMIT_CSM_CLSM2 2
#define MCLIMIT_CSM_CLSM1 3
#define MCLIMIT_CSM_CLSMED 4
#define MCLIMIT_CSM_CLSP1 5
#define MCLIMIT_CSM_CLSP2 6

#define MCLIMIT_CSM_LUMI95 1
#define MCLIMIT_CSM_LUMI3S 2
#define MCLIMIT_CSM_LUMI5S 3

// use to steer the choice of integration method in bayes_heinrich and bayes_heinrich_withexpect

#define CSM_BAYESINTEGRAL_JOEL 0
#define CSM_BAYESINTEGRAL_QUICK 1

#define PVMORPH_MAXBINS 5000
#define CSM_DEBUGPRINT 1

class mclimit_csm: public TObject
{
 public:
   mclimit_csm();
   ~mclimit_csm();

   void print_version();

   void set_datahist(TH1 *,const char *); /* data histogram and channel name */
   // the hypothesis set routines do not make clones of their inputs, they just
   // store pointers to the models. Best practice -- fully define a model (that
   // is, add all templates, before calling these routines and do not update
   // the models before the pseudoexperiments are run.

   void set_null_hypothesis(csm_model *);
   void set_test_hypothesis(csm_model *);
   void set_null_hypothesis_pe(csm_model *);
   void set_test_hypothesis_pe(csm_model *);
   void set_npe(Int_t);      // sets the number of pseudoexperiments to do.  
                             // The default is set in the constructor to 10000
   Int_t get_npe();          // returns the value set in set_npe

   void set_chisquarehistos(TH1 *,TH1 *,TH1 *,TH1 *);
   // set pointers to histograms to accumulate chisquare distributions for
   // null hyp chisquare distrib in null hyp pseudoexperiments
   // test hyp chisquare distrib in null hyp pseudoexperiments
   // null hyp chisquare distrib in test hyp pseudoexperiments
   // test hyp chisquare distrib in test hyp pseudoexperiments

   // these accessor methods are used to control the behavior of MINUIT in the calculation
   // of the test statistic.  Their values are just passed to the csm class instance when
   // the test statistic is computed.

   void  setminuitmaxcalls(Int_t); // max number of calls to MINUIT's minimization function
                                   // default: 500
   Int_t getminuitmaxcalls();
   void  setminosmaxcalls (Int_t); // max number of calls MINOS is allowed to use per parameter
                                   // default: 500
   Int_t getminosmaxcalls ();
   void  setminuitstepsize(Double_t);  // initial step size for MINUIT parameters (default=0.1)
   Double_t getminuitstepsize();
   void setminosflag(bool);   // true: call MINOS. Best to have printing set too if 
                              // you're going to run MINOS.  False:  Do not call MINOS (default) 
   bool getminosflag(); 
   void setprintflag(bool);  // true: let MINUIT print stuff out; FALSE: turn off MINUIT printing (default)
   bool getprintflag();
   void setpxprintflag(bool);    // print out pseudoexperiment results
                             // this flag affects printing test statistic outputs in the
                             // run_pseudoexperiment loop, as well as bayes_heinrich_withexpect
                             // and bayes_heinrich_coverage check
   bool getpxprintflag();
   void setprintnpextremeflag(bool); // turn on printing of nuisance parameters if the null hyp has
                                     // a very signal-like outcome.  True: print, False: don't print.  Default: F
   void setprintextremevalue(Double_t);  // cut on how negative -2lnQ(bg) has to be before triggering
                                         // a nuisance parameter printout

   void run_pseudoexperiments(); // see set_npe() above to determine how many pseudoexperiments to run

   // output retrieval:

   Double_t cls  ();
   Double_t clsw ();   // cls+b computed with null hyp px's reweighted -- good for small cls's
   Double_t clsb ();
   Double_t clsbw();  // cls+b computed with null hyp px's reweighted -- good for small cls's
   Double_t clb  ();
   Double_t omclb();  // "1-clb" computed as a p-value, including the probability of the exact outcome
                      // observed in the data.  Computed with null hypothesis pseudoexperiments
   Double_t omclbw(); // Same as above, but using test hypothesis pseudoexperiments, reweighted
   // with the likelihood ratio to approximate the null hypothesis distribution

   Double_t ts();  // test statistic -- a delta chisquared -- computed for the observed data

   Double_t tsbm2(); // distributions of test statistic in null hyp pseudoexperiments 2 sigma low edge
   Double_t tsbm1(); // 1 sigma low edge
   Double_t tsbmed(); // median test statistic in null hyp pseudoexperiments
   Double_t tsbp1();  // 1 sigma upper edge
   Double_t tsbp2();  // 2 sigma upper edge

   Double_t tssm2();  // distributions of test statistic in test hyp pseudoexperiments 2 sigma low edge
   Double_t tssm1();  // 1 sigma low edge
   Double_t tssmed(); // median test statistic in null hyp pseudoexperiments
   Double_t tssp1();  // 1 sigma upper edge
   Double_t tssp2();  // 2 sigma upper edge
 
   Double_t clsexpbm2();  // Expected cls in null hypothesis -- 2 sigma low edge
   Double_t clsexpbm1();  // Expected cls in null hypothesis -- 2 sigma low edge
   Double_t clsexpbmed(); // Expected cls in null hypothesis -- median
   Double_t clsexpbp1();  // Expected cls in null hypothesis -- 1 sigma upper edge
   Double_t clsexpbp2();  // Expected cls in null hypothesis -- 2 sigma upper edge

   // generalize the above a bit in case CLs is not a monotonic function of -2lnQ

   Double_t gclsexpbm2(); // Expected cls in null hypothesis -- 2 sigma low edge
   Double_t gclsexpbm1(); // Expected cls in null hypothesis -- 2 sigma low edge
   Double_t gclsexpbmed(); // Expected cls in null hypothesis -- median
   Double_t gclsexpbp1(); // Expected cls in null hypothesis -- 1 sigma upper edge
   Double_t gclsexpbp2(); // Expected cls in null hypothesis -- 2 sigma upper edge

   Double_t clsbexpbm2();  // Expected clsb in null hypothesis -- 2 sigma low edge
   Double_t clsbexpbm1();  // Expected clsb in null hypothesis -- 2 sigma low edge
   Double_t clsbexpbmed(); // Expected clsb in null hypothesis -- median
   Double_t clsbexpbp1();  // Expected clsb in null hypothesis -- 1 sigma upper edge
   Double_t clsbexpbp2();  // Expected clsb in null hypothesis -- 2 sigma upper edge


   // computed with null hypothesis px's reweighted -- good for small expected CLs's

   Double_t clsexpbm2w(); // Expected cls in null hypothesis -- 2 sigma low edge
   Double_t clsexpbm1w(); // Expected cls in null hypothesis -- 2 sigma low edge
   Double_t clsexpbmedw();// Expected cls in null hypothesis -- median
   Double_t clsexpbp1w(); // Expected cls in null hypothesis -- 1 sigma upper edge
   Double_t clsexpbp2w(); // Expected cls in null hypothesis -- 2 sigma upper edge

   Double_t clsexpsm2();  // Expected cls in test hypothesis -- 2 sigma low edge
   Double_t clsexpsm1();  // Expected cls in test hypothesis -- 2 sigma low edge
   Double_t clsexpsmed(); // Expected cls in test hypothesis -- median
   Double_t clsexpsp1();  // Expected cls in test hypothesis -- 1 sigma upper edge
   Double_t clsexpsp2();  // Expected cls in test hypothesis -- 2 sigma upper edge
 
   // these accessors below use the CLs definition of CLb which includes the 
   // probability of observing exactly the data outcome
   // (subtracting it from 1 makes 1-CLb computed with these routines omit the
   // probability of observing exactly the data outcome)  Not to be used
   // for discovery significance!

   Double_t clbexpsm2(); // Expected clb in test hypothesis -- 2 sigma low edge
   Double_t clbexpsm1(); // Expected clb in test hypothesis -- 2 sigma low edge
   Double_t clbexpsmed();// Expected clb in test hypothesis -- median
   Double_t clbexpsp1(); // Expected clb in test hypothesis -- 1 sigma upper edge
   Double_t clbexpsp2(); // Expected clb in test hypothesis -- 2 sigma upper edge

   // these accessors below use the p-value definition of 1-CLb which includes the
   // probability of observing exactly the data outcome.  These are computed
   // using null hypothesis px's to compute 1-CLb's.

   Double_t omclbexpsm2(); // Expected clb in test hypothesis -- 2 sigma low edge
   Double_t omclbexpsm1(); // Expected clb in test hypothesis -- 2 sigma low edge
   Double_t omclbexpsmed(); // Expected clb in test hypothesis -- median
   Double_t omclbexpsp1(); // Expected clb in test hypothesis -- 1 sigma upper edge
   Double_t omclbexpsp2(); // Expected clb in test hypothesis -- 2 sigma upper edge

   // Same as above, but use the test hypothesis pseudoexperiments, reweighted with
   // the likelihood ratio, to approximate the distribution of the null hypothesis
   // distribution of the test statistic
   Double_t omclbexpsm2w(); // Expected clb in test hypothesis -- 2 sigma low edge
   Double_t omclbexpsm1w(); // Expected clb in test hypothesis -- 2 sigma low edge
   Double_t omclbexpsmedw(); // Expected clb in test hypothesis -- median
   Double_t omclbexpsp1w(); // Expected clb in test hypothesis -- 1 sigma upper edge
   Double_t omclbexpsp2w(); // Expected clb in test hypothesis -- 2 sigma upper edge

   // these accessors below use the CLs definition of CLb which includes the 
   // probability of observing exactly the data outcome
   // (subtracting it from 1 makes 1-CLb computed with these routines omit the
   // probability of observing exactly the data outcome)
   Double_t clbexpbm2(); // Expected clb in null hypothesis -- 2 sigma low edge
   Double_t clbexpbm1(); // Expected clb in null hypothesis -- 2 sigma low edge
   Double_t clbexpbmed(); // Expected clb in null hypothesis -- median
   Double_t clbexpbp1(); // Expected clb in null hypothesis -- 1 sigma upper edge
   Double_t clbexpbp2(); // Expected clb in null hypothesis -- 2 sigma upper edge

   // these accessors below use the p-value definition of 1-CLb which includes the
   // probability of observing exactly the data outcome
   Double_t omclbexpbm2(); // Expected 1-clb in null hypothesis -- 2 sigma low edge
   Double_t omclbexpbm1(); // Expected 1-clb in null hypothesis -- 2 sigma low edge
   Double_t omclbexpbmed(); // Expected 1-clb in null hypothesis -- median
   Double_t omclbexpbp1(); // Expected 1-clb in null hypothesis -- 1 sigma upper edge
   Double_t omclbexpbp2(); // Expected 1-clb in null hypothesis -- 2 sigma upper edge

   // Same as above, but use the test hypothesis pseudoexperiments, reweighted with
   // the likelihood ratio, to approximate the distribution of the null hypothesis
   // distribution of the test statistic
   Double_t omclbexpbm2w(); // Expected 1-clb in null hypothesis -- 2 sigma low edge
   Double_t omclbexpbm1w(); // Expected 1-clb in null hypothesis -- 2 sigma low edge
   Double_t omclbexpbmedw(); // Expected 1-clb in null hypothesis -- median
   Double_t omclbexpbp1w(); // Expected 1-clb in null hypothesis -- 1 sigma upper edge
   Double_t omclbexpbp2w(); // Expected 1-clb in null hypothesis -- 2 sigma upper edge

   // these three below are computed with test hyp. px's reweighted to get the small
   // tails of the null hyp distribution better modeled with low px stats
   Double_t p2sigmat(); // Probability of a 2-sigma evidence assuming test hyp. is true
   Double_t p3sigmat(); // Probability of a 3-sigma evidence assuming test hyp. is true
   Double_t p5sigmat(); // Probability of a 5-sigma discovery assuming test hyp. is true

   // these are computed with the null hyp px's -- can be unreliable for low px stats.
   // hopefully the 1-CLb p-value is uniformly distributed between 0 and 1, but
   // these routines, as well as omclbexpb* are designed to quantify just how
   // uniform that is.
   Double_t p2sigman(); // Probability of a 2-sigma evidence assuming null hyp. is true
   Double_t p3sigman(); // Probability of a 3-sigma evidence assuming null hyp. is true
   Double_t p5sigman(); // Probability of a 5-sigma discovery assuming null hyp. is true

   Double_t calc_chi2(csm_model *, const TH1* []);  // interface to chisquare calculator
                                                   // using our models

   Double_t weightratio(csm_model *testhyp, csm_model *nullhyp, TH1 *[]); // weight factor
   // for reweighting MC px's using the varied templates.

   // rate limit calculators
   Double_t s95();     // scale factor on signal which is excluded at exactly 95% CL
   Double_t s95m2();   // variation around the median expected s95 in the bg hypothesis, -2 sigma
   Double_t s95m1();   // variation around the median expected s95 in the bg hypothesis, -1 sigma
   Double_t s95med();  // median expected s95 in the background hypothesis
   Double_t s95p1();   // variation around the median expected s95 in the bg hypothesis, +1 sigma
   Double_t s95p2();   // variation around the median expected s95 in the bg hypothesis, +2 sigma

   Double_t lumi95(); // calculates the lumi needed for a median experiment to exclude at 95% CL
                      // what's returned is a multiplicative factor on whatever luminosity was used
                      // to construct the test and null hypotheses.
   Double_t lumi3s(); // calculates the lumi needed for a median experiment to discover at 3 sigma
   Double_t lumi5s(); // calculates the lumi needed for a median experiment to discover at 5 sigma

   void tshists(TH1*,TH1*);  // fills histograms with test statisic values in the pseudoexperiments
                             // (you define the binning).  First histo: test hypothesis, second histo:
                             // null hypothesis

   void plotlnsb(TH1 *mcb_hist, TH1 *mcs_hist, TH1 *data_hist); // make a summed plot of ln(1+s/b) given the 
                             // input histogram pointers.  They are filled in with summed MC 
                             // and data (they are reset first).

   // Call Joel Heinrich's (CDF 7587) genlimit
   // Bayesian limit calculator.  First arg:  credibility level:  e.g., 0.95.  Second arg, 
   // scale factor on signal to produce the limit.  Third arg, uncertainty on the limit scale factor.
   // as with the s95 routines above, this assumes that the signal adds incoherently to the
   // background.  Requires set_test_hypothesis_pe to be called first in order to make the
   // "Prior ensemble".  The size of the prior ensemble is set with set_npe()
   // also call the relevant set_datahist methods too.

   void bayes_heinrich(Double_t beta, Double_t* sflimit, Double_t* unc);

   // Routine to call Joel Heinrich's (CDF 7587) genlimit, a Bayesian limit calculator,
   // but to repeat the calculation for pseudoexpeirments drawn from the null hypothesis
   // (be sure to call both set_test_hypothesis_pe and set_null_hypothesis_pe before using
   // this).  This computes the observed and expected limits.
   // Arguments:  1:  beta (credibility level, e.g. 0.95)
   // Argument 2:  observed limit
   // Agrument 3:  error on observed limit
   // Argument 4:  npx pseudoexperiments to run to compute expected limits
   // Arguments 5-9: Expected limits.  Median +-1, 2 sigma expectations

   void bayes_heinrich_withexpect(Double_t  beta,     // beta (credibility level, e.g. 0.95)
				  Double_t* sflimit,  // observed limit (scale factor)
				  Double_t* unc,      // error on observed limit
				  Int_t     npx,      // N pseudoexperiments to run 
                                  Double_t* sm2,      // median-2*sigma
				  Double_t* sm1,      // median-sigma
				  Double_t* smed,     // median expected limit
				  Double_t* sp1,      // median+sigma
                                  Double_t* sp2);     // median +2*sigma limit

   // Some extra things that bayes_heinrich and bayes_heinrich_withexpect will compute,
   // if requested.  You can get a plot of the posterior PDF by specifying the range
   // over which it is to be evaluated and the point sample density  -- These are initialized
   // to zero by the constructor.  Just specify the beginning and the end of the interval and
   // the step size, and the bayes_posterior vector will be filled in when bayes_heinrich and 
   // bayes_heinrich_withexpect are called.  bayes_interval_end > bayes_interval_begin and
   // bayes_inteval_step > 0 for bayes_posterior to be filled in.

   Double_t bayes_interval_begin;
   Double_t bayes_interval_end;
   Double_t bayes_interval_step;

   // these two vectors are filled if the three paramters above are specified
   // and bayes_heinrich or bayes_heinrich_withexpect is called.  Plot the
   // contents of the first vector versus the values of the second vector.

   vector<Double_t> bayes_posterior;
   vector<Double_t> bayes_posterior_points;

   // Default method is Joel's (CDF 7587), which is option 0 (CSM_BAYESINTEGRAL_JOEL)
   // Option 1 (CSM_BAYESINTEGRAL_QUICK) uses the interval above and performs a quick and
   // dirty numerical sum integral to speed up slow calculations.

   void set_bayes_integration_method(int imethod);
   int get_bayes_integration_method();

   // With bayes_heinrich_withexpect, you also can get a histogram of expected limits on the
   // background-only pseudoexperiments.  Just specify the histogram to receive the entries
   // (it is reset on entry into bayes_heinrich_withexpect and filled in) -- the pointer is
   // initially null, and no filling occurs unless this pointer is set.


   TH1F *bayes_pseudoexperiment_limits;

   void bayes_heinrich_coverage_check(Double_t beta,
                                      Double_t sflimit,
				      Int_t npx,
                                      Double_t* falsex);  // checks the false exclusion rate (coverage)

// accessor methods (rather cobbled onto a method of getting access to a 
// value which is global to mclimit_csm.C, since the Bayesian routines are 
// written in C and not C++, this value wasn't put in a nicer place.

// dlcsn is an exponential factor on the posterior.  If you want the sum of the likelihoods,
// you need to take the posterior and multiply it by exp(-dlcsn)/bnorm, where dlcsn is the
// output of getdlcsn, and bnorm is the output of getbnorm

    double getdlcsn();

    double getdlcsn2d();

    double getbhnorm();


   //---------------------------------------------------------------

   // fit a cross section using Joel Heinrich's routines  -- only testhyp_pe is used here.
   // cross sections are fit as scale factors of the cross section put in the templates
   // uses the data histograms supplied.  This also fills in the vectors
   // bayes_posterior and bayees_posterior_points
   // for now it needs a suggestion of what interval the maximum is expected
   // (to set the scale), so set bayes_interval_begin and bayes_interval_end and bayes_interval_step
   // make bayes_interval_step very small (otherwise there will be some discretization in your answer)

   void bh_xsfit(Double_t *xsfit, Double_t *downerr, Double_t *uperr);

   // draw pseudoexperiments from the fluctuated testhyp and find the cross section
   // from each one using testhyp_pe  Get expected values.

   void bh_xsfit_expect(Int_t npx, Double_t *xsfitavg, Double_t *m2s, 
			Double_t *m1s, Double_t *med, Double_t *p1s, Double_t *p2s);

   // scan over two signals and print out the marginalized posterior for later analysis
   // always assume the two signals are declared in the same order in all channels.
   // If more than two signals are present, the first one in  each channel is called signal 1,
   // and the sum of all others is called signal 2

   void bh_2d_scan(Double_t s1low, Double_t s1high, Double_t ds1,
                   Double_t s2low, Double_t s2high, Double_t ds2);

// This routine below runs pseudoexperiments to get the expected distributions.  It does all the
// calculations that bh_2d_scan does for the real data, but for the pseudoexperiments instead, and prints
// out the 2D posterior distributions, cumulative integrals, and best-fit points for each pseudoexperiment.
// As an added feature, you put in the signal scale factors used in the pseudoexperiment generation
//  (testhyp is used for pseudoexperiments, and testhyp_pe -- yes, I know that seems backwards -- is the
//  model used to compute the posterior), relative to the model used to test as s1true and s2true, and it
// will compute 68% and 95% coverage fractions based on the nearest grid point.  Verbose printout is
// supplied.

   void bh_2d_scan_expect(Int_t npx, Double_t s1true, Double_t s2true,
			  Double_t s1low, Double_t s1high, Double_t ds1,
			  Double_t s2low, Double_t s2high, Double_t ds2);

   TH1* get_datahist(Int_t);
                                        // 2010-01-27 P.Murat
   TH1* get_datahist(const char* ChannelName);

   void printsbd();              // dump all signal, background, candidates

 private:

   csm_model *null_hypothesis;
   csm_model *test_hypothesis;
   csm_model *null_hypothesis_pe;
   csm_model *test_hypothesis_pe;

   TRandom3*   fBayesRandom;

   Int_t nmc;  // number of pseudoexperiments which have been done
               // this is set to zero by the constructor and by
               // anything that modifies the hypotheses, which is the
               // indication that the pseudoexperiments need to be redone.

   Int_t nmc_req;  // number of pseudoexperiments to do

   Int_t recalctsflag; // 1 if we need to redo the data test statistic

   Double_t *tsb;  // test statistic in null hypothesis -- one per pseudoexperiment
   Double_t *tss;  // test statistic in test hypothesis -- one per pseudoexperiment
   Double_t *wtsb;  // weight for converting the corresponding tsb into a s+b px probability
   Double_t *wtss;  // weight for converting the corresponding tss into a b px probability
   Int_t *itsb; // sort order for tsb
   Int_t *itss; // sort order for tss
   Double_t tsd;                     // test statistic for data
   vector<TH1*> datahist;
   vector<char*> dhname;            // channel names for each data histogram -- must match
                                    // the channel names used in the models.

   Double_t s95aux(Int_t); // s95 calculator using a function you pass in.
   Double_t lumipaux(Int_t); // luminosity threshold calculator
   Double_t clsaux(Double_t); // cls, clsb and clb for an arbitrary test statistic
   Double_t clsauxw(Double_t); // cls, clsb and clb for an arbitrary test statistic
   Double_t clsbaux(Double_t);
   Double_t clsbauxw(Double_t);
   Double_t clbaux(Double_t);
   Double_t omclbaux(Double_t);
   Double_t omclbauxw(Double_t);
   Double_t gclsaux(Double_t);

   TH1 *nullnullchisquare;
   TH1 *nulltestchisquare;
   TH1 *testnullchisquare;
   TH1 *testtestchisquare;

   Int_t minuitmaxcalls;
   Int_t minosmaxcalls;
   Double_t minuitstepsize;
   bool minuitprintflag;       // true if we let MINUIT print stuff out.
   bool minosflag;             // true if we are to call MINOS
   bool pxprintflag;           // print out results of pseudoexperiments
   bool extremenpprintflag;    // print out nuisance parameters on extreme null hypothesis pseudoexperiments
   Double_t extremenpprintvalue;  // cut on how negative -2lnQ(bg) has to be before triggering a printout

   int bayesintegralmethod;
   Double_t quickbint(Double_t); // quick and dirty integrator

   ClassDef(mclimit_csm,0)
};

#endif

#ifndef GENLIMIT

typedef struct {
  float e,b;
}EB;

typedef struct {
  float e1,e2,b;
}EB2D;

typedef enum {
  flat=10,
  corr=20
} PRIOR;

double cspdf(double s,double norm,
             int nchan,int nens,const int nobs[],const EB* ens,PRIOR prior);

double cspdf2d(double s1, double s2, double norm,
             int nchan,int nens,const int nobs[],const EB2D* ens, PRIOR prior);

double csint(double s0,int nchan,int nens,const int nobs[],const EB* ens,
	     int* ngl,double xgl[],double lwgl[],
	     PRIOR prior,double* uncertainty);

void csint2(double s1,double s2,
	    int nchan,int nens,const int nobs[],const EB* ens,
	    int* ngl,double xgl[],double lwgl[],PRIOR prior,
	    double* int1,double* int2,
	    double* v11,double* v12,double* v22);

void csint2cut(double s1,double s2,double shi,
	       int nchan,int nens,const int nobs[],const EB* ens,
	       int* ngl,double xgl[],double lwgl[],PRIOR prior,
	       double* int1,double* int2,
	       double* v11,double* v12,double* v22);

void setdlcsn(int nchan, int nens, int nobs[], const EB* ens);

void setdlcsn2d(int nchan, int nens, int nobs[], const EB2D* ens);

double cslimit(double beta,int nchan,int nens,const int nobs[],const EB* ens,
               int* ngl,double xgl[],double lwgl[],
	       PRIOR prior,double* uncertainty);

double cscutlimit(double beta,double smax,
		  int nchan,int nens,const int nobs[],const EB* ens,
		  int* ngl,double xgl[],double lwgl[],
		  PRIOR prior,double* uncertainty);

double galim(double beta,int nchan,int nens,const int nobs[],const EB* ens);

void gausslaguerre(double x[],double lw[],int n,double alpha);

#define GENLIMIT 1

#endif
