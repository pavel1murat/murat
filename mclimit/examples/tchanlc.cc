#include <stddef.h>
#include "TFile.h"
#include "TROOT.h"
#include "mclimit_csm.h"
#include <iostream>
#include "TH1.h"
#include "TRandom3.h"

using namespace std;

// Use mclimit_csm.C on the lvbb channel -- 955 pb-1 with rates and errors
// from CDF 8286 and CDF 8292
// Calculate chisquare -- but do not include shape errors this time around.

// int main(int argc, char **argv)
int tchanlc()
{

  delete gRandom;
  gRandom = (TRandom*) new TRandom3;

#include "preparetchan.h"

	// get a distribution of chisqure for the different hypotheses

        mclimit_csm* mymclimit = (mclimit_csm*) new mclimit_csm();
	mymclimit->set_null_hypothesis(nullhyp);
	mymclimit->set_test_hypothesis(testhyp);
	mymclimit->set_null_hypothesis_pe(nullhyp_pe);
	mymclimit->set_test_hypothesis_pe(testhyp_pe);
	mymclimit->set_datahist(histo[4][0],"tchan");
	mymclimit->set_npe(5000);
	//	mymclimit->run_pseudoexperiments();
	double sflimit,sfunc,sm2,sm1,smed,sp1,sp2;

	/*  stuff for testing the Bayesian posterior printout */
	mymclimit->bayes_interval_begin = 0;
	mymclimit->bayes_interval_end = 10;
	mymclimit->bayes_interval_step = 0.01;
	mymclimit->bayes_pseudoexperiment_limits = new TH1F("PX","Null Hyp PX Limits",100,0.0,10.0);

	// may need more pseudoexperiments for +-1, 2  sigma to be reliable

	mymclimit->bayes_heinrich_withexpect(0.95,&sflimit,&sfunc,
      	        		   1000,&sm2,&sm1,&smed,&sp1,&sp2);

	for (int i=0;i< (int) mymclimit->bayes_posterior.size();i++)
	  {
	    cout << mymclimit->bayes_posterior_points[i] << " " << mymclimit->bayes_posterior[i] << endl;
	  }
	mymclimit->bayes_pseudoexperiment_limits->Print("all");
		
        cout << " " << sflimit << " " << sm2 << " " << sm1 << " " << smed << " " << sp1 << " " << sp2 << endl;
	delete mymclimit;
}
