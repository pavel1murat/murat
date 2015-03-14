#include <stddef.h>
#include "TFile.h"
#include "TROOT.h"
#include "mclimit_csm.h"
#include <iostream>
#include "TH1.h"
#include "TRandom3.h"

using namespace std;

int main(int argc, char **argv)
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

	testhyp_pe->print();
	mymclimit->set_npe(1000000);
	//mymclimit->set_npe(50);
        mymclimit->run_pseudoexperiments();

	cout << "Getting results" << endl;
        Double_t tsobs = mymclimit->ts();
	cout << "ts: " << tsobs << endl;
	cout << "tsbm2: " << mymclimit->tsbm2() << endl;
	cout << "tsbm1: " << mymclimit->tsbm1() << endl;
	cout << "tsbmed: " << mymclimit->tsbmed() << endl;
	cout << "tsbp1: " << mymclimit->tsbp1() << endl;
	cout << "tsbp2: " << mymclimit->tsbp2() << endl;
	cout << "tssm2: " << mymclimit->tssm2() << endl;
	cout << "tssm1: " << mymclimit->tssm1() << endl;
	cout << "tssmed: " << mymclimit->tssmed() << endl;
	cout << "tssp1: " << mymclimit->tssp1() << endl;
	cout << "tssp2: " << mymclimit->tssp2() << endl;

	cout << "CLs: " << mymclimit->cls() << endl;
	cout << "CLb: " << mymclimit->clb() << endl;
	cout << "CLsb: " << mymclimit->clsb() << endl;

	cout << "CLs -2sigma (bkg): " << mymclimit->clsexpbm2() << endl;
	cout << "CLs -1sigma (bkg): " << mymclimit->clsexpbm1() << endl;
	cout << "CLs median  (bkg): " << mymclimit->clsexpbmed() << endl;
	cout << "CLs +1sigma (bkg): " << mymclimit->clsexpbp1() << endl;
	cout << "CLs +2sigma (bkg): " << mymclimit->clsexpbp2() << endl;

	cout << "1-CLb -2sigma (sig): " << mymclimit->omclbexpsm2() << endl;
	cout << "1-CLb -1sigma (sig): " << mymclimit->omclbexpsm1() << endl;
	cout << "1-CLb median  (sig): " << mymclimit->omclbexpsmed() << endl;
	cout << "1-CLb +1sigma (sig): " << mymclimit->omclbexpsp1() << endl;
	cout << "1-CLb +2sigma (sig): " << mymclimit->omclbexpsp2() << endl;

	// these don't work so well yet

	//cout << "1-CLbw -2sigma (sig): " << mymclimit->omclbexpsm2w() << endl;
	//cout << "1-CLbw -1sigma (sig): " << mymclimit->omclbexpsm1w() << endl;
	//cout << "1-CLbw median  (sig): " << mymclimit->omclbexpsmedw() << endl;
	//cout << "1-CLbw +1sigma (sig): " << mymclimit->omclbexpsp1w() << endl;
	//cout << "1-CLbw +2sigma (sig): " << mymclimit->omclbexpsp2w() << endl;


	delete mymclimit;
}
