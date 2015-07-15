//

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TRandom.h"
#include "TRandom2.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TObjArray.h"

#include "limits/TLHChannel.hh"
#include "limits/TLHCalculator.hh"

//-----------------------------------------------------------------------------
TLHCalculator::TLHCalculator(): TNamed() {
  fPValue         = -1;
  fLhData         = 0;
  fListOfChannels = new TObjArray();
  fLhHist = new TH1F("lh_tot","Combined Likelihood"      ,1000,-200,0);
}

//-----------------------------------------------------------------------------
TLHCalculator::~TLHCalculator() {
  delete fListOfChannels;
}


//-----------------------------------------------------------------------------
int TLHCalculator::Init() {
  TLHChannel  *ch;

  int nbgr = fListOfChannels->GetEntriesFast();

  fLhData         = 0;
  for (int i=0; i<nbgr; i++) {
    ch = (TLHChannel*) fListOfChannels->At(i);
    ch->Init();
    fLhData += ch->GetLhData();
  }
  return 0;
}

//-----------------------------------------------------------------------------
// generate one pseudoexperiment , return 0 or 1 if LHR for this 
// pseudoexperiment is greater than LHR(data)
//-----------------------------------------------------------------------------
int TLHCalculator::GeneratePseudoExperiment(double& Likelihood) {
  //  TObject  *o; 
  int       res, success;
  TLHChannel* ch(NULL);
  //  TFile*    f;
  double    lh;
  
  res = 1; // default: passed

  int nch = fListOfChannels->GetEntriesFast();

  Likelihood = 0;

  for (int i=0; i<nch; i++) {
    ch          = (TLHChannel*) fListOfChannels->At(i);
    success     = ch->PseudoExperiment(lh);
    if (success < 0) printf("pseudoexperiment i=%10i failed\n",i);
    Likelihood += lh ;
  }

  if (Likelihood > fLhData) {
    res = 0;
  }
  return res;
}

//-----------------------------------------------------------------------------
// generate one pseudoexperiment , return 0 or 1 if generated P-value is 
// greater than P-value of the data
//-----------------------------------------------------------------------------
int TLHCalculator::GeneratePseudoExperiments(double NExp) {

  double lh;
  int    success;

  fNPassed = 0;
  fNMax    = NExp;
  for (double ipe=0; ipe<fNMax; ipe++) {

					// generate one pseudoexperiment
					// success = 0 or 1;

    success = GeneratePseudoExperiment(lh);

    fLhHist->Fill(lh );
    
    if (success == 1) {
      fNPassed += 1;
      //      printf(" success : ipe = %12.0f\n",ipe);
    }
  }
//-----------------------------------------------------------------------------
// P-value: fNPassed/fNMax
//-----------------------------------------------------------------------------
  fPValue = fNPassed/fNMax;

  printf("total = %12.5e, passed = %12.5e, P-value = %12.5e \n",
	 fNMax, fNPassed, fPValue);

  return 0;
}


//-----------------------------------------------------------------------------
void TLHCalculator::Print(const char* Opt) const {
  printf(" TLHCalculator::Print not implemented yet\n");
}
