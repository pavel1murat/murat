///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "TH1F.h"
#include "TObjString.h"
#include "TCanvas.h"

#include "math.h"
#include "string.h"
#include "limits/TLHChannel_ZZ_2lmet.hh"

ClassImp(TLHChannel_ZZ_2lmet)

//-----------------------------------------------------------------------------
TLHChannel_ZZ_2lmet::TLHChannel_ZZ_2lmet() {
  fLhData = -1;
  //  fH0     = 0;   // probability density for a single event
}

//-----------------------------------------------------------------------------
TLHChannel_ZZ_2lmet::~TLHChannel_ZZ_2lmet() {
}

//-----------------------------------------------------------------------------
// return code = 1 : pseudoexperiment succeeded
//             = 0 : pseudoexperiment failed
//-----------------------------------------------------------------------------
int TLHChannel_ZZ_2lmet::PseudoExperiment(double& Likelihood) {
  int    success;
//-----------------------------------------------------------------------------
// search mode: pseudoexperiment is considered a success if the calculated 
// likelihood is less then the likelihood of the data
//-----------------------------------------------------------------------------
  Likelihood = 0;

  if (Likelihood > fLhData) success = 0;
  else                      success = 1;

  return success;
}



