///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "TH1F.h"
#include "TObjString.h"
#include "TCanvas.h"

#include "math.h"
#include "string.h"
#include "limits/TLHChannel_ZZ_2l2j.hh"

ClassImp(TLHChannel_ZZ_2l2j)

//-----------------------------------------------------------------------------
TLHChannel_ZZ_2l2j::TLHChannel_ZZ_2l2j() {
  fLhData = -1;
  //  fH0     = 0;   // probability density for a single event
}

//-----------------------------------------------------------------------------
TLHChannel_ZZ_2l2j::~TLHChannel_ZZ_2l2j() {
}

//-----------------------------------------------------------------------------
// return code = 1 : pseudoexperiment succeeded
//             = 0 : pseudoexperiment failed
//-----------------------------------------------------------------------------
int TLHChannel_ZZ_2l2j::PseudoExperiment(double& Likelihood) {
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



