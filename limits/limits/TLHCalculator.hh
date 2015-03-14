///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////

#ifndef zzx_limits_TLHCalculator
#define zzx_limits_TLHCalculator

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TRandom.h"
#include "TRandom2.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TPad.h"

#include "TLHChannel.hh"

//-----------------------------------------------------------------------------
class TLHCalculator : public TNamed {
public:

  double  fNMax;
  double  fNPassed;
					// need to make sure fSig and fBgr 
					// have the same number of channels
  double  fNChannels;   
  double  fPValue;
  double  fLhData;

  TObjArray*  fListOfChannels;

  TH1F*       fLhHist;			// combined likelihood
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
public:
  TLHCalculator();
  ~TLHCalculator();

  int Init();
  int GeneratePseudoExperiment(double& Likelihood);

  int GeneratePseudoExperiments(double NExp);

  void AddChannel(TObject* Channel) { fListOfChannels->Add(Channel); }

  TLHChannel*   GetChannel(const char* Name) {
    return (TLHChannel*) fListOfChannels->FindObject(Name);
  }

  virtual void  Print(const char* Opt) const ;

  ClassDef(TLHCalculator,0)
};

#endif
