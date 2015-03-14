#ifndef zzx_limits_TLHChannel
#define zzx_limits_TLHChannel

#include "TNamed.h"
#include "TString.h"
#include "TH1.h"
#include "TObjArray.h"

#include "TLHProcess.hh"

//-----------------------------------------------------------------------------
class TLHChannel: public TNamed {
public:

  TH1*        fProbHist;	        // ! probability dist for a single event

  TH1*        fLhHist;                  // ! likelihood distribution  for PE

  double      fLhData;			// ! data LHR for this channel

  TObjArray*  fListOfProcesses;

  double      fNExpected;

public:

  TLHChannel();
  TLHChannel(const char* Name, const char* Title = "");
  ~TLHChannel();

  double       GetLhData         () { return fLhData   ; }
  double       GetNExpected      () { return fNExpected; }
  TObjArray*   GetListOfProcesses() { return fListOfProcesses; }

  void  AddProcess(TLHProcess* Process) { 
    fListOfProcesses->Add(Process); 
  }
					// to be overloaded

  virtual int Init();
  virtual int PseudoExperiment(double& Likelihood) = 0;

  ClassDef(TLHChannel,0) 
};

#endif
