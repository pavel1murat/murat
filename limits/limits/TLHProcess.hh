#ifndef zzx_limits_TLHProcess
#define zzx_limits_TLHProcess

#include "TNamed.h"
#include "TString.h"
#include "TH1.h"

//-----------------------------------------------------------------------------
class TLHProcess: public TNamed {
public:

  TH1*     fProbHist;			// ! probability distribution for event
  double   fNExpected;			// ! number of expected events
  int      fSignal;                     // ! =1 for the signal, 0 for the bgr

public:

  TLHProcess();
  TLHProcess(const char* Name, const char* Title, double NExpected);
  ~TLHProcess();

  TH1*    GetProbHist () { return fProbHist;  }
  double  GetNExpected() { return fNExpected; }

					// if Module=0, use TFile::Get
					// work with 1D or 2D histograms
  void SetProbHist(const char* Filename, 
		   const char* Module  , 
		   const char* Hist    ,
		   int         NDim    );

  ClassDef(TLHProcess,0) 
};

#endif
