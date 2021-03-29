//
#ifndef __murat_alg_TLogLHR__
#define __murat_alg_TLogLHR__

#include "TRandom3.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TNamed.h"
#include  "TGraph.h"

class TLogLHR : public TNamed {
public:
  enum {
	MaxNx   = 200,			// max Poisson bin
	MaxNy   = 10000,		// max steps in Mu
	NxLogLH = 1000,			// N(bins) in the LogLHr 
  };   

  struct Hist_t {

    TH1D*    fLHPoi[MaxNx];  // i-th hist is a poisson distribuion for <mu> = i
    
    TH1D*    fLHb;
    TH1D*    fLHs;

    TH1D*    fLogLHb;
    TH1D*    fLogLHs;
    TH1D*    fLogLHr;

    TH1D*    fPTail;                     // tail probability
  } fHist;
  
  double   fCL;				// confidence level
  double   fMeanBgr;
  double   fMeanSig;
  
  TRandom3 fRn;

  double   fBestSig  [MaxNx];
  double   fBestProb [MaxNx];

  double   fLhRatio  [MaxNx];

  double   fLHb   [MaxNx];		//
  double   fLHs   [MaxNx];

  // X(Bg) - Poisson random number sampled from Bg     distribution
  // X(Bs) - Poisson random number sampled from Bg+Sig distribution
  
  long int fNExp;	      // N(pseudo experiments) to throw

  int      fDebugLevel;
  int      fIMin;
  int      fIMax;
  int      fNSummed;          // likely, fIMax-fIMin+1
  double   fProb;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  TLogLHR(const char* Name,
		   double      CL,
		   int         DebugLevel = 0);
  
  ~TLogLHR() {};

  void   Init           (double Bgr , double Sig);
  void   InitPoissonDist(double Mean, double* Prob, int NMax);

  void   InitLogLHr(double* Num, double* Denom);

  void   PrintData(const char* Title, char DataType, void* Data, int MaxInd);
  void   PrintProbs(int N);

  double PTail(double LogLHr);

  // plot discovery probability for a given background and signal range
  // calling makes sense only if CL=-1
  // discovery corresponds to prob=50%

  void   DiscoveryProb(double Bgr, double SMin, double SMax, int NSteps= 10);

  // it is up to you to make sure NSteps < 10000
  void ConfInterval(double Bgr, int N, double SMin, double SMax, int NSteps, double* Prob);
  
  ClassDef(TLogLHR,0)
};

#endif
