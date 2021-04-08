//
#ifndef __murat_alg_TFeldmanCousinsA__
#define __murat_alg_TFeldmanCousinsA__

#include "TRandom3.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TNamed.h"
#include  "TGraph.h"

class TFeldmanCousinsA : public TNamed {
public:
  enum {
	MaxNx = 200,                    // max Poisson bin
	MaxNy = 10000,			// max steps in Mu
  };   

  struct Hist_t {
    TH1D*    fBgProb;
    TH1D*    fBsProb;
    TH1D*    fCumBgProb;
    TH1D*    fCumBsProb;
    TH1D*    fProb;                     // the same as fBsProb, but shifted by half-bin
    TH2D*    fBelt;
  } fHist;
  
  double   fCL;
  double   fLog1mCL;			// log(1-fCL)
  double   fMeanBgr;
  double   fMeanSig;
  
  TRandom3 fRn;

  double   fBestSig  [MaxNx];
  double   fBestProb [MaxNx];

  double   fLhRatio  [MaxNx];
  int      fRank     [MaxNx];

  double   fBgProb   [MaxNx];		//
  double   fBsProb   [MaxNx];

  double   fCumBgProb[MaxNx];           // fCumBgProb[i]: P(X<=i) for BGR-only hyp
  double   fCumBsProb[MaxNx];

  double   fFactorial[MaxNx];	        // precalculate to speed things up
  double   fSBelt    [MaxNx][2];	// belt boundaries
  
  // X(Bg) - Poisson random number sampled from Bg     distribution
  // X(Bs) - Poisson random number sampled from Bg+Sig distribution
  
  long int fNExp;	      // N(pseudo experiments) to throw

  int      fDebugLevel;
  int      fIMin;
  int      fIMax;
  int      fNSummed;          // likely, fIMax-fIMin+1
  double   fProb;

  int      fBelt[MaxNy][MaxNx];     // assume 10000 to be large enough for MaxNy
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  TFeldmanCousinsA(const char* Name,
		   double      CL,
		   int         DebugLevel = 0);
  
  ~TFeldmanCousinsA() {};

  void   SetCL          (double CL);
  
  void   Init           (double Bgr, double Sig);
  void   InitPoissonDist(double Mean, double* Prob, double* CumProb, int NMax);

  int    ConstructInterval(double Bgr, double Sig);
  
  int    ConstructBelt    (double Bgr, double SMin, double SMax, int NPoints, double S[], double Belt[][2]);

  void   PrintData(const char* Title, char DataType, void* Data, int MaxInd);
  void   PrintProbs(int N);

  // plot discovery probability for a given background and signal range
  // calling makes sense only if CL=-1
  // discovery corresponds to prob=50%

  void   DiscoveryProb(double MuBB, double SMin, double SMax, int NPoints, double* MuS, double* Prob);

  // make sure NPoints < 10000
  double UpperLimit(double MuB, double SMin, double SMax, int NPoints);
  
  ClassDef(TFeldmanCousinsA,0)
};

#endif
