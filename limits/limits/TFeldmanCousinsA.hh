//
#ifndef __murat_limits_TFeldmanCousinsA__
#define __murat_limits_TFeldmanCousinsA__

#include "TRandom3.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TNamed.h"
#include  "TGraph.h"

#include "murat/limits/analysis.hh"

namespace murat {
  
class TFeldmanCousinsA : public TNamed {
public:
  enum {
	MaxNx =    50,                  // max Poisson bin
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

  struct Belt_t {
    int    fNy;                   // use not to reinitialize
    double fBgr;
    double fSMin;
    double fSMax;
    double fDy;
    int    fCont[MaxNy][MaxNx];  // either 0 or 1 assume 10000 to be large enough for MaxNy
    double fSign[MaxNx][2];      // fSBelt[ix][0] = SMin[ix], fSBelt[ix][1] = SMax[ix]
    int    fIndx[MaxNx][2];      // iymin, iymax
  } fBelt;
  
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
  TFeldmanCousinsA(const char* Name,
		   double      CL,
		   int         DebugLevel = 0);
  
  ~TFeldmanCousinsA();

  void   SetCL          (double CL);

  void   SetNExp        (long int N) { fNExp = N; };
  
  void   SetDebugLevel  (int Level ) { fDebugLevel = Level; };
  
  void   Init           (double Bgr, double Sig, int FillHist = 0);
  
  void   InitPoissonDist(double Mean, double* Prob, double* CumProb, int NMax);

  int    ConstructInterval(double Bgr, double Sig);
  
  int    ConstructBelt    (double Bgr, double SMin, double SMax, int NPoints);

  double Factorial(int Ix) { return fFactorial[Ix]; }

  void   PrintData(const char* Title, char DataType, void* Data, int MaxInd);
  void   PrintProbs(int N);

  // plot discovery probability for a given background and signal range
  // calling makes sense only if CL=-1
  // discovery corresponds to prob=50%

  void   DiscoveryProb(double MuB               , double SMin, double SMax, int NPoints, double* MuS, double* Prob);
  void   DiscoveryProb(murat::analysis* Analysis, double SMin, double SMax, int NPoints, double* MuS, double* Prob);
  
  // plot discovery probability for a given background and signal range
  // calling makes sense only if CL=-1
  // discovery corresponds to NSig = 5
  void   DiscoveryProbMean(double MuB, double SMin, double SMax, int NPoints, double* MuS, double* NSig);

  void   PlotDiscoveryProbMean(double MuB, double Mu2);

  int    SolveFor(double Val, const double* X, const double* Y, int NPoints, double* XVal);

  void   UpperLimit(double MuB, double SMin, double SMax, int NPoints, double* S, double* Prob);
  
  ClassDef(TFeldmanCousinsA,0)
};
}
#endif
