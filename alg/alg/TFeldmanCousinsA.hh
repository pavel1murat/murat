//
#ifndef __murat_alg_TFeldmanCousinsA__
#define __murat_alg_TFeldmanCousinsA__

#include "TRandom3.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TNamed.h"

class TFeldmanCousinsA : public TNamed {
public:
  enum {
	MaxNx = 100,            // max Poisson bin
  } ;   
  
  double   fCL;
  double   fMeanBgr;
  double   fMeanSig;
  
  TRandom3 fRn;

  TH1D*    fBgProbHist;
  TH1D*    fBsProbHist;

  double   fBestSig [MaxNx];
  double   fBestProb[MaxNx];

  double   fLhRatio [MaxNx];
  int      fRank    [MaxNx];

  double   fBgProb  [MaxNx];
  double   fBsProb  [MaxNx];

  double   fFactorial[MaxNx]; // precalculate to speed things up
  
  // X(Bg) - Poisson random number sampled from Bg     distribution
  // X(Bs) - Poisson random number sampled from Bg+Sig distribution
  
  long int fNPE;	      // N(pseudo experiments) to throw

  TH1D*    fProbHist;
  TH1D*    fBeltHist;
  
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
  
  ~TFeldmanCousinsA() {};

  void   Init           (double Bgr, double Sig);
  void   InitPoissonDist(double Mean, double* Prob, int NMax);

  int    ConstructInterval(double Bgr, double Sig);

  void   PrintData(const char* Title, char DataType, void* Data, int MaxInd);
  void   PrintProbs(int N);

  // for a given signal and background, calculates probability of a discovery,
  // where discovery is defined as 
  void   DiscoveryProb(double Bgr, double Sig);

  ClassDef(TFeldmanCousinsA,0)
};

#endif
