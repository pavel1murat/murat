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
	MaxNx   = 20,			// max Poisson bin
	MaxNy   = 10000,		// max steps in Mu
	NxLogLH = 3000,			// N(bins) in the LogLHr 
  };

  struct MData_t {
    int    N;                           // measured number of N
    double mub;                         // poisson mean of the distributon from which N has been sampled
    double lhb;                         // likelihood : exp(-mu)*mu^N/N!
    double log_lhb;                     // log likelihood;
    double mus;
    double lhs;
    double log_lhs;                     // log likelihood;

    double log_lhr() {
      if (lhb > 1.e-12) return log_lhb - log_lhs ;
      else              return -30;
    }
  };

  struct Hist_t {
    
    TH1D*    fLHb;                      // Poisson dist for the 'best' hyp
    TH1D*    fLHs;			// Poisson dist for the hyp with signal
    TH1D*    fLHt;			// for future

    TH1D*    fLogLHb;
    TH1D*    fLogLHs;
    TH1D*    fLogLHr;
    TH1D*    fLogLHw;
    TH1D*    fLogLHn;			// poisson N for a given log(lhr)

    TH1D*    fPTail;                     // tail probability
  } fHist;
  
  double   fCL;				// confidence level
  double   fMuB;
  double   fMuS;
  
  TRandom3 fRn;

  double   fLHb   [MaxNx];		//
  double   fLHs   [MaxNx];
  double   fLHt   [MaxNx];

  MData_t* fData  [MaxNx];              // use instead of histograms

  double   fFactorial[MaxNx];

  // X(Bg) - Poisson random number sampled from Bg     distribution
  // X(Bs) - Poisson random number sampled from Bg+Sig distribution
  
  long int fNExp;	      // N(pseudo experiments) to throw

  int      fDebugLevel;
  double   fProb;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  TLogLHR(const char* Name,
	  double      CL,
	  int         DebugLevel = 0);
  
  ~TLogLHR() {};

  void   Init           (double Bgr , double Sig, int N, double MuBest, MData_t** Data);
  
  void   InitPoissonDist(double Bgr, double Sig, int N, double* Prob, MData_t** Data, int NMax);

  void   PrintData(MData_t** Data, int MaxInd = -1);

  double PTail(MData_t** Data, double LogLHr);
  
  // it is up to you to make sure NSteps < 10000
  void ConfInterval(double Bgr, int NMeas, double SMin, double SMax, int NSteps, double* Prob);

  void MeasInterval(double Bgr, int NMeas, double SMin, double SMax, int NSteps, double* Prob);
  
  ClassDef(TLogLHR,0)
};

#endif
