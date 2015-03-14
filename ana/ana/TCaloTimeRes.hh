///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef __murat_ana_TCaloTimeRes__
#define __murat_ana_TCaloTimeRes__

#include "TH1.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TMath.h"


//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
class TCaloTimeRes: public TObject {
public:

  struct Hist_t {
    TH1F* fWaveform;
    TH1F* fRes;
    TH1F* fPull;
  };

  Hist_t fHist;

  double fTau1;
  double fTau2;

  double fNPePerMeV; // light yield: number of "pe" per MeV

  double fClusterEnergy;
  double fNoise;
  double fSigmaNoise;

  double fTimeBin;
  double fMaxFitTime;   // max range of the fitting time

  TRandom3 fRn;

  TF1*    fPulse;

  int    fNEvents;

  int    fNSamples;
  
  TCaloTimeRes();
  ~TCaloTimeRes();

  static double pulse        (double* X, double* P);
  static double pulse2       (double* X, double* P);
  static double func_integral(double t0, double tmin, double tmax);

  void   time_resolution(int NEvents = 1);
};




#endif
