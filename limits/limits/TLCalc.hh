#ifndef zzx_limits_TLHCalc
#define zzx_limits_TLHCalc

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TRandom.h"
#include "TRandom2.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TPad.h"

// #include "murat/ana/TAnaUtils.hh"
#include "Stntuple/val/stntuple_val_functions.hh"

class TLCalc : public TNamed {
public:
  struct Hist_t {
    TH1*   fMass;      // ! sampled mass distribution
    TH1F*  fLhMass;    // ! mass likelihood;
    TH1*   fPt;        // ! sampled PT distribution
    TH1F*  fLhPt;      // ! mass likelihood;
    TH1F*  fLh2d;      // ! 2D likelihood histogram;
    TH2*   fPtVsMass;  // 
    TH1F*  fNEvents[2];  // ! 
    TH1F*  fNHighMass[2];// ! 
    TH2F*  fNHighClose;  // ! 
  } fHist;

  double fData[8][2];
  int    fNMax;
  int    fMode;          // =1: PYTHIA, =2: MC@NLO+HERWIG
  double fLhMassData;    // ! data likelihood position;
  double fLhPtData;      // ! data likelihood position;
  double fLh2dData;      // ! 2D data likelihood ;
  double fMassWindow;    // !
  double fNMean;         // !
  double fMinMass;       // !

 
  TLCalc(int Mode = 1);
  ~TLCalc();

  int calc_lh_2d   (int NMax = 100000000);
  int calc_lh_mass (int NMax =  10000000);
  int calc_4l_prob (double Mass = 300., int NMax =  10000000);
  int calc_lh_pt   (int NMax =  10000000);
  int calc_lh_pt_hm(int NMax =  10000000);

  void  InitHistograms();

};


#endif
