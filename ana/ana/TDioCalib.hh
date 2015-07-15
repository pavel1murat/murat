///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////

#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TMath.h"
#include "include/Stntuple/val/stntuple_val_functions.hh"


class TDioCalib: public TObject {

  struct Hist_t {
    TH1F*  fDetResponse;
    TH1F*  fRecoMom;		// reconstructed momentum
    TH1F*  fDioMom;
    TH1F*  fDioExpected;
  };

public:
  TGraphErrors* gr;

  Hist_t  fHist;

  // TH1F*         h_dio_orig;
  // TH1F*         h_dio_expected;


  int    kNBins  = 500 ;
  double kEhmin  = 101.;
  double kEhmax  = 106.;


  TString       fDioFileName;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  TDioCalib();
  ~TDioCalib();

  static double f_crystal_ball(double* X, double* P);

  static double f_dio_spectrum(double* X, double* P);

  void          HalfField(double Scale = 1.0);

  double        ReconstructedMomentum(double P);

  double        TrkRecoEff(double P);

  void          InitDetectorResponseFunction();

  ClassDef(TDioCalib,0)

};


