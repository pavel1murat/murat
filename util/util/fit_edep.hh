#ifndef __murat_util_fit_edep_hh__
#define __murat_util_fit_edep_hh__

#include "TH1F.h"
#include "TH1D.h"

//-----------------------------------------------------------------------------
class fit_edep {

public:
  enum {
    kNPanels     =   6,
    kNStraws     =  96,
    kMaxRuns     = 100,                      // currently available: 100003-100033
    kNPar        =   4,
  };

  struct fit_result_t {
    int          rn   {-1};                             // run number
    int          panel{-1};                          // panel
    TH1F*        h_edep;
    double       p[kNPar];
    double       e[kNPar];
    double       chi2d{-1.};
    int          ndf  {-1 };
    std::string  msg;
  };

//-----------------------------------------------------------------------------
// data
//-----------------------------------------------------------------------------
  int          fRunNumber;
  fit_result_t frr   [216];
  TH1D*        h_edep[216];
  TF1*         f_fit;

  fit_edep(int RunNumber = -1) { fRunNumber = RunNumber;};

  void fit_histogram(TH1* Hist, const char* FitOpt, double XMin, double XMax, fit_result_t* Fr);
  int  fit_panels   (int IMin, int IMax, const char* FileName);

  void plot_summary_histograms();

  static double fit_function(double* X, double* P);
};

#endif
