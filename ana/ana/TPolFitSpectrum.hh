//
#ifndef murat_ana_TPolFitSpectrum_hh
#define murat_ana_TPolFitSpectrum_hh


#include "TObject.h"
#include "TH1.h"
#include "TF1.h"
#include "TSpline.h"

class TPolFitSpectrum: public TNamed {
public:
  
  struct Range_t {
    int    fMode;
    double fXMin;
    double fXMax;
  };
    
  TH1F*   fHist;		 
  TF1*    fFunc;	 
  TF1*    fFull;	 
  int     fNPar;
  double  fStep;

  static  TPolFitSpectrum* fgFit;

  static  int      fNRanges;
  static  int      fgMode;	 
  static  double   X0;		 
  Range_t* fRange;
  static  double   fP[100][10];
  static  int      fgDebugMode;
//-----------------------------------------------------------------------------
// functions;  Name = "dio_rc" or "cnv_rc" for DIO and conversions respectively
//-----------------------------------------------------------------------------
   TPolFitSpectrum(const char* Name = "none");
  ~TPolFitSpectrum();

  static double dio_energy_spectrum  (double E);

  static double f_pol_dio            (double* X, double* P);
  static double f_dio_energy_spectrum(double* X, double* P);

  void          fit_pol_dio(int Mode, double XMin, double XMax);

  int           init_hist();

  int           main_pol_dio();

  void          PrintSpectrum(double EMin,  double EMax, double Step);

  ClassDef(TPolFitSpectrum,0)
};

#endif
