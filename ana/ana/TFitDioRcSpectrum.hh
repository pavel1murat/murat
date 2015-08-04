//
#ifndef murat_ana_TFitDioRcSpectrum_hh
#define murat_ana_TFitDioRcSpectrum_hh


#include "TObject.h"
#include "TH1.h"
#include "TF1.h"
#include "TSpline.h"

class TFitDioRcSpectrum: public TObject {
public:
  
  struct Range_t {
    int    fMode;
    double fXMin;
    double fXMax;
  };
    
  TH1F*   Hist;		 
  TF1*    fFunc;	 
  TF1*    fFull;	 
  int     fNPar;

  static  int      fNRanges;
  static  int      fgMode;	 
  static  double   X0;		 
  static  Range_t* fRange;
  static  double   fP[100][10];
  static  int      fgDebugMode;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------  
   TFitDioRcSpectrum();
  ~TFitDioRcSpectrum();

  static double dio_energy_spectrum  (double E);

  static double f_pol_dio            (double* X, double* P);
  static double f_dio_energy_spectrum(double* X, double* P);

  void          fit_pol_dio(int Mode, double XMin, double XMax);

  int           init_hist();

  int           main_pol_dio();

  void          PrintSpectrum(double EMin,  double EMax, double Step);
};

#endif
