#ifndef murat_plot_muon_decay_hh
#define murat_plot_muon_decay_hh

#include "TNamed.h"

class TMuonDecay : public TNamed {
public:
  
  static double   fEmu;
  static double   fMmu;
  static double   fMe;
  static double   fMAl;
  static double   fZAl;
  static double   fAlpha;
  static double   fAPi;

//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
public:
  TMuonDecay();

  static  double  ff        (double* X, double* P);
  static  double  Li2       (double X);
  static  double  R1Arbuzov (double Z);

  static  double  f_li2    (double* X, double* P);
  static  double  f_val_li2(double* X, double* P);

  static  double  muon_spectrum_lo_ross  (double  E);
  static  double  f_muon_spectrum_lo_ross(double* X, double* P);

  static  double  muon_spectrum_lo_arbuzov  (double  E);
  static  double  f_muon_spectrum_lo_arbuzov(double* X, double* P);

  static  double  muon_spectrum_nlo_arbuzov  (double E);
  static  double  f_muon_spectrum_nlo_arbuzov(double* X, double* P);

  static  double  f_muon_spectrum_nlo_kuraev(double* X, double* P);

  
  
  void    plot_li2();
  void    plot_val_li2();

  void    plot_muon_spectrum_lo(const char* Opt = "");

  void    plot_muon_spectrum_nlo(const char* Name, const char* Opt = "");


  ClassDef(TMuonDecay,0)
};


#endif
