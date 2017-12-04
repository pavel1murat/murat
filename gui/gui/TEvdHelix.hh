//-----------------------------------------------------------------------------
#ifndef __murat_gui_TEvdHelix__
#define __murat_gui_TEvdHelix__

#include "TEveLine.h"

class TEvdHelix : public TEveLine {
public:
  double fX0;
  double fY0;
  double fZ0;
  double fD0;
  double fOmega;
  double fTanDip; //
  double fPhi0;
  double fZMin;
  double fZMax;
public:
  TEvdHelix();
  TEvdHelix(double Z0, double D0, double Phi0, double Omega, double TanDip, double ZMin, double ZMax);
  ~TEvdHelix();

  void  GetPointAtS(double S, TVector3* Point);
  void  GetPointAtZ(double Z, TVector3* Point);

  double CosDip() { return 1./sqrt(1.+fTanDip*fTanDip); }
  double SinDip() { return fTanDip/sqrt(1.+fTanDip*fTanDip); }
  
  void  Print(Option_t* Opt) const;	// *MENU*

  ClassDef(TEvdHelix,0);
};

#endif
