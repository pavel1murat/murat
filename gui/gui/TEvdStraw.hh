//-----------------------------------------------------------------------------
#ifndef __murat_gui_TEvdStraw__
#define __murat_gui_TEvdStraw__

#include "TEveGeoShape.h"

#include "murat/gui/TEvdNumerology.hh"

class TEvdStraw : public TEveGeoShape {
public:
  int        fNumber;
  int        fIndex;  // straw index - channel ID
  double     fZ;
  double     fRho;
  double     fNx;
  double     fNy;
  double     fHalfLength;

  TEvdStraw(int I = -1);
  ~TEvdStraw();
  
  void Init(int Index, double rho, double z, double nx, double ny, double HalfLength);

  double   Z    () { return fZ    ; }  // *MENU*
  double   Rho  () { return fRho  ; }  // *MENU*
  int      Index() { return fIndex; }  // *MENU*
  void     Print(Option_t* Opt) const; // *MENU*

  ClassDef(TEvdStraw,0);
};

#endif
