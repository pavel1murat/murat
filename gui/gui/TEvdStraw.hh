//-----------------------------------------------------------------------------
#ifndef __murat_gui_TEvdStraw__
#define __murat_gui_TEvdStraw__

#include "Offline/DataProducts/inc/StrawId.hh"
#include "murat/gui/TEvdNumerology.hh"
#include "murat/gui/TEvdSubdetector.hh"

namespace murat {
class TEvdPanel;

class TEvdStraw : public TEvdSubdetector {
public:
  int        fNumber;
  int        fID;       // straw index - channel ID
  int        fPlane;
  int        fPanel;
  int        fLayer;
  double     fZ;
  double     fRho;
  double     fNx;
  double     fNy;
  double     fHalfLength;

  TEvdStraw(int I = -1);
  ~TEvdStraw();
  
  void Init(int Index, int Plane, int Panel, int Layer, double rho, double z, double nx, double ny, double HalfLength);

  double   Z    () { return fZ    ; }  // *MENU*
  double   Rho  () { return fRho  ; }  // *MENU*
  int      ID   () { return fID   ; }  // *MENU*
  void     Print(Option_t* Opt) const; // *MENU*

  ClassDef(murat::TEvdStraw,0);
};
}
#endif
