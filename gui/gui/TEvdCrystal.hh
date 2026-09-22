//-----------------------------------------------------------------------------
#ifndef __murat_gui_TEvdCrystal__
#define __murat_gui_TEvdCrystal__

#include "TGeoVolume.h"

namespace murat {

//-----------------------------------------------------------------------------
class TEvdCrystal: public TGeoVolume {
public:

  TEvdCrystal(const char* Name, TGeoShape* Shape, TGeoMedium* Medium);

  virtual void Print (Option_t* Opt = "") const override;  // *MENU*
  virtual void PrintA()                   const ;          // *MENU*

  ClassDefOverride(murat::TEvdCrystal,0)
};

}
#endif
