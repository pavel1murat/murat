//-----------------------------------------------------------------------------
#ifndef __murat_gui_TEvdCrystal__
#define __murat_gui_TEvdCrystal__

#include "TGeoVolume.h"
#include "Stntuple/obj/TCaloHit.hh"

namespace murat {

//-----------------------------------------------------------------------------
class TEvdCrystal: public TGeoVolume {
public:

  TObjArray* fListOfHits;           // pointers to TCaloHits ==NOT OWNED==

  float      fEDep;                 // total deposited energy

  TEvdCrystal(const char* Name, TGeoShape* Shape, TGeoMedium* Medium);

  TObjArray* ListOfHits() { return fListOfHits; }

  void         AddHit(TCaloHit* Hit);
//-----------------------------------------------------------------------------
// overloaded functions of TObject
//-----------------------------------------------------------------------------
  virtual void Clear (Option_t* Opt = "") override;
  virtual void Print (Option_t* Opt = "") const override;  // *MENU*

  virtual void PrintA()                   const ;          // *MENU*

  ClassDefOverride(murat::TEvdCrystal,0)
};

}
#endif
