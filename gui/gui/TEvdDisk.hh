//-----------------------------------------------------------------------------
#ifndef __murat_gui_TEvdDisk__
#define __murat_gui_TEvdDisk__

#include "Offline/CalorimeterGeom/inc/Disk.hh"
#include "murat/gui/TEvdCrystal.hh"
#include "murat/gui/TEvdSubdetector.hh"

namespace murat {

//-----------------------------------------------------------------------------
// fMu2eDisk is owned by the Mu2e calorimeter
//-----------------------------------------------------------------------------
class TEvdDisk: public TEvdSubdetector {
public:
  mu2e::Disk*  fMu2eDisk;

  TObjArray*   fListOfCrystals;

  TEvdDisk();
  TEvdDisk(const char* Name, TGeoShape* Shape, TGeoMedium* Medium);

  ~TEvdDisk();

  void         AddCrystal(TEvdCrystal* Crystal);

  TObjArray*   ListOfCrystals() { return fListOfCrystals; }

  TEvdCrystal* Crystal(int I) { return (TEvdCrystal*) fListOfCrystals->At(I); }

  virtual void Print (Option_t* Opt = "") const override;  // *MENU*
  
  ClassDefOverride(murat::TEvdDisk,0)
};

}
#endif
