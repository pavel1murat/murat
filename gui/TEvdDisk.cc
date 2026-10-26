///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <format>

#include "murat/gui/TEvdDisk.hh"

ClassImp(murat::TEvdDisk)

namespace murat {

TEvdDisk::TEvdDisk() : TEvdSubdetector() {
  fMu2eDisk       = nullptr;
  fListOfCrystals = nullptr;
}
//-----------------------------------------------------------------------------
TEvdDisk::TEvdDisk(const char* Name, TGeoShape* Shape, TGeoMedium* Medium):
  TEvdSubdetector(Name,Shape,Medium) {
  fMu2eDisk       = nullptr;
  fListOfCrystals = new TObjArray();
  fListOfCrystals->SetOwner(kFALSE);
}

//-----------------------------------------------------------------------------
TEvdDisk::~TEvdDisk() {
  if (fListOfCrystals) {
    delete fListOfCrystals;
  }
}

//-----------------------------------------------------------------------------
void TEvdDisk::AddCrystal(TEvdCrystal* Crystal) {
  // GetVolume()->Add
  fListOfCrystals->Add(Crystal);
}
  
//-----------------------------------------------------------------------------
void TEvdDisk::Print(Option_t* Opt) const {
  std::cout << std::format("TEvdDisk emoe \n");
}

}
