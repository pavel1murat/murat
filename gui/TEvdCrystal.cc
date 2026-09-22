///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <format>

#include "murat/gui/TEvdCrystal.hh"

ClassImp(murat::TEvdCrystal)

namespace murat {
//-----------------------------------------------------------------------------
TEvdCrystal::TEvdCrystal(const char* Name, TGeoShape* Shape, TGeoMedium* Medium):
  TGeoVolume(Name,Shape,Medium) {
}

//-----------------------------------------------------------------------------
void TEvdCrystal::Print(Option_t* Opt) const {
  std::cout << std::format("emoe\n");
}

//-----------------------------------------------------------------------------
void TEvdCrystal::PrintA() const {
  std::cout << std::format("emoe AAAAAAA\n");
}

}
