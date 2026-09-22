///////////////////////////////////////////////////////////////////////////////

#include "murat/gui/TEvdSubdetector.hh"

ClassImp(murat::TEvdSubdetector)

namespace murat {
//-----------------------------------------------------------------------------
TEvdSubdetector::TEvdSubdetector() : TGeoVolume() {
}

TEvdSubdetector::TEvdSubdetector(const char* Name,  TGeoShape* Shape, TGeoMedium* Medium)
  : TGeoVolume(Name,Shape,Medium) {
}

  TEvdSubdetector::~TEvdSubdetector() {
  }

  int TEvdSubdetector::InitGeometry(const char* Fn) {
    return 0;
  }

  int TEvdSubdetector::InitEvent() {
    return 0;
  }
}
