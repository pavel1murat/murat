///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "murat/gui/TEvdCrvSector.hh"

ClassImp(murat::TEvdCrvSector)

namespace murat {
//-----------------------------------------------------------------------------
TEvdCrvSector::TEvdCrvSector(): TEvdSubdetector() {
  fListOfModules = new TObjArray();
}

//-----------------------------------------------------------------------------
TEvdCrvSector::~TEvdCrvSector() {
  fListOfModules->Delete();
  delete fListOfModules;
}

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
int TEvdCrvSector::InitGeometry(const char* Fn) {
  return 0;
}

//-----------------------------------------------------------------------------
void TEvdCrvSector::Print(Option_t* Opt) const {
  printf("TEvdCrvSector: emoe\n");
}

}
