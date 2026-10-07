///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "murat/gui/TEvdCrvSector.hh"

ClassImp(murat::TEvdCrvSector)

namespace murat {
//-----------------------------------------------------------------------------
  TEvdCrvSector::TEvdCrvSector(const char* Name, int CopyNumber):
    TEvdSubdetector(Name,CopyNumber)
  {
    fNHits         = 0;
    fListOfModules = new TObjArray();
  }

//-----------------------------------------------------------------------------
TEvdCrvSector::~TEvdCrvSector() {
}

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
int TEvdCrvSector::InitGeometry(const char* Fn) {
  return 0;
}

//-----------------------------------------------------------------------------
void TEvdCrvSector::Clear(Option_t* Opt) {
  int n = GetNModules();
  for (int i=0; i<n; i++) {
    TEvdCrvModule* module = GetModule(i);
    module->Clear(Opt);
  }
  fNHits = 0;
}
  
//-----------------------------------------------------------------------------
void TEvdCrvSector::Print(Option_t* Opt) const {
  printf("TEvdCrvSector: emoe\n");
}

}
