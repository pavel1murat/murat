///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "murat/gui/TEvdCrvLayer.hh"

ClassImp(murat::TEvdCrvLayer)

namespace murat {
//-----------------------------------------------------------------------------
  TEvdCrvLayer::TEvdCrvLayer(const char* Name, int CopyNumber): TEvdSubdetector(Name,CopyNumber) {
  fNHits = 0;
  fListOfCounters = new TObjArray();
}

//-----------------------------------------------------------------------------
TEvdCrvLayer::~TEvdCrvLayer() {
  fListOfCounters->Delete();
  delete fListOfCounters;
}

//-----------------------------------------------------------------------------
void TEvdCrvLayer::Clear(Option_t* Opt) {
  int n = GetNCounters();
  for (int i=0; i<n; i++) {
    TEvdCrvCounter* counter = GetCounter(i);
    counter->Clear(Opt);
  }
  fNHits = 0;
}
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
int TEvdCrvLayer::InitGeometry(const char* Fn) {
  return 0;
}

//-----------------------------------------------------------------------------
void TEvdCrvLayer::Print(Option_t* Opt) const {
  printf("TEvdCrvLayer: emoe\n");
}

}
