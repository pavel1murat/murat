///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "murat/gui/TEvdCrvLayer.hh"

ClassImp(murat::TEvdCrvLayer)

namespace murat {
//-----------------------------------------------------------------------------
// counters are TGeoVolumeAssemblies and are registered with and deleted by TGeoManager
//-----------------------------------------------------------------------------
  TEvdCrvLayer::TEvdCrvLayer(const char* Name, int CopyNumber): TEvdSubdetector(Name,CopyNumber) {
  fNHits = 0;
  fListOfCounters = new TObjArray();
  fListOfCounters->SetOwner(kFALSE);
}

//-----------------------------------------------------------------------------
TEvdCrvLayer::~TEvdCrvLayer() {
  delete fListOfCounters;
}

//-----------------------------------------------------------------------------
void TEvdCrvLayer::Clear(Option_t* Opt) {
  int n = GetNCounters();
  for (int i=0; i<n; i++) {
    TEvdCrvCounter* counter = GetCounter(i);
    if (counter->GetNRecoPulses() > 0) {
      counter->Clear(Opt);
    }
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
