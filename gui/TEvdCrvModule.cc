///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "murat/gui/TEvdCrvModule.hh"

ClassImp(murat::TEvdCrvModule)

namespace murat {
//-----------------------------------------------------------------------------
  TEvdCrvModule::TEvdCrvModule(const char* Name, int CopyNumber): TEvdSubdetector(Name,CopyNumber) {
    fNHits        = 0;
    fListOfLayers = new TObjArray();
  }

//-----------------------------------------------------------------------------
  TEvdCrvModule::~TEvdCrvModule() {
    fListOfLayers->Delete();
    delete fListOfLayers;
  }

//-----------------------------------------------------------------------------
void TEvdCrvModule::Clear(Option_t* Opt) {
  int n = GetNLayers();
  for (int i=0; i<n; i++) {
    TEvdCrvLayer* layer = GetLayer(i);
    layer->Clear(Opt);
  }
  fNHits = 0;
}

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
int TEvdCrvModule::InitGeometry(const char* Fn) {
  return 0;
}

//-----------------------------------------------------------------------------
void TEvdCrvModule::Print(Option_t* Opt) const {
  printf("TEvdCrvModule: emoe\n");
}

}
