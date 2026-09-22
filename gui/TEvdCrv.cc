///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "Offline/ConfigTools/inc/SimpleConfig.hh"

#include "murat/gui/TEvdCrv.hh"

#include "Offline/GeometryService/inc/CosmicRayShieldMaker.hh"

#include "murat/gui/TEvdManagerA.hh"

ClassImp(murat::TEvdCrv)

namespace murat {
//-----------------------------------------------------------------------------
TEvdCrv::TEvdCrv(): TEvdSubdetector() {
}

//-----------------------------------------------------------------------------
TEvdCrv::TEvdCrv(const char* Fn): TEvdSubdetector() {
}

//-----------------------------------------------------------------------------
TEvdCrv::~TEvdCrv() {
}

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
int TEvdCrv::InitGeometry(const char* Fn) {
  mu2e::CosmicRayShieldMaker crv_maker(mu2e::SimpleConfig(Fn),0);
  fCrvPtr    = crv_maker.getCosmicRayShieldPtr();

  fNSectors = fCrvPtr->getCRSScintillatorShields().size();
  
  // and initialize the sectors

  auto vm = TEvdManagerA::Instance(); // a

  for (int i=0; i<fNSectors; i++) {
    TEvdCrvSector* sector = Sector(i);
    vm->AddSubdetector(sector);
  }
  
  return 0;
}

//-----------------------------------------------------------------------------
void TEvdCrv::Print(Option_t* Opt) const {
  printf("TEvdCrv: emoe\n");
}

}
