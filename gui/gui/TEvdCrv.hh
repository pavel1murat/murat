//-----------------------------------------------------------------------------
#ifndef __murat_gui_TEvdCrv__
#define __murat_gui_TEvdCrv__

#include "TObjArray.h"
#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"

#include "murat/gui/TEvdCrvSector.hh"

#include "murat/gui/TEvdSubdetector.hh"

namespace murat {

//-----------------------------------------------------------------------------  
class TEvdCrv: public TEvdSubdetector {
public:
  std::unique_ptr<mu2e::CosmicRayShield> fCrvPtr;
  // each sector ahs its mother volume, the CRV as a whole - doesn't
  // sectors are placed directly into a global TOP
  // a sector is a subdetector

  int        fNSectors;
 
  TEvdCrv();
  ~TEvdCrv();

  int GetNSectors() { return fNSectors; }

  int GetCounterLocation(int               Sbid   ,
                          TEvdCrvSector*&  Sector ,
                          TEvdCrvModule*&  Module ,
                          TEvdCrvLayer*&   Layer  ,
                          TEvdCrvCounter*& Counter);

  // CRV 'subdetectors' are sectors
  TEvdCrvSector* GetSector(int I) { return (TEvdCrvSector*) fListOfSubdetectors->At(I); }

  virtual int    InitGeometry(const char* Fn) override;
  virtual int    InitEvent() override;

  // virtual void   Draw (Option_t* Opt = "") override ;
  virtual void   Print(Option_t* Opt = "") const override ;

  ClassDefOverride(murat::TEvdCrv,0)
};
}
#endif
