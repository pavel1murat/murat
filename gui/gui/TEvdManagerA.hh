#ifndef __draw_calo_geometry__
#define __draw_calo_geometry__

#include "TObject.h"
#include "TGeoVolume.h"
#include "TGeoManager.h"

#include "murat/gui/TEvdSubdetector.hh"

namespace murat {

//-----------------------------------------------------------------------------  
class TEvdManagerA : public TNamed {
public:
  TGeoVolume*                            fTop;

  TGeoManager* fGeoManager;
  TObjArray    fListOfViews;                // multiple views
  TObjArray    fListOfSubdetectors;         // each     view

  int          fDisplayCalorimeter;
  int          fDisplayCrv;
  int          fDisplayTracker;

private:
  TEvdManagerA(const char* Fcl);        // configuration file name
  
public:
  static TEvdManagerA* Instance(const char* Fcl = "");

  TGeoManager* GetGeoManager() { return fGeoManager; }
 
  void AddSubdetector(TEvdSubdetector* Sd) ;

  virtual void Draw(Option_t* Opt = "") override;
  
  ClassDefOverride(murat::TEvdManagerA,0)
};
  
}
#endif

