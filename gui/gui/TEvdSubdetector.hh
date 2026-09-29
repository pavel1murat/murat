//-----------------------------------------------------------------------------
#ifndef __murat_gui_TEvdSubdetector__
#define __murat_gui_TEvdSubdetector__

#include "TGeoVolume.h"
#include "TGeoShape.h"
#include "TGeoMedium.h"

namespace murat {
//-----------------------------------------------------------------------------
class TEvdSubdetector : public TGeoVolumeAssembly {
public:
  TString     fName;
  int         fCopyNumber;                 // for a bar 
  TGeoVolume* fTopVolume = {nullptr};      // if null, add 'this' to the geo tree

  TObjArray*  fListOfHits;                 // if null pointer, then no hits

  TObjArray*  fListOfSubdetectors;   

  TEvdSubdetector();
  TEvdSubdetector(const char* Name);
  
  TEvdSubdetector(const char* Name,  TGeoShape* Shape, TGeoMedium* Medium,
                  int Color = kGreen+2, int Transparency = 0);
  
  ~TEvdSubdetector();
  
  void AddSubdetector(TEvdSubdetector* sd);
  
                                        // to hide the inheritance
  
  const char* GetName() const { return fName.Data(); }
  TGeoVolume* GetVolume() { return this; }

  int         CopyNumber() { return fCopyNumber; }

  virtual int InitGeometry(const char* Filename); // needed
  
  virtual int InitEvent(); // maybe .. there could be multiple ways of initializing

  TGeoVolume* TopVolume() { return fTopVolume; }
  
  ClassDefOverride(TEvdSubdetector,0);
};
}
#endif
