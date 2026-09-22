//-----------------------------------------------------------------------------
#ifndef __murat_gui_TEvdSubdetector__
#define __murat_gui_TEvdSubdetector__

#include "TGeoVolume.h"
#include "TGeoShape.h"
#include "TGeoMedium.h"

namespace murat {
//-----------------------------------------------------------------------------
class TEvdSubdetector : public TGeoVolume {
public:
  int         fCopyNumber;                    // for a bar 
  TGeoVolume* fTopVolume = {nullptr}; // if null, add this to the geo tree

  TObjArray* fListOfHits;  // if null pointer, then no hits

  TEvdSubdetector();
  TEvdSubdetector(const char* Name,  TGeoShape* Shape, TGeoMedium* Medium);
  ~TEvdSubdetector();

                                        // to hide the inheritance
  TGeoVolume* GetVolume() { return this; }

  int         CopyNumber() { return fCopyNumber; }

  virtual int InitGeometry(const char* Filename); // needed
  
  virtual int InitEvent(); // maybe .. there could be multiple ways of initializing

  TGeoVolume* TopVolume() { return fTopVolume; }
  
  ClassDef(TEvdSubdetector,0);
};
}
#endif
