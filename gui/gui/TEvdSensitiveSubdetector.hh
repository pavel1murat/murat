//-----------------------------------------------------------------------------
#ifndef __murat_gui_TEvdSensitiveSubdetector__
#define __murat_gui_TEvdSensitiveSubdetector__

#include "TGeoVolume.h"
#include "TGeoShape.h"
#include "TGeoMedium.h"

#include "murat/gui/TEvdSensitiveSubdetector.hh"

namespace murat {
//-----------------------------------------------------------------------------
class TEvdSensitiveSubdetector : public TGeoVolume {
public:
  TString     fName;
  int         fCopyNumber;                 // for a bar

  int         f_EventNumber{-1};
  int         f_SubrunNumber{-1};
  int         f_RunNumber{-1};

  TObjArray*  fListOfHits;                 // if null pointer, then no hits

  TEvdSensitiveSubdetector();
  TEvdSensitiveSubdetector(const char* Name, int CopyNumber = 1);
  
  TEvdSensitiveSubdetector(const char* Name,  TGeoShape* Shape, TGeoMedium* Medium,
                  int Color = kGreen+2, int Transparency = 0);
  
  ~TEvdSensitiveSubdetector();
  
  void        AddSubdetector(TEvdSensitiveSubdetector* sd, TGeoMatrix* Matrix = nullptr);
   
  TGeoVolume* GetVolume () { return this; }

  TObjArray*  GetListOfSubdetectors() { return fListOfSubdetectors; }

  int         CopyNumber() { return fCopyNumber; }

  bool        Initialized();
  virtual int InitGeometry(const char* Filename); // needed
  
  virtual int InitEvent(); // maybe .. there could be multiple ways of initializing

  void        SetCopyNumber(int N) { fCopyNumber = N; }

  //  TGeoVolume* TopVolume() { return fTopVolume; }
//-----------------------------------------------------------------------------
// overloaded functions of TObject
//-----------------------------------------------------------------------------
                                        // to hide the inheritance
  const char* GetName() const override { return fName.Data(); }
  
  
  ClassDefOverride(TEvdSensitiveSubdetector,0);
};
}
#endif
