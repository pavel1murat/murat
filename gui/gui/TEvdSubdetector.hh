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

  int         f_EventNumber {-1};
  int         f_SubrunNumber{-1};
  int         f_RunNumber   {-1};

  TObjArray*  fListOfHits;                 // if null pointer, then no hits

  TObjArray*  fListOfSubdetectors;   

  TEvdSubdetector();
  TEvdSubdetector(const char* Name, int CopyNumber = 1);
  
  TEvdSubdetector(const char* Name,  TGeoShape* Shape, TGeoMedium* Medium,
                  int Color = kGreen+2, int Transparency = 0);
  
  ~TEvdSubdetector() override;
  
  void        AddSubdetector(TEvdSubdetector* sd, TGeoMatrix* Matrix = nullptr);
   
  TGeoVolume* GetVolume () { return this; }

  TObjArray*  GetListOfSubdetectors() { return fListOfSubdetectors; }

  int         CopyNumber() { return fCopyNumber; }

  bool        Initialized();
  virtual int InitGeometry(const char* Filename); // needed
  
  virtual int InitEvent(); // maybe .. there could be multiple ways of initializing

  void        SetCopyNumber(int N) { fCopyNumber = N; }
  void        SetInitialized();

  //  TGeoVolume* TopVolume() { return fTopVolume; }
//-----------------------------------------------------------------------------
// overloaded functions of TObject
//-----------------------------------------------------------------------------
                                        // to hide the inheritance
  const char* GetName() const override { return fName.Data(); }
  
  
  ClassDefOverride(TEvdSubdetector,0);
};
}
#endif
