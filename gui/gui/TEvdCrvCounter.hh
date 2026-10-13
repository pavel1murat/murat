//-----------------------------------------------------------------------------
#ifndef __murat_gui_TEvdCrvCounter__
#define __murat_gui_TEvdCrvCounter__

#include "TGeoVolume.h"
#include "Stntuple/obj/TCrvRecoPulse.hh"

namespace murat {
//-----------------------------------------------------------------------------
// inherits from TGeoVolume, not TGeoVolumeAssembly,
// thus - not a subdetector
//-----------------------------------------------------------------------------
class TEvdCrvCounter: public TGeoVolume {
public:
  int         fCopyNumber;
  int         fSectorType;
  int         fSectorID;
  int         fModuleID;
  int         fLayerID;
  int         fSbid;
  TObjArray*  fListOfHits;           // pointers to TCrvRecoPulses ==NOT OWNED==
  TGeoMatrix* fGeoMatrix;           // ==NOT OWNED==

  TEvdCrvCounter(const char* Name, TGeoShape* Shape, TGeoMedium* Medium);

  TObjArray* ListOfHits() { return fListOfHits; }

  int CopyNumber() { return fCopyNumber; }

  int GetSectorType() { return fSectorType;   }
  int GetSectorID  () { return fSectorID; }
  int GetModuleID  () { return fModuleID; }
  int GetLayerID   () { return fLayerID;  }
  int GetSbid      () { return fSbid; }

  void AddRecoPulse(TCrvRecoPulse* Hit);

  int GetNRecoPulses() { return fListOfHits->GetEntriesFast(); }

  TGeoMatrix*  GetGeoMatrix() { return fGeoMatrix; }

  void SetCopyNumber(int N) { fCopyNumber = N; }

  TCrvRecoPulse* GetRecoPulse(int I) { return (TCrvRecoPulse*) fListOfHits->At(I); }
//-----------------------------------------------------------------------------
// overloaded functions of TObject
//-----------------------------------------------------------------------------
  virtual void Clear (Option_t* Opt = "") override;
  virtual void Print (Option_t* Opt = "") const override;  // *MENU*

  virtual void PrintA()                   const ;          // *MENU*

  ClassDefOverride(murat::TEvdCrvCounter,0)
};

}
#endif
