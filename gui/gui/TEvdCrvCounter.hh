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
  int        fCopyNumber;
  TObjArray* fListOfHits;           // pointers to TCrvRecoPulses ==NOT OWNED==

  TEvdCrvCounter(const char* Name, TGeoShape* Shape, TGeoMedium* Medium);

  TObjArray* ListOfHits() { return fListOfHits; }

  int CopyNumber() { return fCopyNumber; }

  void AddRecoPulse(TCrvRecoPulse* Hit);

  int GetNRecoPulses() { return fListOfHits->GetEntriesFast(); }

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
