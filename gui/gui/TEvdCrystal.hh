//-----------------------------------------------------------------------------
#ifndef __murat_gui_TEvdCrystal__
#define __murat_gui_TEvdCrystal__

#include "TGeoVolume.h"
#include "Stntuple/obj/TCaloHit.hh"

class TCanvas;

namespace murat {

//-----------------------------------------------------------------------------
class TEvdCrystal: public TGeoVolume {
public:
  int        fCid;
  TObjArray* fListOfHits;           // pointers to TCaloHits ==NOT OWNED==
  float      fEDep;                 // total deposited energy
  TObjArray* fListOfHistograms;

  static TCanvas* fgCanvasWf;           // to plot waveforms

  TEvdCrystal(int Index, TGeoShape* Shape, TGeoMedium* Medium);

  TObjArray* ListOfHits() { return fListOfHits; }

  int        Cid  () { return fCid; }
  void       AddHit(TCaloHit* Hit);
  int        NHits() { return fListOfHits->GetEntriesFast(); }
  TCaloHit*  Hit(int I) { return (TCaloHit*) fListOfHits->UncheckedAt(I); }
//-----------------------------------------------------------------------------
// overloaded functions of TObject
//-----------------------------------------------------------------------------
  virtual void Clear (Option_t* Opt = "") override;

  void PrintHits    ();          // *MENU*
  void PlotWaveforms();          // *MENU*
  

  ClassDefOverride(murat::TEvdCrystal,0)
};

}
#endif
