//-----------------------------------------------------------------------------
#ifndef __murat_gui_TEvdCrvLayer__
#define __murat_gui_TEvdCrvLayer__

#include "murat/gui/TEvdCrvCounter.hh"
#include "murat/gui/TEvdSubdetector.hh"

namespace murat {

//-----------------------------------------------------------------------------  
class TEvdCrvLayer: public TEvdSubdetector {
public:
  // each sector has its mother volume, the CRV as a whole - doesn't
  // sectors are placed directly into a global TOP
  // formally, a layer is a subdetector
  // layer's fListOfSubdetectors will be empty
  
  int        fNHits;
  int        fFirstCounter;
  TObjArray* fListOfCounters;
  
  TEvdCrvLayer (const char* Name, int CopyNumber = 1);
  ~TEvdCrvLayer();

  int             FirstCounter() { return fFirstCounter; }
  int             GetNCounters() { return fListOfCounters->GetEntriesFast(); }
  int             GetNHits    () { return fNHits; }

  TEvdCrvCounter* GetCounter(int I) { return (TEvdCrvCounter*) fListOfCounters->At(I); }

  void AddCounter(TEvdCrvCounter* Counter) { fListOfCounters->Add(Counter); }

  void            IncrementNHits() { fNHits++; }

  virtual int     InitGeometry(const char* Fn) override;
  //  int InitGeometry(const char* Fn);
  virtual void    Clear(Option_t* Opt = "") override;
  virtual void    Print(Option_t* Opt = "") const override;

  ClassDefOverride(murat::TEvdCrvLayer,0)
};

}
#endif
