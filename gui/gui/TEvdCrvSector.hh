//-----------------------------------------------------------------------------
#ifndef __murat_gui_TEvdCrvSector__
#define __murat_gui_TEvdCrvSector__

#include "murat/gui/TEvdCrvModule.hh"
#include "murat/gui/TEvdSubdetector.hh"

namespace murat {

//-----------------------------------------------------------------------------  
class TEvdCrvSector: public TEvdSubdetector {
public:
  // each sector has its mother volume, the CRV as a whole - doesn't
  // sectors are placed directly into a global TOP
  // a sector is a subdetector
  
  TObjArray* fListOfModules;
  int        fNHits;
  //  int        fNModules;
  int        fNCounters;
  int        fFirstCounter;
  
  TEvdCrvSector (const char* Name, int CopyNumber = 1);
  ~TEvdCrvSector();

  void           AddModule(TObject* Module) { fListOfModules->Add(Module); }

  int            FirstCounter() { return fFirstCounter; }
  int            GetNCounters() { return fNCounters; }
  int            GetNModules () { return fListOfModules->GetEntriesFast(); }
  int            GetNHits    () { return fNHits; }
  TEvdCrvModule* GetModule(int I) { return (TEvdCrvModule*) fListOfModules->At(I); }

  void           IncrementNHits() { fNHits++; }
  
  virtual int    InitGeometry(const char* Fn) override;
  //  int InitGeometry(const char* Fn);

  virtual void   Clear(Option_t* Opt = "") override;
  virtual void   Print(Option_t* Opt = "") const override;

  ClassDefOverride(murat::TEvdCrvSector,0)
};
}
#endif
