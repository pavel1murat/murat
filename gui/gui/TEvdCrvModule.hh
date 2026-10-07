//-----------------------------------------------------------------------------
#ifndef __murat_gui_TEvdCrvModule__
#define __murat_gui_TEvdCrvModule__

#include "murat/gui/TEvdCrvLayer.hh"
#include "murat/gui/TEvdSubdetector.hh"

namespace murat {

//-----------------------------------------------------------------------------  
class TEvdCrvModule: public TEvdSubdetector {
public:
  // each sector has its mother volume, the CRV as a whole - doesn't
  // sectors are placed directly into a global TOP
  // a sector is a subdetector
  
  int           fNHits;
  int           fFirstCounter;
  int           fNCounters;
  TObjArray*    fListOfLayers;
  
  TEvdCrvModule (const char* Name, int CopyNumber = 1);
  ~TEvdCrvModule();

  int           FirstCounter() { return fFirstCounter; }
  int           GetNCounters() { return fNCounters; }
  int           GetNHits    () { return fNHits; }
  int           GetNLayers() { return fListOfLayers->GetEntriesFast(); }

  void          AddLayer(TEvdCrvLayer* Layer) { fListOfLayers->Add(Layer); }
  
  TEvdCrvLayer* GetLayer(int I) { return (TEvdCrvLayer*) fListOfLayers->At(I); }

  void          IncrementNHits() { fNHits++; }
  
  virtual int   InitGeometry(const char* Fn) override;
  //  int InitGeometry(const char* Fn);

  virtual void  Clear(Option_t* Opt="") override;
  virtual void  Print(Option_t* Opt="") const override;

  ClassDefOverride(murat::TEvdCrvModule,0)
};
}
#endif
