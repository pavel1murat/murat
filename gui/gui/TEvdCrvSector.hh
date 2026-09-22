//-----------------------------------------------------------------------------
#ifndef __murat_gui_TEvdCrvSector__
#define __murat_gui_TEvdCrvSector__

#include "murat/gui/TEvdSubdetector.hh"

namespace murat {

//-----------------------------------------------------------------------------  
class TEvdCrvSector: public TEvdSubdetector {
public:
  // each sector has its mother volume, the CRV as a whole - doesn't
  // sectors are placed directly into a global TOP
  // a sector is a subdetector
  
  TObjArray* fListOfModules;
  
  TEvdCrvSector();
  ~TEvdCrvSector();

  virtual int InitGeometry(const char* Fn) override;
  //  int InitGeometry(const char* Fn);

  virtual void Print(Option_t* Opt="") const override;

  ClassDef(murat::TEvdCrvSector,0)
};
}
#endif
