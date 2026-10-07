///////////////////////////////////////////////////////////////////////////////
// a view 
///////////////////////////////////////////////////////////////////////////////
#ifndef __murat_gui_TEvdCaloView_hh__
#define __murat_gui_TEvdCaloView_hh__

#include "TMarker.h"
#include "TNamed.h"
#include "TObjArray.h"
#include "TGaxis.h"
#include "TVector3.h"

#include "murat/gui/TEvdView.hh"

namespace murat {
  
class TEvdCaloView: public TEvdView {

public:
  TEvdCaloView(const char* Name, int Type = -1, int Index = -1); 

  TEvdCaloView(const char* Name, int Type, int Index, const char* Title);

  virtual ~TEvdCaloView();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// overloaded functions of TObject
//-----------------------------------------------------------------------------
  virtual void Draw(Option_t* Opt = "") override;
  
  // virtual void  Paint               (Option_t* option = "") override;
  // virtual void  Print               (Option_t* option = "") const override;  // *MENU* 

  //  ClassDefOverride(TEvdView,0)
};
}
#endif
