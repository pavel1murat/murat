///////////////////////////////////////////////////////////////////////////////
// a view 
///////////////////////////////////////////////////////////////////////////////
#ifndef __murat_gui_TEvdCrvView_hh__
#define __murat_gui_TEvdCrvView_hh__

#include "TMarker.h"
#include "TNamed.h"
#include "TObjArray.h"
#include "TGaxis.h"
#include "TVector3.h"

#include "murat/gui/TEvdView.hh"

namespace murat {

class TEvdCrvView: public TEvdView {

public:
  TEvdCrvView(const char* Name, TGLViewer::ECameraType Type, int Index = -1); 

  TEvdCrvView(const char* Name, TGLViewer::ECameraType Type, int Index, const char* Title);

  virtual ~TEvdCrvView();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// overloaded functions of TObject
//-----------------------------------------------------------------------------
//  virtual void Draw(Option_t* Opt = "") override;
  
  // virtual void  Paint               (Option_t* option = "") override;
  // virtual void  Print               (Option_t* option = "") const override;  // *MENU* 

  //  ClassDefOverride(TEvdView,0)
};
}
#endif
