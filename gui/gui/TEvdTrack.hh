///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef __murat_gui_TEvdTrack_hh__
#define __murat_gui_TEvdTrack_hh__

#include "THelix.h"

//-----------------------------------------------------------------------------
class TEvdTrack : public THelix {
public:
  TString fName;

  
  TEvdTrack(const char* Name, double* X0, double* V0, double W, double* ZRange);
  virtual const char* GetName() const;

   ClassDef(TEvdTrack,0)
};

#endif
