///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef __murat_gui_TEvdTrack_hh__
#define __murat_gui_TEvdTrack_hh__

#include "TPolyLine3D.h"

//-----------------------------------------------------------------------------
class TEvdTrack : public TPolyLine3D {
public:
  TString fName;

  TEvdTrack(const char* Name, double* X0, double* V0, double W, double* ZRange);
  
  virtual const char* GetName() const override;

  virtual void Print(Option_t* Opt = "") const cverride;

  ClassDefOverride(TEvdTrack,0)
};

#endif
