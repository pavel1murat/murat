///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef __murat_gui_TEvdTrajectory_hh__
#define __murat_gui_TEvdTrajectory_hh__

#include "TPolyLine3D.h"

//-----------------------------------------------------------------------------
class TEvdTrajectory : public TPolyLine3D {
public:
  TString fName;
  
  TEvdTrajectory(const char* Name, int NPoints, float* X, float* Y, float* Z);
  virtual const char* GetName() const;

   ClassDef(TEvdTrajectory,0)
};

#endif
