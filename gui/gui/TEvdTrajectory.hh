///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef __murat_gui_TEvdTrajectory_hh__
#define __murat_gui_TEvdTrajectory_hh__

#include "TPolyLine3D.h"

class TStrTrack;

//-----------------------------------------------------------------------------
class TEvdTrajectory : public TPolyLine3D {
public:
  TString      fName;
  TStrTrack*   fTrack{nullptr};     // ==NOT OWNED==

  TEvdTrajectory(const char* Name, TStrTrack* Track, int NPoints, float* X, float* Y, float* Z);
  virtual const char* GetName() const;

   ClassDef(TEvdTrajectory,0)
};

#endif
