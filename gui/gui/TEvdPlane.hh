//-----------------------------------------------------------------------------
#ifndef __murat_gui_TEvdPlane__
#define __murat_gui_TEvdPlane__

#include "murat/gui/TEvdPanel.hh"
#include "murat/gui/TEvdSubdetector.hh"

namespace murat {

  class TEvdTracker;

//-----------------------------------------------------------------------------
class TEvdPlane : public TEvdSubdetector {
public:
  int          fNumber;
  TEvdPanel*   fPanel[6];
  TEvdTracker*  fTracker;

  TEvdPlane(TEvdTracker* Tracker, int PlaneNumber = -1);

  TEvdPanel*  Panel(int I)   { return fPanel[I]; }
  int         Number() const { return fNumber;   }

  ClassDef(murat::TEvdPlane,0);
};
}
#endif
