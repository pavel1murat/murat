//-----------------------------------------------------------------------------
#ifndef __murat_gui_TEvdTracker__
#define __murat_gui_TEvdTracker__

#include "TEveElement.h"

#include "murat/gui/TEvdNumerology.hh"
#include "murat/gui/TEvdPlane.hh"
#include "murat/gui/TEvdStation.hh"

//-----------------------------------------------------------------------------  
class TEvdTracker: public TEveElementList {
public:
  TEvdStation*    fStation[kNStations];

  TEvdTracker();

  TEvdPanel* Panel(int Station, int Plane, int Panel) {
    return fStation[Station]->fPlane[Plane]->fPanel[Panel];
  }

  int InitGeometry(const char* Fn);

  ClassDef(TEvdTracker,0)
};

#endif
