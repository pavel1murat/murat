//-----------------------------------------------------------------------------
#ifndef __murat_gui_TEvdTracker__
#define __murat_gui_TEvdTracker__

#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "murat/gui/TEvdStation.hh"
#include "murat/gui/TEvdSubdetector.hh"

namespace murat {
//-----------------------------------------------------------------------------  
class TEvdTracker: public TEvdSubdetector {
public:
  std::unique_ptr<mu2e::Tracker> fTrkPtr;
  
  TEvdStation*    fStation[kNStations];

  TEvdTracker();

  TEvdPanel* Panel(int Station, int Plane, int Panel) {
    return fStation[Station]->fPlane[Plane]->fPanel[Panel];
  }

  int InitGeometry(const char* Fn);

  ClassDef(murat::TEvdTracker,0)
};
}
#endif
