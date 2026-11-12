//-----------------------------------------------------------------------------
#ifndef __murat_gui_TEvdTracker__
#define __murat_gui_TEvdTracker__

#include "Offline/TrackerGeom/inc/Tracker.hh"

#include "murat/gui/TEvdStation.hh"
#include "murat/gui/TEvdPanel.hh"
#include "murat/gui/TEvdSubdetector.hh"

namespace murat {
//-----------------------------------------------------------------------------  
class TEvdTracker: public TEvdSubdetector {
public:
  std::unique_ptr<mu2e::Tracker> fTrkPtr;
  
  TEvdStation*    fStation[kNStations];
  int             fDisplayPanels{0};       // 1: individual panel volumes
  int             fDisplayStraws{1};       // 1: straw volumes inside the panels
  int             fPanelTransparency{97};  // ROOT scale: 0 opaque, 100 invisible

  TEvdTracker();

  TEvdPanel*   Panel(int Station, int Plane, int Panel);

  virtual int  InitEvent   ()               override;
  virtual int  InitGeometry(const char* Fn) override;
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  virtual void Print(Option_t* Opt = "") const override;

  ClassDefOverride(murat::TEvdTracker,0)
};
}
#endif
