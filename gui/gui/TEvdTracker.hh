//-----------------------------------------------------------------------------
#ifndef __murat_gui_TEvdTracker__
#define __murat_gui_TEvdTracker__

#include "murat/gui/TEvdNumerology.hh"
#include "murat/gui/TEvdPanel.hh"

//-----------------------------------------------------------------------------
class TEvdPlane {
public:
  int         fNumber;
  TEvdPanel*  fPanel[6];

  TEvdPlane(int I = -1) {
    fNumber = I;
    for (int i=0; i<kNPanels; i++) {
      fPanel[i]         = new TEvdPanel(i);
      fPanel[i]->fPlane = this;
    }
  }
  
  TEvdPanel*  Panel(int I) { return fPanel[I]; }

  int         Number() const { return fNumber; }
};

//-----------------------------------------------------------------------------
class TEvdStation {
public:
  int         fNumber;
  TEvdPlane*  fPlane[2];

  TEvdStation(int I = -1) {
    fNumber = I;
    for (int i=0; i<2; i++) fPlane[i] = new TEvdPlane(i);
  }

  TEvdPlane*  Plane(int I) { return fPlane[I]; }
};

//-----------------------------------------------------------------------------  
class TEvdTracker {
public:
  TEvdStation*    fStation[kNStations];

  TEvdTracker() {
    for (int i=0; i<kNStations; i++) fStation[i] = new TEvdStation(i);
  }

  TEvdPanel* Panel(int Station, int Plane, int Panel) {
    return fStation[Station]->fPlane[Plane]->fPanel[Panel];
  }
};

#endif
