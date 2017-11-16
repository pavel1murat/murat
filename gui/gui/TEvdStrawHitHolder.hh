//-----------------------------------------------------------------------------
#ifndef __murat_gui_TEvdStrawHitHolder__
#define __murat_gui_TEvdStrawHitHolder__

#include "TEveGeoShape.h"

#include "murat/gui/TEvdNumerology.hh"

#include "TNamed.h"

class TEvdPanelStrawHitHolder: public TEveElement, public TNamed {
public:
  int              fNumber;                      // panel number

  TEvdPanelStrawHitHolder (int Number = -1);
  ~TEvdPanelStrawHitHolder();
  
  void Clear(Option_t* Opt = ""); 
  
  ClassDef(TEvdPanelStrawHitHolder,0)

};

class TEvdPlaneStrawHitHolder: public TEveElement, public TNamed {
public:
  int                       fNumber;     // plane number
  TEvdPanelStrawHitHolder*  fPanelHitHolder[kNPanels];

  TEvdPlaneStrawHitHolder(int Number = -1);
  ~TEvdPlaneStrawHitHolder();
  
  void Clear(Option_t* Opt = ""); 
  
  ClassDef(TEvdPlaneStrawHitHolder,0)
};


class TEvdStationStrawHitHolder: public TEveElement, public TNamed {
public:
  int                       fNumber;                      // station number
  TEvdPlaneStrawHitHolder*  fPlaneHitHolder[2];

  TEvdStationStrawHitHolder(int Number = -1);
  ~TEvdStationStrawHitHolder();
  
  void Clear(Option_t* Opt = ""); 
  
  ClassDef(TEvdStationStrawHitHolder,0)
};


class TEvdStrawHitHolder : public TEveElementList {
public:
  TEvdStationStrawHitHolder*  fStationHitHolder[kNStations];

  TEvdStrawHitHolder();
  ~TEvdStrawHitHolder();

  void Clear(Option_t* Opt = "");


  TEvdPanelStrawHitHolder* Panel(int Station, int Plane, int Panel) {
    return fStationHitHolder[Station]->fPlaneHitHolder[Plane]->fPanelHitHolder[Panel];
  }
  
  ClassDef(TEvdStrawHitHolder,0)    // Tracker Straw Hits, structured
};

#endif
