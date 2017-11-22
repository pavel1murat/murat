///////////////////////////////////////////////////////////////////////////////

#include "TGeoTube.h"
#include "murat/gui/TEvdTracker.hh"

ClassImp(TEvdTracker)

//-----------------------------------------------------------------------------
TEvdTracker::TEvdTracker(): TEveElementList ("tracker") {
  for (int i=0; i<kNStations; i++) {
    fStation[i] = new TEvdStation(i);
    AddElement(fStation[i]);
  }
  
  // TGeoTube* tube = new TGeoTube(300,700,1550);
  // SetShape(tube);
  
  SetRnrSelfChildren(false,true);
}


