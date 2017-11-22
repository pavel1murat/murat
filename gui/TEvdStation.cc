///////////////////////////////////////////////////////////////////////////////

#include "murat/gui/TEvdStation.hh"

ClassImp(TEvdStation)

//-----------------------------------------------------------------------------
TEvdStation::TEvdStation(int I): TEveElementList(Form("station_%02i",I),Form("station_%02i",I)) {
  fNumber = I;
  for (int i=0; i<2; i++) {
    fPlane[i] = new TEvdPlane(i);
    AddElement(fPlane[i]);
  }
  SetRnrSelfChildren(false,true);
}

