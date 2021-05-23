///////////////////////////////////////////////////////////////////////////////

#include "murat/gui/TEvdTrack.hh"

ClassImp(TEvdTrack)

//-----------------------------------------------------------------------------
TEvdTrack::TEvdTrack(const char* Name, double* X0, double* V0, double W, double* ZRange):
THelix(X0,V0,W) , 
  fName(Name)
{
  SetLineWidth(2);
  SetRange(ZRange);
}

//-----------------------------------------------------------------------------
const char* TEvdTrack::GetName() const {
  return fName.Data();
}
