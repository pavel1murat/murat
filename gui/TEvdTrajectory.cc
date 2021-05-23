///////////////////////////////////////////////////////////////////////////////

#include "murat/gui/TEvdTrajectory.hh"

ClassImp(TEvdTrajectory)

//-----------------------------------------------------------------------------
TEvdTrajectory::TEvdTrajectory(const char* Name, int NPoints, float* X, float* Y, float* Z)
: TPolyLine3D(NPoints,X,Y,Z), 
  fName(Name)
{
  SetLineColor(kRed+1);
  SetLineWidth(1);
}

//-----------------------------------------------------------------------------
const char* TEvdTrajectory::GetName() const {
  return fName.Data();
}
