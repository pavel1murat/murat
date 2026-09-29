///////////////////////////////////////////////////////////////////////////////

#include "murat/gui/TEvdTrack.hh"

ClassImp(TEvdTrack)

//-----------------------------------------------------------------------------
TEvdTrack::TEvdTrack(const char* Name, double* X0, double* V0, double W, double* ZRange):
TPolyLine3D(),
  fName(Name)
{
  SetLineWidth(2);

  // figure out the points... set polyline
  // either : SetPoint(Int_t point, Double_t x, Double_t y, Double_t z)
  // or     : SetPolyLine(Int_t n, Float_t *p, Option_t *option = "");
}

//-----------------------------------------------------------------------------
const char* TEvdTrack::GetName() const {
  return fName.Data();
}
