///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////

#include "murat/gui/TEvdHelix.hh"
#include "TVector3.h"

ClassImp(TEvdHelix)

//-----------------------------------------------------------------------------
TEvdHelix::TEvdHelix(): TEveLine("helix") {
}

//-----------------------------------------------------------------------------
// omega = 1/R (signed),
// phi0 gives the direction of the particle in the point of closest approach
// D0 : Y coordinate of the point of closest approach in a rotated reference
//      frame with Y axis pointing towards the center of the helix in XY plane
// tandip : Pz/Pxy
//-----------------------------------------------------------------------------
TEvdHelix::TEvdHelix(double Z0, double D0, double Phi0, double Omega, double TanDip,
		     double ZMin, double ZMax): TEveLine("helix")
{
  fD0     = D0;
  fPhi0   = Phi0;
  fTanDip = TanDip;
  fOmega  = Omega;
  fZ0     = Z0;
  fZMin   = ZMin;
  fZMax   = ZMax;

  fX0     =  (1/fOmega+fD0)*sin(fPhi0);
  fY0     =  -(1/fOmega+fD0)*cos(fPhi0);

  double dphi = (ZMax-ZMin)*Omega/TanDip;
//-----------------------------------------------------------------------------
// angular step between points - 0.05
//-----------------------------------------------------------------------------
  int npt = dphi/0.05+1;
  TVector3 pos;
  for (int i=0; i<npt; i++) {
    double z = ZMin+i*(ZMax-ZMin)/(npt-1);
    GetPointAtZ(z,&pos);
    SetNextPoint(pos.x(),pos.y(),pos.z());
  }
}

//-----------------------------------------------------------------------------
TEvdHelix::~TEvdHelix() {
}

//-----------------------------------------------------------------------------
// point at a given S from Z0
//-----------------------------------------------------------------------------
void TEvdHelix::GetPointAtS(double S, TVector3* Point) {

  double cdip = CosDip();
  double sdip = fTanDip*cdip;

  double phi  = fPhi0 + cdip*S*fOmega;
  double cphi = cos(phi);
  double sphi = sin(phi);
  // double sphi0 = sin(fPhi0);
  // double cphi0 = cos(fPhi0);

  // double x     =  (sphi - sphi0)/fOmega - fD0*sphi0;
  // double y     = -(cphi - cphi0)/fOmega + fD0*cphi0;

  double x     = fX0+sphi/fOmega;
  double y     = fY0-cphi/fOmega;
  double z     = fZ0+S*sdip;

  Point->SetXYZ(x,y,z);
}

//-----------------------------------------------------------------------------
// point at a given Z
//-----------------------------------------------------------------------------
void TEvdHelix::GetPointAtZ(double Z, TVector3* Point) {
  double cdip  = CosDip();
  double sdip  = fTanDip*cdip;
  double ds    = (Z-fZ0)/sdip;
  GetPointAtS(ds,Point);
}

//-----------------------------------------------------------------------------
void TEvdHelix::Print(Option_t* Opt) const {
  // printf("straw: fIndex: %5i P:P:L:S = %02i:%1i:%1i:%02i rho:%8.2f halfL:%8.2f nx:%8.4f ny: %8.4f \n",
  // 	 fID,fPlane,fPanel,fLayer,fNumber,fRho,fHalfLength,fNx,fNy);
}
