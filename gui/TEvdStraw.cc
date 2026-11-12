///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////

#include "murat/gui/TEvdStraw.hh"
#include <string>
#include "TGeoTube.h"
#include "TGeoNode.h"
#include "TEveTrans.h"

ClassImp(murat::TEvdStraw)

namespace murat {
//-----------------------------------------------------------------------------
TEvdStraw::TEvdStraw(int I): TEvdSubdetector() {
  fNumber = I;
}

//-----------------------------------------------------------------------------
TEvdStraw::TEvdStraw(const char* Name, double InnerRadius, double OuterRadius,
                     double HalfLength, TGeoMedium* Medium, bool Visible, int Transparency):
  TEvdSubdetector(Name, new TGeoTube((std::string(Name)+"_tube").c_str(), InnerRadius, OuterRadius, HalfLength),
                  Medium, kGray+1, Transparency)
{
  fNumber     = -1;
  fID         = -1;
  fPlane      = -1;
  fPanel      = -1;
  fLayer      = -1;
  fZ          = 0.;
  fRho        = 0.;
  fNx         = 0.;
  fNy         = 0.;
  fHalfLength = HalfLength;
  // the solid is the volume of the first node: hide that one, unless the straw is visible
  GetNode(0)->GetVolume()->SetVisibility(Visible);
  SetVisibility(Visible);   // the straw assembly too
}

//-----------------------------------------------------------------------------
TEvdStraw::~TEvdStraw() {
}

//-----------------------------------------------------------------------------
// rotation is related to the direction of the view...
//-----------------------------------------------------------------------------
void TEvdStraw::Init(int ID, int Plane, int Panel, int Layer, double Rho, double Z, double nx, double ny, double HalfLength) {

  fID         = ID;
  fPlane      = Plane;
  fPanel      = Panel;
  fLayer      = Layer;
  fZ          = Z;
  fRho        = Rho;
  fNx         = nx;
  fNy         = ny;
  fHalfLength = HalfLength;

  // TGeoTube* tube = new TGeoTube(2.45,2.5,fHalfLength);
  // SetShape(tube);
}

//-----------------------------------------------------------------------------
void TEvdStraw::Print(Option_t* Opt) const {
  printf("straw: fIndex: %5i P:P:L:S = %02i:%1i:%1i:%02i rho:%8.2f halfL:%8.2f nx:%8.4f ny: %8.4f \n",
	 fID,fPlane,fPanel,fLayer,fNumber,fRho,fHalfLength,fNx,fNy);
}

}
