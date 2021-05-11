///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////

#include "murat/gui/TEvdStraw.hh"
#include "TGeoTube.h"
#include "TEveTrans.h"

ClassImp(murat::TEvdStraw)

namespace murat {
//-----------------------------------------------------------------------------
TEvdStraw::TEvdStraw(int I): TEveGeoShape("a") {
  fNumber = I;
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

  TGeoTube* tube = new TGeoTube(2.45,2.5,fHalfLength);
  SetShape(tube);
}

//-----------------------------------------------------------------------------
void TEvdStraw::Print(Option_t* Opt) const {
  printf("straw: fIndex: %5i P:P:L:S = %02i:%1i:%1i:%02i rho:%8.2f halfL:%8.2f nx:%8.4f ny: %8.4f \n",
	 fID,fPlane,fPanel,fLayer,fNumber,fRho,fHalfLength,fNx,fNy);
}

}
