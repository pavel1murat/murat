///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////

#include "murat/gui/TEvdStraw.hh"
#include "TGeoTube.h"
#include "TEveTrans.h"

ClassImp(TEvdStraw)

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
void TEvdStraw::Init(int Index, double Rho, double Z, double nx, double ny, double HalfLength) {

  fIndex      = Index;
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
  printf("straw : fNumber : %5i,  fIndex: %5i\n",fNumber,fIndex);
}
