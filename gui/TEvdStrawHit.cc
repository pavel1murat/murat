///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////

#include "murat/gui/TEvdStrawHit.hh"
#include "TGeoSphere.h"
#include "TGeoTube.h"
#include "TEveTrans.h"
#include "TGeoMatrix.h"

ClassImp(TEvdStrawHit)

//-----------------------------------------------------------------------------
TEvdStrawHit::TEvdStrawHit(): TEveGeoShape("") {
  //  fNumber = I;
}

//-----------------------------------------------------------------------------
TEvdStrawHit::~TEvdStrawHit() {
}

//-----------------------------------------------------------------------------
// rotation is related to the direction of the view...
//-----------------------------------------------------------------------------
void TEvdStrawHit::Init(int Index, int StrawIndex, int Plane, int Panel, int Layer,
			int Straw, // within the layer
			int DeltaID,
			float Time, float Dt, float EDep, float WDist, float WRes,
			float X, float Y, float Z) {

  fIndex      = Index;
  fStrawIndex = StrawIndex;
  fPlane      = Plane;
  fPanel      = Panel;
  fLayer      = Layer;
  fStraw      = Straw;
  fDeltaID    = DeltaID;
  fTime       = Time;
  fDt         = Dt;
  fEDep       = EDep;
  fWDist      = WDist;
  fWRes       = WRes;
  fX          = X;
  fY          = Y;
  fZ          = Z;
  
  // TGeoHMatrix m;
  // Double_t sscale[3] = { 1., fWRes/2.5, 1. };
  // m.SetScale(sscale);
  // SetTransMatrix(m);
  
  // TGeoSphere* shape = new TGeoSphere(0,2.5);
  // SetShape(shape);

  TGeoTube* shape = new TGeoTube(0,2.5,5.0);
  SetShape(shape);

  SetMainColor(2);
  SetMainTransparency(0);

  fErrorBars = new TEveGeoShape("errbar");
  fErrorBars->SetShape(new TGeoTube(0,0.5,fWRes));
  fErrorBars->SetMainColor(2);
  fErrorBars->SetMainTransparency(70);
}

//-----------------------------------------------------------------------------
void TEvdStrawHit::Print(Option_t* Opt) const {
  printf("hit: index: %5i, strawindex: %5i P:P:L:S %02i:%1i:%1i:%02i time: %8.2f dt: %8.2f edep: %8.3f x: %8.3f y: %8.3f z: %9.3f wdist: %8.3f wres: %8.3f \n",
	 fIndex,fStrawIndex,fPlane,fPanel,fLayer,fStraw,fTime,fDt,fEDep,fX,fY,fZ,fWDist,fWRes);
}
