///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "murat/gui/TEvdStraw.hh"
#include "murat/gui/TEvdPanel.hh"
#include "murat/gui/TEvdTracker.hh"
#include "TGeoXtru.h"
#include "TEveTrans.h"

ClassImp(TEvdPanel)

//-----------------------------------------------------------------------------
TEvdPanel::TEvdPanel(int I): TEveGeoShape() {
  fNumber = I;
  SetElementName(Form("panel_%03i",I));
  for (int i=0; i<kNStraws; i++) {
    fStraw[i] = new TEvdStraw(i);
    AddElement(fStraw[i]);
  }

  //-----------------------------------------------------------------------------
  // define dimensions
  //-----------------------------------------------------------------------------
  const int npt(33);

  double yy[npt] = {
    370.000, 380.000, 390.000, 400.000, 410.000, 420.000, 430.000, 440.000, 450.000, 460.000,
    470.000, 480.000, 490.000, 500.000, 510.000, 520.000, 530.000, 540.000, 550.000, 560.000,
    570.000, 580.000, 590.000, 600.000, 610.000, 620.000, 630.000, 640.000, 650.000, 660.000,
    670.000, 680.000, 690.000
  };

  double xx[npt] = {
    605.970, 599.750, 593.296, 586.600, 579.655, 572.451, 564.978, 557.225, 549.181, 540.833,
    532.165, 523.163, 513.809, 504.083, 493.964, 483.425, 472.440, 460.977, 448.999, 436.463,
    423.320, 409.512, 394.968, 379.605, 363.318, 345.977, 327.414, 307.409, 285.657, 261.725,
    234.947, 204.206, 167.332
  };

  double zx[2*npt], zy[2*npt];

  for (int i=0; i<npt; i++) { zy[i]     = yy[i]      ; zx[i]     = -xx[i]; }
  for (int i=0; i<npt; i++) { zy[npt+i] = yy[npt-i-1]; zx[npt+i] = xx[npt-i-1]; }

  TGeoXtru* xtru = new TGeoXtru(2);

  xtru->DefinePolygon(2*npt,zx,zy);
 
  xtru->DefineSection(0,-10.,0,0,1.);
  xtru->DefineSection(1, 10.,0,0,1.);

  SetShape(xtru);
}


//-----------------------------------------------------------------------------
// needs to be called after  TEvemanager has been initialized
// at this point the wire coordinates are already known
//-----------------------------------------------------------------------------
void TEvdPanel::InitGeometry() {

  double zpanel = (fStraw[0]->Z() + fStraw[1]->Z())/2;
  this->RefMainTrans().SetPos(0,0,zpanel);
  this->RefMainTrans().RotatePF(1,2,fPhi);

  SetMainTransparency(95);
  SetMainColor(kGray+1);
  //-----------------------------------------------------------------------------
  // by default, straws are drawn along the Z axis
  // panels are orthogonal to the Z axis, synchronize
  //-----------------------------------------------------------------------------
  for (int i=0; i<kNStraws; i++) {
    TEvdStraw* s = fStraw[i];
   
    double dz = s->Z()-zpanel;
    s->RefMainTrans().SetPos(s->Rho(),0,dz);
    
    s->RefMainTrans().RotateLF(2,3,TMath::Pi()/2);
    s->RefMainTrans().RotatePF(1,2,TMath::Pi()/2);

    s->SetMainTransparency(80);
    s->SetMainColor(kRed-9);
  }

  SetRnrSelfChildren(true,false);
}


//-----------------------------------------------------------------------------
void TEvdPanel::InitStraw(int straw, int StrawID, int Plane, int Panel, int Layer,
			  double rho, double z,
			  double nx, double ny, double half_length)
{
  fStraw[straw]->Init(StrawID,Plane,Panel,Layer,rho,z,nx,ny,half_length);
}

//-----------------------------------------------------------------------------
TEvdPanel::~TEvdPanel() {
}


//-----------------------------------------------------------------------------
void TEvdPanel::Print(Option_t* Opt) const {
  printf("panel: plane: %02i index: %1i nx: %8.5f ny: %8.5f phi: %8.5f z: %9.3f\n",
	 fPlane->Number(),fNumber,fNx,fNy,fPhi,Z());
}
