///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "murat/gui/TEvdStraw.hh"
#include "murat/gui/TEvdPanel.hh"
#include "murat/gui/TEvdTracker.hh"
#include "TGeoXtru.h"
#include "TGeoTube.h"
#include <algorithm>
#include <format>
#include "TGeoBoolNode.h"
#include "TGeoMatrix.h"
#include "Offline/TrackerGeom/inc/Panel.hh"
#include "Offline/TrackerGeom/inc/Straw.hh"
#include "TGeoCompositeShape.h"
#include "TMath.h"

ClassImp(murat::TEvdPanel)

namespace murat {
//-----------------------------------------------------------------------------
namespace {
//-----------------------------------------------------------------------------
// sector of the disk minus the triangle with the origin and the two sector ends as vertices
//-----------------------------------------------------------------------------
TGeoShape* panelShape(const char* Name, double Rout, double Dz, double Phi1, double Phi2) {
  std::string name(Name);
  TGeoTubeSeg* sector = new TGeoTubeSeg((name+"_sector").c_str(), 0., Rout, Dz, Phi1, Phi2);

  double xv[3] = {0., Rout*std::cos(Phi1*TMath::DegToRad()), Rout*std::cos(Phi2*TMath::DegToRad())};
  double yv[3] = {0., Rout*std::sin(Phi1*TMath::DegToRad()), Rout*std::sin(Phi2*TMath::DegToRad())};
  TGeoXtru* triangle = new TGeoXtru(2);
  triangle->DefinePolygon(3, xv, yv);
  triangle->DefineSection(0, -Dz-1.);   // extended beyond the sector: no coincident faces
  triangle->DefineSection(1,  Dz+1.);

  TGeoSubtraction* sub = new TGeoSubtraction(sector, triangle, nullptr, nullptr);
  return new TGeoCompositeShape(name.c_str(), sub);
}
}

//-----------------------------------------------------------------------------
TEvdPanel::TEvdPanel(const char* Name, double Rout, double Dz, double Phi1, double Phi2,
                     TGeoMedium* Medium, int Color, int Transparency):
  TEvdSubdetector(Name, panelShape(Name,Rout,Dz,Phi1,Phi2), Medium, Color, Transparency)
{
  fNumber = -1;
  fPlane  = nullptr;
  fNx     = 0.;
  fNy     = 0.;
  fPhi    = 0.;
}


//-----------------------------------------------------------------------------
// needs to be called after  TEvemanager has been initialized
// at this point the wire coordinates are already known
//-----------------------------------------------------------------------------
int TEvdPanel::InitGeometry(const char* Fn) {

  // double zpanel = (fStraw[0]->Z() + fStraw[1]->Z())/2;
  // this->RefMainTrans().SetPos(0,0,zpanel);
  // this->RefMainTrans().RotatePF(1,2,fPhi);

  // SetMainTransparency(95);
  // SetMainColor(kGray+1);
  // //-----------------------------------------------------------------------------
  // // by default, straws are drawn along the Z axis
  // // panels are orthogonal to the Z axis, synchronize
  // //-----------------------------------------------------------------------------
  // for (int i=0; i<kNStraws; i++) {
  //   TEvdStraw* s = fStraw[i];
   
  //   double dz = s->Z()-zpanel;
  //   s->RefMainTrans().SetPos(s->Rho(),0,dz);
    
  //   s->RefMainTrans().RotateLF(2,3,TMath::Pi()/2);
  //   s->RefMainTrans().RotatePF(1,2,TMath::Pi()/2);

  //   s->SetMainTransparency(80);
  //   s->SetMainColor(kRed-9);
  // }

  // SetRnrSelfChildren(true,false);
  return 0;
}


//-----------------------------------------------------------------------------
void TEvdPanel::InitStraw(int straw, int StrawID, int Plane, int Panel, int Layer,
			  double rho, double z,
			  double nx, double ny, double half_length)
{
  fStraw[straw]->Init(StrawID,Plane,Panel,Layer,rho,z,nx,ny,half_length);
}

//-----------------------------------------------------------------------------
void TEvdPanel::AddStraws(const mu2e::Panel& Panel, double StrawInnerRadius, double StrawOuterRadius,
                          double ZCenter, TGeoMedium* Medium) {
  std::string panelName(GetName());
  for (int is = 0; is < (int) Panel.nStraws(); ++is) {
    const mu2e::Straw& straw = Panel.getStraw(is);
    std::string name = std::format("{}_straw_{}", panelName, is);

    murat::TEvdStraw* evd_straw = new murat::TEvdStraw(name.c_str(), StrawInnerRadius, StrawOuterRadius,
                                                       straw.halfLength(), Medium);
    auto m = straw.getMidPoint();
    auto d = straw.strawDirection();
    const auto& sid = straw.id();
    evd_straw->Init(sid.asUint16(), sid.getPlane(), sid.getPanel(), sid.getLayer(),
                    std::hypot(m.x(),m.y()), m.z(), d.x(), d.y(), straw.halfLength());

    // the tube axis follows the straw direction
    double theta = std::acos(std::clamp(d.z(), -1., 1.))*TMath::RadToDeg();
    double phi   = std::atan2(d.y(), d.x())*TMath::RadToDeg();
    TGeoRotation* rot = new TGeoRotation(name.c_str(), phi, theta, 0.);

    AddNode(evd_straw, is+1, new TGeoCombiTrans(m.x(), m.y(), m.z()-ZCenter, rot));
  }
}

//-----------------------------------------------------------------------------
TEvdPanel::~TEvdPanel() {
}


//-----------------------------------------------------------------------------
void TEvdPanel::Print(Option_t* Opt) const {
  printf("panel: plane: %02i index: %1i nx: %8.5f ny: %8.5f phi: %8.5f z: %9.3f\n",
	 fPlane->Number(),fNumber,fNx,fNy,fPhi,Z());
}

}
