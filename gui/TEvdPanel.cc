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
#include "TGeoManager.h"

ClassImp(murat::TEvdPanel)

namespace murat {
//-----------------------------------------------------------------------------
namespace {
//-----------------------------------------------------------------------------
// sector of the disk minus the triangle with the origin and the two sector ends as vertices
//-----------------------------------------------------------------------------
TGeoShape* panelShape(const char* Name, double Rin, double Rout, double Dz, double Phi1, double Phi2) {
  std::string name(Name);
  TGeoTubeSeg* sector = new TGeoTubeSeg((name+"_sector").c_str(), Rin, Rout, Dz, Phi1, Phi2);

  double xv[3] = {0., Rout*std::cos(Phi1*TMath::DegToRad()), Rout*std::cos(Phi2*TMath::DegToRad())};
  double yv[3] = {0., Rout*std::sin(Phi1*TMath::DegToRad()), Rout*std::sin(Phi2*TMath::DegToRad())};
  TGeoXtru* triangle = new TGeoXtru(2);
  triangle->DefinePolygon(3, xv, yv);
  triangle->DefineSection(0, -Dz-1.);   // extended beyond the sector: no coincident faces
  triangle->DefineSection(1,  Dz+1.);

  TGeoSubtraction* sub = new TGeoSubtraction(sector, triangle, nullptr, nullptr);
  return new TGeoCompositeShape(name.c_str(), sub);
}

//-----------------------------------------------------------------------------
std::string panelName(int Plane, int Panel) {
  return std::format("TrackerPlane_{}_Panel_{}", Plane, Panel);
}

}


//-----------------------------------------------------------------------------
// geometry of the panel volume, from straws 0, 1 and 2 of the panel as TEvdTracker places them:
// straw 0 is the innermost, longest straw, at rInnermostWire; straw 1 is the neighbour in the
// other layer, straw 2 the next straw of layer 0; the layers are symmetric about the panel
// centre, and all straws of a panel lie on one ray (the panel azimuth)
//-----------------------------------------------------------------------------
PanelGeometry panelGeometry(const TEvdTracker* Tracker, int Plane, int Panel) {
  const double kPanelPhi = 120.;                  // degrees, as the tracker support parts
  const double rOut      = Tracker->fTrkPtr->strawProperties().strawOuterRadius();

  const mu2e::Panel& panel = Tracker->fTrkPtr->planes().at(Plane).getPanel(Panel);
  const auto& s0 = panel.getStraw(0).getMidPoint();
  const auto& s1 = panel.getStraw(1).getMidPoint();
  const auto& s2 = panel.getStraw(2).getMidPoint();

  double r0           = std::hypot(s0.x(),s0.y());
  double pitch        = std::hypot(s2.x(),s2.y()) - r0;          // straw spacing within a layer
  double shift        = std::hypot(s1.x(),s1.y()) - r0;          // offset of layer 1
  const int nPerLayer = 96/2;                                    // 96 straws per panel

  double phiC = std::atan2(s0.y(),s0.x())*TMath::RadToDeg();

  PanelGeometry g;
  g.rin  = r0 - rOut;                                              // innermost straw - straw radius
  g.rout = r0 + (nPerLayer-1)*pitch + std::max(shift,0.) + rOut;   // outermost straw + straw radius
  g.dz   = 0.5*std::abs(s1.z()-s0.z()) + rOut;                     // thin in z
  g.phi1 = phiC - 0.5*kPanelPhi;
  g.phi2 = phiC + 0.5*kPanelPhi;
  g.zc   = 0.5*(s0.z()+s1.z());
  return g;
}

//-----------------------------------------------------------------------------
TEvdPanel::TEvdPanel(TEvdTracker* Tracker, int Plane, int Panel, bool CreateStraws):
  TEvdPanel(panelName(Plane,Panel).c_str(), panelGeometry(Tracker,Plane,Panel),
            gGeoManager->GetMedium("TrackerStraw"), kGray)
{
  // tracker local frame -> parent frame, as in TTrkVisNode::trackPoints
  const double kDX = -3904.;
  const double kDZ = 22935.;

  const PanelGeometry g = panelGeometry(Tracker,Plane,Panel);
  TGeoMedium* medium    = gGeoManager->GetMedium("TrackerStraw");
  const double rOut     = Tracker->fTrkPtr->strawProperties().strawOuterRadius();

  Tracker->AddNode(this,Plane*100+Panel+1,new TGeoTranslation(kDX,0.,kDZ+g.zc));
  if (CreateStraws) {
    const mu2e::Panel& panel = Tracker->fTrkPtr->planes().at(Plane).getPanel(Panel);
    AddStraws(panel, Tracker->fTrkPtr->strawProperties().strawInnerRadius(), rOut, g.zc, medium);
  }
  Tracker->fListOfSubdetectors->Add(this);
}

//-----------------------------------------------------------------------------
TEvdPanel::TEvdPanel(const char* Name, const PanelGeometry& Geometry,
                     TGeoMedium* Medium, int Color):
  TEvdSubdetector(Name, panelShape(Name,Geometry.rin,Geometry.rout,Geometry.dz,Geometry.phi1,Geometry.phi2),
                  Medium, Color, 0)
{
  fNumber = -1;
  fPlane  = nullptr;
  fNx     = 0.;
  fNy     = 0.;
  fPhi    = 0.;
}


//-----------------------------------------------------------------------------
// the transparency is set on the solid, the first node of the panel volume
//-----------------------------------------------------------------------------
void TEvdPanel::SetPanelTransparency(int Transparency) {
  GetNode(0)->GetVolume()->ResetTransparency(Transparency);   // per volume: the medium is shared
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
                          double ZCenter, TGeoMedium* Medium, bool Visible, int Transparency) {
  std::string panelName(GetName());
  for (int is = 0; is < (int) Panel.nStraws(); ++is) {
    const mu2e::Straw& straw = Panel.getStraw(is);
    std::string name = std::format("{}_straw_{}", panelName, is);

    murat::TEvdStraw* evd_straw = new murat::TEvdStraw(name.c_str(), StrawInnerRadius, StrawOuterRadius,
                                                       straw.halfLength(), Medium, Visible, Transparency);
    auto m = straw.getMidPoint();
    auto d = straw.strawDirection();
    const auto& sid = straw.id();
    evd_straw->Init(sid.asUint16(), sid.getPlane(), sid.getPanel(), sid.getLayer(),
                    std::hypot(m.x(),m.y()), m.z(), d.x(), d.y(), straw.halfLength());

    // the tube axis follows the straw direction
    double theta = std::acos(std::clamp(d.z(), -1., 1.))*TMath::RadToDeg();
    double phi   = std::atan2(d.y(), d.x())*TMath::RadToDeg();
    // Rz(phi)*Ry(theta) maps the tube axis (z) onto the direction (sin theta cos phi, sin theta sin phi, cos theta)
    TGeoRotation* rot = new TGeoRotation();
    rot->RotateY(theta);
    rot->RotateZ(phi);

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
