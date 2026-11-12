//-----------------------------------------------------------------------------
// vis node of the single-panel view
//-----------------------------------------------------------------------------
#include <algorithm>
#include <cmath>
#include <format>
#include <string>

#include "TColor.h"
#include "TGeoBBox.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMatrix.h"
#include "TGeoMedium.h"
#include "TGeoTube.h"
#include "TGeoVolume.h"
#include "TMath.h"

#include "Stntuple/obj/TComboHit.hh"
#include "Stntuple/obj/TComboHitBlock.hh"

#include "murat/gui/TEvdManager.hh"
#include "murat/gui/TEvdTracker.hh"
#include "murat/gui/TComboHitEvd.hh"
#include "murat/gui/TPanelVisNode.hh"

namespace {
//-----------------------------------------------------------------------------
// tracker_v7.txt parameters, the same as in murat/work/test_panel.C
//-----------------------------------------------------------------------------
  const double kDriftVelocity      = 0.060;    // mm/ns, linear drift
  const double kShellHalfThickness = 0.05;     // mm, the hit surface is drawn 100 um thick
  const double kStrawR             = 2.5235;   // tracker.strawOuterRadius
  const double kStrawWallThickness = 0.25;     // display thickness of the wall, 100 um (tracker.strawWallThickness is 15 um)
  const double kStrawGap           = 1.203;    // tracker.strawGap
  const double kRInnermostWire     = 380.;     // tracker.rInnermostWire
  const double kManifoldHL0        = 3.125;    // tracker.manifoldHalfLengths[0]
  const int    kNManifolds         = 48;       // tracker.manifoldsPerEnd

//-----------------------------------------------------------------------------
// the tracker frame in the global Mu2e frame, the same offset as in TEvdTracker,
// TTrkVisNode and TComboHitEvd
//-----------------------------------------------------------------------------
  const double kTrackerX = -3904.;
  const double kTrackerY =     0.;
  const double kTrackerZ = 22935.;

//-----------------------------------------------------------------------------
// the medium with that name if the geometry has it, otherwise a display-only one
//-----------------------------------------------------------------------------
  TGeoMedium* findOrMakeMedium(const char* Name, int Id, double Z, double A, double Density) {
    if (gGeoManager) {
      if (TGeoMedium* medium = gGeoManager->GetMedium(Name)) return medium;
    }
    auto* material = new TGeoMaterial(Name, A, Z, Density);
    return new TGeoMedium(Name, Id, material);
  }

//-----------------------------------------------------------------------------
// radial distance from the tracker axis of the origin of the panel frame
//-----------------------------------------------------------------------------
  double panelOriginRadius() {
    const double rSpacing       = 2.*kStrawR + kStrawGap;
    const double layerHalfShift = 0.25*rSpacing;
    const double yMin = kRInnermostWire - kStrawR;
    const double yMax = kRInnermostWire + 2.*kManifoldHL0*(kNManifolds-1) + 2.*layerHalfShift + kStrawR;
    return 0.5*(yMin + yMax);
  }

//-----------------------------------------------------------------------------
// azimuth and z centre of the panel, from its straws as TEvdTracker places them;
// Flip is true if the layer-0 straws are at larger z than the centre (the panel z-side
// alternates from panel to panel in TrackerMaker); returns false if the panel is not in the tracker
//-----------------------------------------------------------------------------
  bool panelCenter(const mu2e::Tracker* Trk, int Plane, int Panel, double& PhiC, double& Zc, bool& Flip) {
    if (Trk == nullptr) return false;

    const auto& planes = Trk->planes();
    if ((Plane < 0) or (Plane >= (int) planes.size())) return false;
    const mu2e::Plane& plane = planes.at(Plane);
    if ((Panel < 0) or (Panel >= (int) plane.nPanels())) return false;
    const mu2e::Panel& panel = plane.getPanel(Panel);

    double xsum(0.), ysum(0.), zmin(1.e9), zmax(-1.e9);
    double zl0(0.);
    int    nl0(0);
    for (int is=0; is<(int) panel.nStraws(); ++is) {
      const auto& straw = panel.getStraw(is);
      const auto& mid   = straw.getMidPoint();
      xsum += mid.x();
      ysum += mid.y();
      zmin  = std::min(zmin,mid.z());
      zmax  = std::max(zmax,mid.z());
      if (straw.id().getLayer() == 0) {
        zl0 += mid.z();
        ++nl0;
      }
    }

    Zc   = 0.5*(zmax+zmin);
    PhiC = std::atan2(ysum,xsum);
    Flip = (nl0 > 0) and (zl0/nl0 > Zc);
    return true;
  }

//-----------------------------------------------------------------------------
// the panel of the tracker (a TEvdPanel with its straws), placed in the panel frame: X along
// the straw, Y radially outwards, Z perpendicular to the panel; the tracker axis is at the
// origin of the TEvdPanel frame, the panel frame origin is at radius panelOriginRadius()
//-----------------------------------------------------------------------------
  TGeoVolume* makePanelAssembly(const murat::TEvdTracker* Tracker, int Plane, int Panel) {
    static int nAssemblies(0);
    std::string prefix = std::format("panel_{}",nAssemblies++);
    auto* assembly = new TGeoVolumeAssembly(prefix.c_str());
    if (Tracker == nullptr) return assembly;

    const murat::PanelGeometry g = murat::panelGeometry(Tracker,Plane,Panel);
    const mu2e::Panel& panel = Tracker->fTrkPtr->planes().at(Plane).getPanel(Panel);
    TGeoMedium* medium       = gGeoManager->GetMedium("TrackerStraw");

    murat::TEvdPanel* evdPanel = new murat::TEvdPanel((prefix+"_panel").c_str(),g,medium,kGray);
    evdPanel->SetPanelTransparency(95);
    // the straws of the panel view: visible, with their own transparency
    evdPanel->AddStraws(panel,Tracker->fTrkPtr->strawProperties().strawOuterRadius()-kStrawWallThickness,
                        Tracker->fTrkPtr->strawProperties().strawOuterRadius(),g.zc,medium,true,0);

    // from the TEvdPanel frame to the panel frame: translate the axis to -r0 along the panel
    // azimuth, then rotate the azimuth to the panel frame's Y direction
    const double phiC = 0.5*(g.phi1+g.phi2);
    const double r0   = panelOriginRadius();
    TGeoTranslation shift(-r0*std::cos(phiC*TMath::DegToRad()),-r0*std::sin(phiC*TMath::DegToRad()),0.);
    TGeoRotation    rot;
    rot.RotateZ(90.-phiC);

    TGeoHMatrix placement(rot);
    placement.Multiply(&shift);
    assembly->AddNode(evdPanel,1,new TGeoHMatrix(placement));
    return assembly;
  }
}

namespace murat {

//-----------------------------------------------------------------------------
TPanelVisNode::TPanelVisNode(const char* Name) : TEvdVisNode(Name) {
  // build the geometry - that can be done in standalone mode
  fAssembly = nullptr;                        // built for the panel, see Draw

  fListOfHits = new TObjArray();
  fListOfHits->SetOwner(kTRUE);
}

//-----------------------------------------------------------------------------
TPanelVisNode::~TPanelVisNode() {
  fListOfHits->Delete();
  delete fListOfHits;
}

//-----------------------------------------------------------------------------
// all combo hits of the panel: the same plane and panel as the selected hit
//-----------------------------------------------------------------------------
int TPanelVisNode::InitEvent() {
  auto vm = TEvdManager::Instance();
  fChBlock = (TComboHitBlock*) vm->FindDataBlock("ChBlock");
  if (fChBlock == nullptr) return -1;

  f_RunNumber    = fChBlock->RunNumber   ();
  f_EventNumber  = fChBlock->EventNumber ();
  f_SubrunNumber = fChBlock->SubrunNumber();

  fListOfHits->Delete();
  int nch = fChBlock->NHits();
  for (int i=0; i<nch; i++) {
    TComboHit* ch = fChBlock->Hit(i);
    if (ch->Sid().getPlane() != fPlane) continue;
    if (ch->Sid().getPanel() != fPanel) continue;
    // drift radius from the drift time, linear drift (drift time in ns); the hit is a thin
    // cylindrical surface of that radius, concentric with the wire
    double rdrift = std::clamp(kDriftVelocity*ch->DrTime(),0.,kStrawR-kStrawWallThickness);
    fListOfHits->Add(new TComboHitEvd(ch,std::max(0.,rdrift-kShellHalfThickness),
                                      rdrift+kShellHalfThickness,
                                      findOrMakeMedium("Vacuum",1,0.,0.,0.)));
  }
  return 0;
}

//-----------------------------------------------------------------------------
bool TPanelVisNode::Initialized() {
  auto vm = TEvdManager::Instance();
  auto cb = (TComboHitBlock*) vm->FindDataBlock("ChBlock");
  if (cb == nullptr) return false;

  return ((cb->EventNumber () == f_EventNumber ) and
          (cb->SubrunNumber() == f_SubrunNumber) and
          (cb->RunNumber   () == f_RunNumber   ));
}

//-----------------------------------------------------------------------------
// the panel geometry (fAssembly, independent of the tracker geometry), placed with the
// transform the panel has in the tracker, so that the hits line up; and the hits
//-----------------------------------------------------------------------------
void TPanelVisNode::Draw(Option_t* Opt) {
  auto vm = TEvdManager::Instance();
  TEvdTracker* tracker = dynamic_cast<TEvdTracker*>(vm->FindSubdetector("TRACKER"));

  const mu2e::Tracker* trk = (tracker) ? tracker->fTrkPtr.get() : nullptr;

  // the geometry of the selected panel, rebuilt when the panel changes
  if ((fAssembly == nullptr) or (fPlane != fBuiltPlane) or (fPanel != fBuiltPanel)) {
    fAssembly   = makePanelAssembly(tracker,fPlane,fPanel);
    fBuiltPlane = fPlane;
    fBuiltPanel = fPanel;
  }

  static int nviews = 0;
  std::string viewName = std::format("panel_view_{}", nviews++);
  TGeoVolumeAssembly* view = new TGeoVolumeAssembly(viewName.c_str());
  // the view is the panel frame: the panel and the hits are in it, no placement
  view->AddNode(fAssembly, 1);

  // the panel frame in the tracker frame: the hits are expressed in it, so they move with the panel
  double phiC(0.), zc(0.);
  bool   flip(false);
  panelCenter(trk,fPlane,fPanel,phiC,zc,flip);
  const double r0    = panelOriginRadius();
  const double theta = phiC - 0.5*TMath::Pi();      // rotation of the panel frame about z
  const double ct    = std::cos(theta);
  const double st    = std::sin(theta);
  const double ox    = r0*std::cos(phiC);           // origin of the panel frame, tracker frame
  const double oy    = r0*std::sin(phiC);

  // the hits: a red tube along the wire, radius 0.5 * straw radius, half-length the
  // resolution along the wire, centred on the hit position
  for (int i=0; i<fListOfHits->GetEntriesFast(); i++) {
    // the hit volume (a TComboHitEvd) and its combo hit
    TComboHitEvd* hitVolume = (TComboHitEvd*) fListOfHits->UncheckedAt(i);
    TComboHit* ch = hitVolume->fCh;
    if (ch == nullptr) continue;

    // tube axis (z) -> x, then x -> the wire direction in the xy plane
    // the centre of the surface is on the wire: the straw midpoint, moved along the wire to the
    // hit's coordinate along the wire; the hit itself is on the drift circle, off the wire
    double cx(ch->X()), cy(ch->Y()), cz(ch->Z());
    if (trk) {
      const auto& mid = trk->getStraw(ch->Sid()).getMidPoint();
      const double along = (ch->X()-mid.x())*ch->Ux() + (ch->Y()-mid.y())*ch->Uy();
      cx = mid.x() + along*ch->Ux();
      cy = mid.y() + along*ch->Uy();
      cz = mid.z();
    }

    // the centre in the panel frame: origin at (ox,oy,zc), local X along the straw
    const double dx = cx-ox;
    const double dy = cy-oy;
    const double lx =  ct*dx + st*dy;
    const double ly = -st*dx + ct*dy;
    const double lz = cz-zc;

    // the wire direction in the panel frame: the tube axis (z) -> x -> the wire direction
    TGeoRotation wire;
    wire.RotateY(90.);
    wire.RotateZ((std::atan2(ch->Uy(),ch->Ux())-theta)*TMath::RadToDeg());
    TGeoHMatrix local(TGeoTranslation(lx,ly,lz));
    local.Multiply(&wire);
    view->AddNode(hitVolume,i+1,new TGeoHMatrix(local));
  }

  view->Draw(Opt);
}

//-----------------------------------------------------------------------------
// azimuth of the panel (radians) in the tracker, the direction of its local Y axis
//-----------------------------------------------------------------------------
double TPanelVisNode::Azimuth() {
  auto vm = TEvdManager::Instance();
  TEvdTracker* tracker = dynamic_cast<TEvdTracker*>(vm->FindSubdetector("TRACKER"));

  double phiC(0.), zc(0.);
  bool   flip(false);
  panelCenter((tracker) ? tracker->fTrkPtr.get() : nullptr,fPlane,fPanel,phiC,zc,flip);
  return phiC;
}

}
