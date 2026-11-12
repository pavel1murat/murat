///////////////////////////////////////////////////////////////////////////////
// a tracker plane: the plane volume and its panels, placed in the tracker
///////////////////////////////////////////////////////////////////////////////

#include <format>
#include <string>

#include "TColor.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoTube.h"
#include "TMath.h"

#include "murat/gui/TEvdPlane.hh"
#include "murat/gui/TEvdTracker.hh"

ClassImp(murat::TEvdPlane)

namespace murat {
//-----------------------------------------------------------------------------
TEvdPlane::TEvdPlane(TEvdTracker* Tracker, int PlaneNumber):
  TEvdSubdetector(std::format("TrackerPlane_{}",PlaneNumber).c_str(),
                  new TGeoTube(0.,700.,10.),gGeoManager->GetMedium("TrackerStraw"),kAzure+1,95),
  fNumber(PlaneNumber),
  fTracker(Tracker)
{
  // tracker local frame -> parent frame, as in TTrkVisNode::trackPoints
  const double kDX = -3904.;
  const auto& plane = fTracker->fTrkPtr->planes().at(PlaneNumber);

  // z range of the plane: the straws of panels 0 and 1, plus half a straw
  const mu2e::Panel& p0 = plane.getPanel(0);
  double z0 = p0.getStraw(0).getMidPoint().z();
  double z1 = p0.getStraw(1).getMidPoint().z();
  double zmin = z0-0.5;
  double zmax = z1+0.5;
  if (z1 < z0) {
    zmin = z1 - 0.5;
    zmax = z0 - 0.5;
  }

  const mu2e::Panel& p1 = plane.getPanel(1);
  double z01 = p1.getStraw(0).getMidPoint().z();
  double z11 = p1.getStraw(1).getMidPoint().z();
  double zmin1 = z01-0.5;
  double zmax1 = z11+0.5;
  if (z11 < z01) {
    zmin1 = z11 - 0.5;
    zmax1 = z01 - 0.5;
  }

  if (zmax1 > zmax) {
    zmax = zmax1;
  }
  else {
    zmin = zmin1;
  }

  // the plane in the tracker
  const double z0Tracker = fTracker->fTrkPtr->g4Tracker()->z0();
  fTracker->AddNode(this,PlaneNumber+1,new TGeoTranslation(kDX,0.,z0Tracker+(zmin+zmax)/2.));
  fTracker->fListOfSubdetectors->Add(this);

  // individual panel volumes: each panel places itself in the tracker
  if (fTracker->fDisplayPanels == 1) {
    for (int ipa = 0; ipa < (int) plane.nPanels(); ++ipa) {
      fPanel[ipa] = new TEvdPanel(fTracker,PlaneNumber,ipa,fTracker->fDisplayStraws == 1);
      fPanel[ipa]->SetPanelTransparency(fTracker->fPanelTransparency);
    }
  }
}

}
