///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <format>

#include <algorithm>
#include <array>
#include "TMath.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoTube.h"
#include "TGeoBoolNode.h"
#include "TGeoCompositeShape.h"
#include "TGeoXtru.h"
#include "murat/gui/TEvdTracker.hh"
#include "murat/gui/TEvdPanel.hh"

#include "Offline/TrackerGeom/inc/Plane.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/GeometryService/inc/TrackerMaker.hh"

#include "murat/gui/TEvdSubdetector.hh"
#include "murat/gui/TEvdManager.hh"

ClassImp(murat::TEvdTracker)

namespace {

TGeoMedium* trackerStrawMedium() {
  auto* gm = gGeoManager;

  if (auto* medium = gm->GetMedium("TrackerStraw")) {
    return medium;
  }

  // Display-only material.
  auto* material = new TGeoMaterial("TrackerStrawMaterial",
                                    1.0,    // density
                                    6.0,    // effective Z
                                    1.0);   // radiation length placeholder

  return new TGeoMedium("TrackerStraw", 200, material);
}
  /*
//-----------------------------------------------------------------------------
std::unique_ptr<TGeoCombiTrans> strawTransform(const mu2e::Straw& straw) {
  const auto& p = straw.origin();
  const auto  d = straw.direction().unit();

  // TGeoTube is oriented along its local z axis.  Construct a rotation
  // which maps local z onto the straw direction.
  const double theta = std::acos(
      std::max(-1.0, std::min(1.0, d.z())));

  const double phi = std::atan2(d.y(), d.x());

  auto* rotation = new TGeoRotation;
  rotation->RotateY(theta * 180.0 / M_PI);
  rotation->RotateZ(phi   * 180.0 / M_PI);

  return std::make_unique<TGeoCombiTrans>(
      p.x(), p.y(), p.z(), rotation);
}
  */
} // namespace



namespace murat {

//-----------------------------------------------------------------------------
TEvdTracker::TEvdTracker(): TEvdSubdetector("TRACKER") {
  fName       = "TRACKER";
  // fTopVolume  = this;
  fCopyNumber    = 1;
  fDisplayPanels = 0;
}

//-----------------------------------------------------------------------------
  TEvdPanel* TEvdTracker::Panel(int Station, int Plane, int Panel) {
    return fStation[Station]->fPlane[Plane]->fPanel[Panel];
  }

//-----------------------------------------------------------------------------
int TEvdTracker::InitGeometry(const char* geomFile) {
  std::cout << std::format("-- START: TEvdTracker::InitGeometry\n");
  if (!gGeoManager) {
    std::cerr << "TEvdTracker::InitGeometry: "
              << "no current TGeoManager\n";
    return 1;
  }

  if (fTrkPtr) {
    return 0;                         // already initialized
  }

  // Build the Mu2e tracker data geometry.
  mu2e::SimpleConfig config(geomFile);
  mu2e::TrackerMaker maker(config);

  fTrkPtr = maker.getTrackerPtr();

  if (!fTrkPtr) {
    std::cerr << "TEvdTracker::InitGeometry: "
              << "tracker geometry is null\n";
    return 1;
  }

  // the straw medium is created here, the planes use it
  trackerStrawMedium();

  // each plane builds its own geometry: the plane volume, its panels and their placement
  int nplanes = fTrkPtr->planes().size();
  for (int ip = 0; ip < nplanes; ++ip) {
    new TEvdPlane(this,ip);
  }

  std::cout << std::format("-- END: TEvdTracker::InitGeometry\n");
  return 0;
}


  //-----------------------------------------------------------------------------
  int TEvdTracker::InitEvent() {
    //    auto vm = TEvdManager::Instance();

    if (Initialized()) return 0;

    SetInitialized();
  
    return 0;
  }

  //-----------------------------------------------------------------------------
  void TEvdTracker::Print(Option_t* Opt) const {
    std::string opt(Opt);
    if (opt == "") {
      std::cout << std::format("TEvdTracker::Print\n");
    }
  }
}
