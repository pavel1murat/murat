///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <format>

#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoTube.h"
#include "murat/gui/TEvdTracker.hh"

#include "Offline/TrackerGeom/inc/Plane.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/GeometryService/inc/TrackerMaker.hh"

#include "murat/gui/TEvdSubdetector.hh"

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

} // namespace



namespace murat {

//-----------------------------------------------------------------------------
TEvdTracker::TEvdTracker(): TEvdSubdetector("TRACKER") {
  fName       = "TRACKER";
  // fTopVolume  = this;
  fCopyNumber = 1;
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

  auto tg4 = fTrkPtr->g4Tracker();

  if (!fTrkPtr) {
    std::cerr << "TEvdTracker::InitGeometry: "
              << "tracker geometry is null\n";
    return 1;
  }

  auto* medium = trackerStrawMedium();

  /*
   * TEvdTracker itself should have been constructed as a logical
   * TEvdSubdetector volume, for example:
   *
   *   TEvdTracker::TEvdTracker()
   *     : TEvdSubdetector("Tracker") {}
   *
   * Therefore all tracker children are attached directly to this object.
   */

  const auto& planes = fTrkPtr->planes();

  for (std::size_t ip = 0; ip < planes.size(); ++ip) {
    const mu2e::Plane& plane = planes.at(ip);

    const std::string planeName = std::format("TrackerPlane_{}", ip);

    const CLHEP::Hep3Vector& plane_origin  = plane.origin();

    double plane_radius = 700.;
    double plane_dz2    = 10.;
    
    auto* sd_plane = new TEvdSubdetector(planeName.c_str(),
                                         new TGeoTube(0,plane_radius,plane_dz2),medium,kAzure+1,50);

    const mu2e::Panel* panel = &plane.getPanel(0);
    
    double z0 = panel->getStraw(0).getMidPoint().z();
    double z1 = panel->getStraw(1).getMidPoint().z();
    double zmin = z0-0.5;
    double zmax = z1+0.5;
    if (z1 < z0) {
      zmin = z1 - 0.5;
      zmax = z0 - 0.5;
    }

    const mu2e::Panel* p1 = &plane.getPanel(1);
    
    double z01 = p1->getStraw(0).getMidPoint().z();
    double z11 = p1->getStraw(1).getMidPoint().z();
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

    // const auto& panels = plane.panels();

    // for (std::size_t ipa = 0; ipa < panels.size(); ++ipa) {
    //   const auto* panel = panels.at(ipa);

    //   if (!panel) {
    //     continue;
    //   }

    //   const std::string panelName = std::format("TrackerPlane_{}_Panel_{}", ip, ipa);

    //   auto* panelVolume = new TEvdSubdetector(panelName.c_str());

    //   panelVolume->SetLineColor(kGreen + 1);
    //   panelVolume->SetFillColor(kGreen + 1);
    //   panelVolume->SetTransparency(40);

    //   // const auto& straws = panel->straws();

    //   // for (std::size_t ist = 0; ist < straws.size(); ++ist) {
    //   //   const auto* straw = straws.at(ist);

    //   //   if (!straw) {
    //   //     continue;
    //   //   }

    //   //   const auto& props = fTrkPtr->strawProperties();

    //   //   const std::string strawName = std::format("TrackerStraw_{}_{}_{}",ip, ipa, ist);

    //   //   auto* strawShape = new TGeoTube(strawName.c_str(),
    //   //                                   props.strawInnerRadius(),
    //   //                                   props.strawOuterRadius(),
    //   //                                   straw->halfLength());

    //   //   auto* strawVolume = new TEvdSubdetector(strawName.c_str(), strawShape, medium);

    //   //   strawVolume->SetLineColor(kGray + 1);
    //   //   strawVolume->SetFillColor(kGray + 1);
    //   //   strawVolume->SetTransparency(0);

    //   //   auto transform = strawTransform(*straw);

    //   //   panelVolume->AddNode(
    //   //     strawVolume,
    //   //     static_cast<int>(ist) + 1,
    //   //     transform.release());
    //   // }

    //   // panelVolume->AddSubdetectorChildrenToListIfNeeded();
    //   planeVolume->AddNode(panelVolume, static_cast<int>(ipa) + 1, new TGeoTranslation());

    //   fListOfSubdetectors->Add(panelVolume);
    // }

    //planeVolume->AddSubdetectorChildrenToListIfNeeded();
    AddNode(sd_plane, static_cast<int>(ip) + 1, new TGeoTranslation(-3904.,0,tg4->z0()+(zmin+zmax)/2.));

    fListOfSubdetectors->Add(sd_plane);
  }

  std::cout << std::format("-- END: TEvdTracker::InitGeometry\n");
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
