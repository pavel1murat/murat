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

    // const CLHEP::Hep3Vector& plane_origin  = plane.origin();

    double plane_radius = 700.;
    double plane_dz2    = 10.;
    
    auto* sd_plane = new TEvdSubdetector(planeName.c_str(),
                                         new TGeoTube(0,plane_radius,plane_dz2),medium,kAzure+1,95);

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

    // individual panel volumes: thin sections of tubes spanning 120 degrees in phi,
    // as the tracker support parts in TrackerMaker (phiCoverage = 120)
    if (fDisplayPanels == 1) {
      // tracker local frame -> parent frame, as in TTrkVisNode::trackPoints
      const double kDX = -3904.;
      const double kDZ = 22935.;
      const double kPanelPhi = 120.;                  // degrees
      double rOut = fTrkPtr->strawProperties().strawOuterRadius();

      for (int ipa = 0; ipa < (int) plane.nPanels(); ++ipa) {
        const mu2e::Panel& panel = plane.getPanel(ipa);
        if (panel.nStraws() == 0) continue;

        // radial, z and azimuthal extent of the straws of the panel
        double rmin = 1.e9, rmax = -1.e9, zmin = 1.e9, zmax = -1.e9;
        double xsum = 0., ysum = 0.;
        for (int is = 0; is < (int) panel.nStraws(); ++is) {
          const auto& straw = panel.getStraw(is);
          double x = straw.getMidPoint().x();
          double y = straw.getMidPoint().y();
          double z = straw.getMidPoint().z();
          double r = std::hypot(x,y);
          rmin = std::min(rmin,r);  rmax = std::max(rmax,r);
          zmin = std::min(zmin,z);  zmax = std::max(zmax,z);
          xsum += x; ysum += y;
        }

        double rout    = rmax + rOut;
        double dz      = 0.5*(zmax-zmin) + rOut;          // thin in z
        double zc      = 0.5*(zmax+zmin);
        double phiC    = std::atan2(ysum,xsum)*TMath::RadToDeg();
        double phi1    = phiC - 0.5*kPanelPhi;
        double phi2    = phiC + 0.5*kPanelPhi;

        std::string panelName = std::format("TrackerPlane_{}_Panel_{}", ip, ipa);
        murat::TEvdPanel* panelVolume = new murat::TEvdPanel(panelName.c_str(), rout, dz, phi1, phi2,
                                                             medium, kGray, fPanelTransparency);
        AddNode(panelVolume, static_cast<int>(ip*100+ipa+1), new TGeoTranslation(kDX, 0., kDZ + zc));
        if (fDisplayStraws == 1) {
          panelVolume->AddStraws(panel, fTrkPtr->strawProperties().strawInnerRadius(), rOut, zc, medium);
        }
        fListOfSubdetectors->Add(panelVolume);
      }
    }

    AddNode(sd_plane, static_cast<int>(ip) + 1, new TGeoTranslation(-3904.,0,tg4->z0()+(zmin+zmax)/2.));

    fListOfSubdetectors->Add(sd_plane);
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
