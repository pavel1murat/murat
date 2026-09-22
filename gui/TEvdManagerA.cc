///////////////////////////////////////////////////////////////////////////////
// calorimeter has its mother and places two disks into it
//-----------------------------------------------------------------------------
#include "Offline/GeometryService/inc/CosmicRayShieldMaker.hh"
#include "Offline/GeometryService/inc/TrackerMaker.hh"

#include <TCanvas.h>
#include <TGeoBBox.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TGeoMatrix.h>
#include <TGeoNode.h>
#include <TGeoTube.h>
#include <TGeoVolume.h>
#include <TH2F.h>
#include <TObject.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TVirtualPad.h>
#include <Buttons.h>
#include <TString.h>

#include <format>
#include <cstdlib>
#include <iostream>
#include <string>

#include "TGeoMatrix.h"

#include "murat/gui/TEvdManagerA.hh"
#include "murat/gui/TEvdCrv.hh"
#include "murat/gui/TEvdCalorimeter.hh"
#include "murat/gui/TEvdTracker.hh"

ClassImp(murat::TEvdManagerA)
// ============================================================================
// Global state used by the interactive callback
// ============================================================================
// std::unique_ptr<mu2e::DiskCalorimeter> gCalo = nullptr;

// ============================================================================
// Extract disk ID and crystal copy number from a TGeo path
//
// Expected path, for example:
//
//   /TOP_1/Disk_0_1/Crystal_37
//
// Disk copy number is not used because the disk ID is encoded in Disk_<id>.
// Crystal copy number is the local crystal index + 1.
// ============================================================================
namespace murat {

// ============================================================================
// Main entry point
// The DiskCalorimeter must already be initialized.
// ============================================================================
void TEvdManagerA::AddSubdetector(TEvdSubdetector* Sd) {
  // at this point internal subdetector tree of Sd is already built
  
  fListOfSubdetectors.Add(Sd);

  TGeoVolume* vol = Sd->TopVolume();
  if (vol) {
    auto top = fGeoManager->GetTopNode();
    // all detectors are different ... this place needs work
    
    top->GetVolume()->AddNode(vol,Sd->CopyNumber(),new TGeoTranslation());
  }

}

//-----------------------------------------------------------------------------
void TEvdManagerA::Draw(Option_t* Opt) {
  auto top = fGeoManager->GetTopNode();  // minimize geometry knowledge by Vis manger 
  top->Draw("ogl");
  
  gSystem->ProcessEvents();
}

//-----------------------------------------------------------------------------
// initial configuration - in FCL
//-----------------------------------------------------------------------------
TEvdManagerA::TEvdManagerA(const char* Fn) {

  //  auto gm = TGeoManagerA::Instance(); // defines its own TOP

  fGeoManager          = new TGeoManager("mu2e_calo", "Mu2e Geometry");

  auto vacuum_material = new TGeoMaterial("Vacuum", 0.0, 0.0, 0.0);

  fGeoManager->AddMaterial(vacuum_material);
  //  fGeoManager->AddMedium  (vacuumMedium);
// --------------------------------------------------------------------------
// World
// --------------------------------------------------------------------------
  auto vacuum_medium = new TGeoMedium  ("Vacuum", 1, vacuum_material);
  TGeoVolume* top = fGeoManager->MakeBox("TOP",vacuum_medium,30000.,30000.,30000);
  fGeoManager->SetTopVolume(top);

  if (fDisplayCalorimeter) {
    auto calo = new TEvdCalorimeter("murat/fcl/geom_extracted.txt"); // includes geometry initialization
    AddSubdetector(calo);
  }
  if (fDisplayCrv) {
    auto crv = new TEvdCrv("murat/fcl/geom_extracted.txt");
    AddSubdetector(crv);
  }
// --------------------------------------------------------------------------
// Close geometry, ready to display
// --------------------------------------------------------------------------
  fGeoManager->CloseGeometry();
  //  return rc;
}

//-----------------------------------------------------------------------------
TEvdManagerA* TEvdManagerA::Instance(const char* Fcl) {
  static TEvdManagerA* instance(nullptr);
  if (instance == nullptr) {
    instance = new TEvdManagerA(Fcl);
  }
  return instance;
}

}
