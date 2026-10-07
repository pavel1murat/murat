///////////////////////////////////////////////////////////////////////////////
// calorimeter has its mother and places two disks into it
//-----------------------------------------------------------------------------
#include "cetlib/filepath_maker.h"
#include "fhiclcpp/ParameterSet.h"

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
#include "TString.h"
#include "TEnv.h"

#include <format>
#include <cstdlib>
#include <iostream>
#include <string>

#include "TGeoMatrix.h"

#include "murat/gui/TEvdManager.hh"
#include "murat/gui/TEvdCrv.hh"
#include "murat/gui/TEvdCalorimeter.hh"
#include "murat/gui/TEvdTracker.hh"

ClassImp(murat::TEvdManager)
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
void TEvdManager::AddSubdetector(TEvdSubdetector* Sd) {
  // at this point internal subdetector tree of Sd is already built
  
  auto top = fGeoManager->GetTopNode();
  // all detectors are different ... this place needs work
  // the subdetector is positioned in a global reference frame
  top->GetVolume()->AddNode(Sd,Sd->CopyNumber(),new TGeoTranslation());

  fListOfSubdetectors->Add(Sd);

}

// //-----------------------------------------------------------------------------
// void TEvdManager::Draw(Option_t* Opt) {
//   auto top = fGeoManager->GetTopNode();  // minimize geometry knowledge by Vis manger 
//   top->Draw("ogl");
  
//   gSystem->ProcessEvents();
// }

//-----------------------------------------------------------------------------
// initial configuration - in FCL
// the name of the config file - in .rootrc (EvdManager.ConfigFile)
//-----------------------------------------------------------------------------
TEvdManager::TEvdManager() {

  std::string config_fn = gEnv->GetValue("EvdManager.ConfigFcl","murat/fcl/evd_config.fcl");

  cet::filepath_lookup policy("FHICL_FILE_PATH");

  auto const pset    = fhicl::ParameterSet::make(config_fn,policy);

  auto evd_config = pset.get<fhicl::ParameterSet>("evd_config");

  fGeometryFile       = evd_config.get<std::string>("geometryFile");
  fDisplayCalorimeter = evd_config.get<bool>       ("displayCalorimeter");
  fDisplayCrv         = evd_config.get<bool>       ("displayCrv");
  fDisplayTracker     = evd_config.get<bool>       ("displayTracker");

  fListOfViews        = new TObjArray();
  fListOfNodes        = new TObjArray();
  fListOfSubdetectors = new TObjArray();
  fListOfDataBlocks   = new TMap();
  fListOfDataBlocks->SetOwnerKeyValue(kTRUE, kFALSE);   // owns only keys, not values
  
  fGeoManager         = new TGeoManager("mu2e_geo", "Mu2e Geometry");
  
  auto vacuum_material = new TGeoMaterial("Vacuum", 0.0, 0.0, 0.0);

  // fGeoManager->AddMaterial(vacuum_material);
  //  fGeoManager->AddMedium  (vacuumMedium); // cant add a medium...
// --------------------------------------------------------------------------
// World
// --------------------------------------------------------------------------
  auto *vacuum_medium = new TGeoMedium("Vacuum"    , 1, vacuum_material);
  auto *world_shape   = new TGeoBBox  ("WorldShape",30000.,30000.,30000);
  TGeoVolume *world   = new TGeoVolume("World"     , world_shape, vacuum_medium);
  
  fGeoManager->SetTopVolume(world);
}
  
//-----------------------------------------------------------------------------
// nodes have unique names, each of them could be a non-trivial object
//-----------------------------------------------------------------------------
int TEvdManager::AddNode(TEvdVisNode* Node) {
  int rc(0);
  int n = fListOfNodes->GetEntriesFast();

  bool found = false;
  for (int i=0; i<n; i++) {
    TEvdVisNode* node = GetNode(i);
    if (strcmp(Node->GetName(),node->GetName()) == 0) {
      found = true;
      break;
    }
  }

  if (not found) {
    fListOfNodes->Add(Node);
  }
  
  return rc;
}

//-----------------------------------------------------------------------------
// views are unique, make sure not adding a view the second time
//-----------------------------------------------------------------------------
int TEvdManager::AddView(TEvdView* View) {
  int  rc   (0);
  bool view_found(false);
  
  int nviews = GetNViews();
  for (int i=0; i<nviews; i++) {
    TEvdView* v = GetView(i);
    if (v == View) {
      view_found = true;
      break;
    }
  }
  if (view_found) return rc;
  // new view
  fListOfViews->Add(View); 
  
  int nvn = View->GetNNodes();

  int nnodes = GetNNodes();
  for (int i=0; i<nvn; i++) {
    TEvdVisNode* n = View->GetNode(i);
    
    bool node_found = false;
    
    for (int j=0; j<nnodes; j++) {
      TEvdVisNode* node = GetNode(j);
      if (node == n) {
        node_found = true;
        break;
      }
    }

    if (not node_found) {
      // new node
      fListOfNodes->Add(n);
    }
  }
  
  return rc;
}
  
//-----------------------------------------------------------------------------
// views are unique, make sure not adding a view the second time
//-----------------------------------------------------------------------------
int TEvdManager::DisplayEvent() {
  int rc(0);
  
  // 1. update all nodes with the current event data
  int nnodes = GetNNodes();
  for (int i=0; i<nnodes; i++) {
    TEvdVisNode* node = GetNode(i);
    node->InitEvent();
  }

  // 2. redraw all currently open views

  int nviews = GetNViews();
  for (int i=0; i<nviews; i++) {
    TEvdView* view = GetView(i);
    if (view->IsOpen()) {
      view->Update();
    }
  }
  
  return rc;
}

//-----------------------------------------------------------------------------
// nodes have unique names, each of them could be a non-trivial object
//-----------------------------------------------------------------------------
TEvdSubdetector* TEvdManager::FindSubdetector(const char* Name) {

  TEvdSubdetector* found(nullptr);
  
  int n = fListOfSubdetectors->GetEntriesFast();

  for (int i=0; i<n; i++) {
    TEvdSubdetector* sd = GetSubdetector(i);
    if (strcmp(sd->GetName(),Name) == 0) {
      found = sd;
      break;
    }
  }

  return found;
}

//-----------------------------------------------------------------------------
// each view has enough information to initialize itself
//-----------------------------------------------------------------------------
  int TEvdManager::InitEvent() {
  int rc(0);
  int nn = GetNNodes();
  for (int i=0; i<nn; i++) {
    TEvdVisNode* node = GetNode(i);
    //    if (node->Initialized()) continue;
    node->InitEvent();
  }
  return rc;
}



//-----------------------------------------------------------------------------
int TEvdManager::InitGeometry() {

  if (fDisplayCalorimeter) {
    auto calo = new TEvdCalorimeter();
    calo->InitGeometry(fGeometryFile.data()); // includes geometry initialization
    AddSubdetector(calo);
  }
  if (fDisplayCrv) {
    auto crv = new TEvdCrv();
    crv->InitGeometry(fGeometryFile.data());
    AddSubdetector(crv);
  }
  if (fDisplayTracker) {
    auto tracker = new TEvdTracker();
    tracker->InitGeometry(fGeometryFile.data());
    AddSubdetector(tracker);
  }
// --------------------------------------------------------------------------
// Close geometry, ready to display
// --------------------------------------------------------------------------
  fGeoManager->CloseGeometry();
  return 0;
}

//-----------------------------------------------------------------------------
TEvdManager* TEvdManager::Instance() {
  static TEvdManager* instance(nullptr);
  if (instance == nullptr) {
    instance = new TEvdManager();
  }
  return instance;
}

}
