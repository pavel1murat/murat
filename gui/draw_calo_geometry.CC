//
#include "Offline/GeometryService/inc/DiskCalorimeterMaker.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Offline/CalorimeterGeom/inc/Disk.hh"
#include "Offline/CalorimeterGeom/inc/DiskInfo.hh"
#include "Offline/CalorimeterGeom/inc/Crystal.hh"

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

#include "murat/gui/draw_calo_geometry.hh"

ClassImp(murat::TEvdCrystal)
ClassImp(murat::TEvdDisk)
ClassImp(murat::TEvdCalorimeter)
ClassImp(murat::TEvdSubdetector)
ClassImp(murat::TEvdManagerA)
// ============================================================================
// Global state used by the interactive callback
// ============================================================================
// std::unique_ptr<mu2e::DiskCalorimeter> gCalo = nullptr;

static const mu2e::DiskCalorimeter*    gCaloForPicking = nullptr;
static TCanvas*                        gCaloCanvas = nullptr;

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
  
static bool decodeCaloPath(const char* path, int& diskID, int& crystalCopy) {
  diskID      = -1;
  crystalCopy = -1;

  if (!path) return false;

  std::string p(path);
  std::size_t begin = 0;

  while (begin < p.size()) {
    std::size_t end = p.find('/', begin);
    if (end == std::string::npos) end = p.size();

    const std::string token = p.substr(begin, end - begin);

    // Disk node: Disk_<disk ID>_<copy number>
    if (token.rfind("Disk_", 0) == 0) {
      const std::size_t first = 5;
      const std::size_t last = token.find('_', first);

      if (last != std::string::npos) {
        diskID = std::atoi(token.substr(first, last - first).c_str());
      }
    }

    // Crystal node: Crystal_<copy number>
    if (token.rfind("Crystal_", 0) == 0) {
      const std::size_t first = 8;
      crystalCopy = std::atoi(token.substr(first).c_str());
    }

    begin = end + 1;
  }

  return diskID >= 0 && crystalCopy > 0;
}

// ============================================================================
// Print the selected crystal
// ============================================================================

static void printSelectedCrystal(int diskID, int crystalCopy) {

  std::cout << std::format("-- START printSelectedCrystal\n");
  
  if (!gCaloForPicking) return;

  if (diskID < 0 || diskID >= static_cast<int>(gCaloForPicking->nDisks())) {
    std::cerr << "Invalid disk ID: " << diskID << '\n';
    return;
  }

  const mu2e::Disk& disk = gCaloForPicking->disk(diskID);

  // Crystal nodes are created with copy number localID + 1.
  const int localID = crystalCopy - 1;

  if ((localID < 0) or (localID >= static_cast<int>(disk.nCrystals()))) {
    std::cerr << "Invalid crystal copy number: "
              << crystalCopy << " for disk " << diskID << '\n';
    return;
  }

  const mu2e::Crystal& crystal = disk.crystal(localID);

  std::cout << "\n========================================\n";
  std::cout << "Selected calorimeter crystal\n";
  std::cout << "disk ID       : " << diskID << '\n';
  std::cout << "local crystal : " << localID << '\n';
  std::cout << "global ID     : " << crystal.ID() << '\n';
  std::cout << "========================================\n";

  /*
   * If mu2e::Crystal derives from TObject and implements Print(), use:

       crystal.Print();

   * For example, if Crystal is declared as:

       class Crystal : public TObject

   * then the following works:
   */
  // const TObject* rootCrystal = dynamic_cast<const TObject*>(&crystal);

  // if (rootCrystal) {
  //   rootCrystal->Print();
  // } else {
    /*
     * The Crystal class shown in the question currently does not derive
     * from TObject.  This fallback is useful until it does.
     */
    std::cout << "local position : "
              << crystal.localPosition() << '\n';
    std::cout << "global position: "
              << crystal.position() << '\n';
    std::cout << "size           : "
              << crystal.size() << '\n';
    //  }

  std::cout << "========================================\n\n";
}

// ============================================================================
// Canvas event callback
//
// Right mouse button pressed in the OpenGL canvas invokes this function.
// ============================================================================
void Mu2eCaloProcessedEvent(Int_t event, Int_t px, Int_t py, TObject*) {
  
  std::cout << std::format("-- START Mu2eCaloProcessedEvent event:{}\n",event);
  
  if (event != kButton3Down) return;
  if (!gCaloForPicking     ) return;
  if (!gGeoManager         ) return;

  /*
   * In the OpenGL viewer, ROOT updates the current TGeo node during
   * mouse picking.  GetPath() returns the selected node path.
   */
  TGeoNode* node = gGeoManager->GetCurrentNode();

  std::cout << std::format("-- Mu2eCaloProcessedEvent node:{}\n",(void*) node);
                                                    
  if (!node) {
    std::cout << "[calo] No geometry node selected\n";
    return;
  }

  const char* path = gGeoManager->GetPath();

  if (!path) {
    std::cout << "[calo] No geometry path available\n";
    return;
  }
  std::cout << std::format("-- Mu2eCaloProcessedEvent path:{}\n",path);

  int diskID = -1;
  int crystalCopy = -1;

  if ( !decodeCaloPath(path, diskID, crystalCopy) ) {
    std::cout << "[calo] Selected node is not a crystal\n";
    std::cout << "[calo] Path: " << path << '\n';
    return;
  }

  std::cout << "[calo] Selected path: " << path << '\n';

  printSelectedCrystal(diskID, crystalCopy);

  (void)px;
  (void)py;
}

// ============================================================================
// Main entry point
// The DiskCalorimeter must already be initialized.
// ============================================================================
void TEvdCalorimeter::Display(bool drawCrystals, int maxCrystalsPerDisk, bool drawXYProjections) {
  fTop->Draw("ogl");
  gSystem->ProcessEvents();
}

//-----------------------------------------------------------------------------
TEvdCalorimeter::TEvdCalorimeter(const char* Fn) {
  // initialize the calorimeter geometry, get ready to display
  // int rc(0);
  
// --------------------------------------------------------------------------
// offline calorimeter
// --------------------------------------------------------------------------
  // auto config = mu2e::SimpleConfig(Fn);
  mu2e::DiskCalorimeterMaker calo_maker(mu2e::SimpleConfig(Fn),0);
  
  auto calo       = calo_maker.calorimeterPtr();
  fCaloPtr        = std::move(calo);
  gCaloForPicking = fCaloPtr.get();
//-----------------------------------------------------------------------------
// tracker and CRV
//-----------------------------------------------------------------------------
  mu2e::CosmicRayShieldMaker crv_maker(mu2e::SimpleConfig(Fn),0);
  fCrvPtr    = crv_maker.getCosmicRayShieldPtr();

  mu2e::SimpleConfig config(Fn);
  
  mu2e::TrackerMaker*  trk_maker = new mu2e::TrackerMaker(config);
  fTrkPtr     = trk_maker->getTrackerPtr();
// --------------------------------------------------------------------------
// Determine world dimensions
// --------------------------------------------------------------------------
  double maxRadius = 0.0;
  double maxZ      = 0.0;
  int ndisks       = fCaloPtr->nDisks();
  
  for (int idisk=0; idisk<ndisks; ++idisk) {
    const auto& info = fCaloPtr->disk(idisk).diskInfo();
    maxRadius        = std::max(maxRadius, info.outerEnvelopeR());
    maxZ             = std::max(maxZ, std::abs(info.origin().z()) + info.size().z());
  }

  std::cout << std::format("maxRadius:{:10} maxZ:{:10}\n",maxRadius,maxZ);
  if (maxRadius <= 0.0) maxRadius = 1000.0;
  if (maxZ      <= 0.0) maxZ      = 3000.0;
// --------------------------------------------------------------------------
// ROOT materials and media
// --------------------------------------------------------------------------
  fGeoManager = new TGeoManager("mu2e_calo", "Mu2e Geometry");

  auto vacuumMaterial  = new TGeoMaterial("Vacuum", 0.0, 0.0, 0.0);
  auto vacuumMedium    = new TGeoMedium  ("Vacuum", 1, vacuumMaterial);

  auto diskMaterial    = new TGeoMaterial("DiskMaterial", 26.98, 13.0, 2.70);
  auto diskMedium      = new TGeoMedium  ("DiskMedium"  , 2, diskMaterial);

  auto crystalMaterial = new TGeoMaterial("CrystalMaterial", 259.81, 54.0, 4.51);
  auto crystalMedium   = new TGeoMedium  ("CrystalMedium"  , 3, crystalMaterial);
// --------------------------------------------------------------------------
// World
// --------------------------------------------------------------------------
  fTop = fGeoManager->MakeBox("TOP",vacuumMedium,1.5 * maxRadius,1.5 * maxRadius,1.5 * maxZ);
  fGeoManager->SetTopVolume(fTop);

// --------------------------------------------------------------------------
// Crystal volume
//
// A single reusable volume is sufficient.
// Each AddNode call creates an individual crystal placement.
// --------------------------------------------------------------------------
  const auto& crystal = fCaloPtr->disk(0).crystal(0);

  auto shape = new TGeoBBox("CrystalShape",
                            0.5 * crystal.size().x()-0.1,
                            0.5 * crystal.size().y()-0.1,
                            0.5 * crystal.size().z()
                            );

  auto crystalVolume = new TEvdCrystal("Crystal",shape,crystalMedium);

  crystalVolume->SetLineColor(kOrange + 1);
  crystalVolume->SetFillColor(kOrange + 1);
  crystalVolume->SetTransparency(0);
// --------------------------------------------------------------------------
// Draw disks and crystals
// --------------------------------------------------------------------------
  for (unsigned idisk=0; idisk<fCaloPtr->nDisks(); ++idisk) {
    const mu2e::Disk& disk = fCaloPtr->disk(idisk);
    const auto& info       = disk.diskInfo();

    const double halfThickness = 0.5 * info.size().z();

    auto diskShape = new TGeoTube(Form("DiskShape_%u", idisk),
                                  info.innerEnvelopeR(),
                                  info.outerEnvelopeR(),
                                  halfThickness
                                  );

    auto diskVolume = new TEvdDisk(Form("Disk_%u", idisk),diskShape,diskMedium);

    diskVolume->SetLineColor(kAzure + 2);
    diskVolume->SetFillColor(kAzure + 2);
    diskVolume->SetTransparency(50);

    /*
     * DiskInfo::toGlobal(local) is:

         global = origin + inverseRotation * local

     * Therefore inverseRotation is the local-to-global rotation.
     */
    const auto& rotation = info.inverseRotation();

    Double_t matrix[9] = {
      rotation.xx(), rotation.xy(), rotation.xz(),
      rotation.yx(), rotation.yy(), rotation.yz(),
      rotation.zx(), rotation.zy(), rotation.zz()
    };

    auto geoRotation = new TGeoRotation;
    geoRotation->SetName(Form("DiskRotation_%u", idisk));
    geoRotation->SetMatrix(matrix);

    auto diskTransform = new TGeoCombiTrans(
      info.origin().x(),
      info.origin().y(),
      info.origin().z(),
      geoRotation
    );

    // Disk copy number is idisk + 1.
    fTop->AddNode(diskVolume, idisk + 1,diskTransform);

    int nPlaced = 0;
    int ncr = disk.nCrystals();
    
    for (int i=0; i<ncr; ++i) {
      const auto& crystal = disk.crystal(i);

      /*
       * Crystal::localPosition() is the front-face position. 
         - how do we know that?
       * TGeoBBox is centered on its local origin, so shift by +z/2.
       */
      
      const double x = crystal.localPosition().x();
      const double y = crystal.localPosition().y();
      const double z = crystal.localPosition().z() + 0.5 * crystal.size().z();

      std::cout << std::format("crystal i:{:4d} x:{:10.3f} y:{:10.3f} z:{:10.3f}\n",i,x,y,z);
//-----------------------------------------------------------------------------
// Copy number is local crystal ID + 1.  This is decoded by the right-click callback.
//-----------------------------------------------------------------------------
      diskVolume->AddNode(crystalVolume,i+1,new TGeoTranslation(x,y,z));
      ++nPlaced;
    }
    std::cout << std::format("idisk:{} nPlaced:{:4d}\n",idisk,nPlaced);
  }
// --------------------------------------------------------------------------
// Close and display OpenGL geometry
// --------------------------------------------------------------------------
  fGeoManager->CloseGeometry();
  //  return rc;
}

//-----------------------------------------------------------------------------
void TEvdCrystal::Print(Option_t* Opt) const {
  std::cout << std::format("emoe\n");
}

//-----------------------------------------------------------------------------
void TEvdCrystal::PrintA() const {
  std::cout << std::format("emoe AAAAAAA\n");
}
  
//-----------------------------------------------------------------------------
void TEvdDisk::Print(Option_t* Opt) const {
  std::cout << std::format("TEvdDisk: emoe AAAAAAA\n");
}
  
}
