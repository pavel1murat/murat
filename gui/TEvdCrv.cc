///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include <format>

#include "TGeoBBox.h"

#include "Offline/ConfigTools/inc/SimpleConfig.hh"

#include "murat/gui/TEvdCrv.hh"

#include "Offline/GeometryService/inc/CosmicRayShieldMaker.hh"

#include "murat/gui/TEvdManager.hh"

ClassImp(murat::TEvdCrv)

//-----------------------------------------------------------------------------
namespace {
  TGeoMedium* crvMedium() {
    // The manager must already exist and be current.
    auto* gm = gGeoManager;
    
    auto* medium = gm->GetMedium("CRVScintillator");
    if (medium) {
      return medium;
    }
    
    // Display-only material.  Replace with the desired optical/material
    // properties if required.
    auto* material = new TGeoMaterial("CRVScintillator",
                                      1.032,   // density
                                      12.0,    // effective Z
                                      0.0);    // radiation length placeholder
    
    return new TGeoMedium("CRVScintillator", 100, material);
  }
}

namespace murat {
  //-----------------------------------------------------------------------------
  TEvdCrv::TEvdCrv(): TEvdSubdetector() {
    // The assembly is the top volume owned by TEvdCrv.
    fTopVolume = new TGeoVolumeAssembly("CRV");
    fName = "CRV";
  }

  //-----------------------------------------------------------------------------
  TEvdCrv::TEvdCrv(const char* Fn): TEvdSubdetector() {
    // The assembly is the top volume owned by TEvdCrv.
    fTopVolume = new TGeoVolumeAssembly("CRV");
    fName = "CRV";
    InitGeometry(Fn);
  }
  
  //-----------------------------------------------------------------------------
  TEvdCrv::~TEvdCrv() {
  }

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
// int TEvdCrv::InitGeometry(const char* Fn) {
//   mu2e::CosmicRayShieldMaker crv_maker(mu2e::SimpleConfig(Fn),0);
//   fCrvPtr    = crv_maker.getCosmicRayShieldPtr();

//   fNSectors = fCrvPtr->getCRSScintillatorShields().size();
  
//   // and initialize the sectors

//   auto vm = TEvdManager::Instance(); // a

//   for (int i=0; i<fNSectors; i++) {
//     TEvdCrvSector* sector = Sector(i);
//     vm->AddSubdetector(sector);
//   }
  
//   return 0;
// }

//-----------------------------------------------------------------------------
int TEvdCrv::InitGeometry(const char* geomFile) {
  if (!gGeoManager) {
    std::cerr << "TEvdCrv::InitGeometry: no current TGeoManager\n";
    return 1;
  }

  if (fCrvPtr) {
    return 0;                         // already initialized
  }

  mu2e::SimpleConfig config(geomFile);
  mu2e::CosmicRayShieldMaker maker(config, 0);

  fCrvPtr = maker.getCosmicRayShieldPtr();
  if (!fCrvPtr) {
    std::cerr << "TEvdCrv::InitGeometry: CRV geometry is null\n";
    return 1;
  }

  auto* medium = crvMedium();

  const auto& shields = fCrvPtr->getCRSScintillatorShields();
  fNSectors = static_cast<int>(shields.size());

  for (int is = 0; is < fNSectors; ++is) {
    auto* sector = new TEvdCrvSector;

    const auto& shield = shields.at(is);
    
    const std::string sectorName = shield.getName();

    auto* sectorAssembly = new TGeoVolumeAssembly(sectorName.c_str());

    // Add all scintillator bars.  Positions in the Mu2e geometry are global,
    // so the sector assembly is placed at the origin.
    for (const auto& module : shield.getCRSScintillatorModules()) {
      for (const auto& layer : module.getLayers()) {
        for (const auto& barPtr : layer.getBars()) {
          const auto& bar = *barPtr;

          const auto& h = bar.getHalfLengths();
          if (h.size() != 3) {
            std::cerr << "Invalid half-length vector for CRV bar\n";
            return 1;
          }
          
          const std::string barName = std::format("CRV_{}_{}_{}_{}_{}",
                                                  sectorName,
                                                  bar.id().getShieldNumber(),
                                                  bar.id().getModuleNumber(),
                                                  bar.id().getLayerNumber (),
                                                  bar.id().getBarNumber   ());

          auto* shape   = new TGeoBBox       (barName.c_str(),h[0],h[1],h[2]);
          auto* counter = new TEvdSubdetector(barName.c_str(), shape, medium);

          // teh bar color will depend on whether the bar has hits
          counter->SetLineColor(kCyan + 1);
          counter->SetFillColor(kCyan + 1);
          counter->SetTransparency(0);

          const auto& p = bar.getPosition();

          // CRSScintillatorBar::getHalfLengths() is in world x/y/z order,
          // hence no rotation is needed here.
          sectorAssembly->AddNode(counter,
                                  bar.id().getBarNumber() + 1,
                                  new TGeoTranslation(p.x(), p.y(), p.z()));
        }
      }
    }

    // // Optional mechanical components.
    // for (const auto& module : shield.getCRSScintillatorModules()) {
    //   for (const auto& sheet : module.getAluminumSheets()) {
    //     const auto& h = sheet.getHalfLengths();
    //     const auto& p = sheet.getPosition();

    //     auto* shape = new TGeoBBox("CRVAluminumSheetShape",
    //                                h[0], h[1], h[2]);

    //     // For a display-only geometry, use the scintillator medium unless
    //     // dedicated aluminum media are defined.
    //     auto* volume = new TGeoVolume("CRVAluminumSheet",
    //                                    shape, medium);

    //     sectorAssembly->AddNode(
    //       volume,
    //       1,
    //       new TGeoTranslation(p.x(), p.y(), p.z()));
    //   }

    //   for (const auto& absorber : module.getAbsorberLayers()) {
    //     const auto& h = absorber.getHalfLengths();
    //     const auto& p = absorber.getPosition();

    //     auto* shape = new TGeoBBox("CRVAbsorberShape",
    //                                h[0], h[1], h[2]);

    //     auto* volume = new TGeoVolume("CRVAbsorber",
    //                                    shape, medium);

    //     sectorAssembly->AddNode(
    //       volume,
    //       1,
    //       new TGeoTranslation(p.x(), p.y(), p.z()));
    //   }
    // }
    
    fTopVolume->AddNode(sectorAssembly, is + 1, new TGeoTranslation);

    // Keep the wrapper object alive through the CRV object.
    fListOfSubdetectors->Add(sector);
  }

  fCopyNumber = 1;

  return 0;
}


//-----------------------------------------------------------------------------
void TEvdCrv::Print(Option_t* Opt) const {
  printf("TEvdCrv: emoe\n");
}

}
