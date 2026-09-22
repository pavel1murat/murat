///////////////////////////////////////////////////////////////////////////////
#include "Offline/GeometryService/inc/DiskCalorimeterMaker.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Offline/CalorimeterGeom/inc/Disk.hh"
#include "Offline/CalorimeterGeom/inc/DiskInfo.hh"
#include "Offline/CalorimeterGeom/inc/Crystal.hh"

#include "murat/gui/TEvdManagerA.hh"
#include "murat/gui/TEvdCalorimeter.hh"
#include "TGeoBBox.h"
#include "TGeoTube.h"
#include "root/TGeoMatrix.h"

#include <format>

ClassImp(murat::TEvdCalorimeter)

namespace murat {
//-----------------------------------------------------------------------------
TEvdCalorimeter::TEvdCalorimeter(): TEvdSubdetector() {
}

//-----------------------------------------------------------------------------
TEvdCalorimeter::TEvdCalorimeter(const char* Fn): TEvdSubdetector() {
}

//-----------------------------------------------------------------------------
int TEvdCalorimeter::InitGeometry(const char* Fn) {
  int rc(0);
  
  auto vm = TEvdManagerA::Instance();
// --------------------------------------------------------------------------
// offline calorimeter
// --------------------------------------------------------------------------
  // auto config = mu2e::SimpleConfig(Fn);
  mu2e::DiskCalorimeterMaker calo_maker(mu2e::SimpleConfig(Fn),0);
  
  auto calo       = calo_maker.calorimeterPtr();
  fCaloPtr        = std::move(calo);
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
  auto gm     = vm->GetGeoManager();
  //   auto vacuum = gm->GetMedium("Vacuum");

  auto diskMaterial    = new TGeoMaterial("DiskMaterial", 26.98, 13.0, 2.70);
  auto diskMedium      = new TGeoMedium  ("DiskMedium"  , 2, diskMaterial);

  auto crystalMaterial = new TGeoMaterial("CrystalMaterial", 259.81, 54.0, 4.51);
  auto crystalMedium   = new TGeoMedium  ("CrystalMedium"  , 3, crystalMaterial);
// --------------------------------------------------------------------------
// Crystal volume
//
// A single reusable volume is sufficient.
// Each AddNode call creates an individual crystal placement.
// all crystals are the same
// --------------------------------------------------------------------------
  const mu2e::Crystal* mu2e_cr0 = &fCaloPtr->disk(0).crystal(0);

  auto shape = new TGeoBBox("CrystalShape",
                            0.5 * mu2e_cr0->size().x()-0.1,
                            0.5 * mu2e_cr0->size().y()-0.1,
                            0.5 * mu2e_cr0->size().z()     );
// --------------------------------------------------------------------------
// Draw disks and crystals
// --------------------------------------------------------------------------
  for (unsigned idisk=0; idisk<fCaloPtr->nDisks(); ++idisk) {
    const mu2e::Disk& mu2e_disk = fCaloPtr->disk(idisk);
    
    const auto& info            = mu2e_disk.diskInfo();

    const double halfThickness  = 0.5 * info.size().z();

    auto diskShape = new TGeoTube(Form("DiskShape_%u", idisk),
                                  info.innerEnvelopeR(),
                                  info.outerEnvelopeR(),
                                  halfThickness
                                  );

    auto disk = new TEvdDisk(Form("Disk_%u", idisk),diskShape,diskMedium);
    
    disk->GetVolume()->SetLineColor(kAzure + 2);
    disk->GetVolume()->SetFillColor(kAzure + 2);
    disk->GetVolume()->SetTransparency(50);
//----------------------------------------------------------------------------------------
//  * DiskInfo::toGlobal(local) is:
//
//         global = origin + inverseRotation * local
//
// Therefore inverseRotation is the local-to-global rotation.
//-----------------------------------------------------------------------------
    const auto& rotation = info.inverseRotation();

    Double_t matrix[9] = {
      rotation.xx(), rotation.xy(), rotation.xz(),
      rotation.yx(), rotation.yy(), rotation.yz(),
      rotation.zx(), rotation.zy(), rotation.zz()
    };

    auto geoRotation = new TGeoRotation;
    geoRotation->SetName(Form("DiskRotation_%u", idisk));
    geoRotation->SetMatrix(matrix);

    auto diskTransform = new TGeoCombiTrans(info.origin().x(),
                                            info.origin().y(),
                                            info.origin().z(),
                                            geoRotation);

    // Disk copy number is idisk + 1.
    gm->GetTopNode()->GetVolume()->AddNode(disk, idisk + 1,diskTransform);

    int nPlaced = 0;
    int ncr     = mu2e_disk.nCrystals();
    
    for (int i=0; i<ncr; ++i) {
      
      const mu2e::Crystal* mu2e_cr_i = &mu2e_disk.crystal(i);

      int index  = ncr*idisk+i;   // assume both diss to have the same number of crystals
      // every crystal has unique name - tied to readout channels
      TEvdCrystal* crystal = new TEvdCrystal(Form("crystal_%04i",index),shape,crystalMedium);

      // crystal color will depend on whether the crystal has hits
      crystal->SetLineColor(kOrange + 1);
      crystal->SetFillColor(kOrange + 1);
      crystal->SetTransparency(0);
      /*
       * Crystal::localPosition() is the front-face position - how do we know that?
       * TGeoBBox is centered on its local origin, so shift by +z/2.
       */
      
      const double x = mu2e_cr_i->localPosition().x();
      const double y = mu2e_cr_i->localPosition().y();
      const double z = mu2e_cr_i->localPosition().z() + 0.5 * mu2e_cr_i->size().z();

      std::cout << std::format("crystal i:{:4d} x:{:10.3f} y:{:10.3f} z:{:10.3f}\n",i,x,y,z);
//-----------------------------------------------------------------------------
// Copy number is local crystal ID + 1.  This is decoded by the right-click callback.
// crystal inherits from TGeoVolume
//-----------------------------------------------------------------------------
      disk->GetVolume()->AddNode(crystal,i+1,new TGeoTranslation(x,y,z));
      disk->AddCrystal(crystal);
      ++nPlaced;
    }
    std::cout << std::format("idisk:{} nPlaced:{:4d}\n",idisk,nPlaced);
  }
  return rc;
}

//-----------------------------------------------------------------------------
int TEvdCalorimeter::InitEvent() {
  std::cout << std::format("%s emoe AAAAAAA\n",__func__);
  return 0;
}

//-----------------------------------------------------------------------------
void TEvdCalorimeter::Draw(Option_t* Opt) {
  std::cout << std::format("DRAW :::: emoe AAAAAAA\n");
}

//-----------------------------------------------------------------------------
void TEvdCalorimeter::Print(Option_t* Opt) const {
  std::cout << std::format("emoe AAAAAAA\n");
}

}
