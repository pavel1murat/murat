#


#include "murat/gui/TEvdManager.hh"
#include "murat/gui/TEvdCaloVisNode.hh"
#include "murat/gui/TEvdCalorimeter.hh"


// ClassImp(TEvdCaloVisNode)

namespace murat {

//-----------------------------------------------------------------------------
  TEvdCaloVisNode::TEvdCaloVisNode(const char* Name) {
    
    auto vm = TEvdManager::Instance();  // has to be initialized at this point
    
    TEvdCalorimeter* calo = (TEvdCalorimeter*) vm->FindSubdetector("CALO");

    // prepare flat view of the two disks
    fAssembly = new TGeoVolumeAssembly("calo_disks");
    fAssembly->AddNode(calo->fDisk[0],1, new TGeoTranslation(-800, 0, 0));
    fAssembly->AddNode(calo->fDisk[1],1, new TGeoTranslation( 800, 0, 0));
  }

//-----------------------------------------------------------------------------
  TEvdCaloVisNode::~TEvdCaloVisNode() {
  }

//-----------------------------------------------------------------------------
  int TEvdCaloVisNode::InitEvent() {
    int rc(0);

    auto vm = TEvdManager::Instance();  // has to be initialized at this point
    
    TEvdCalorimeter* calo = (TEvdCalorimeter*) vm->FindSubdetector("CALO");

    // reset colors
    for (int i=0; i<2; i++) {
      TEvdDisk* disk = calo->Disk(i);
      int ncr = disk->NCrystals();
      for (int icr=0; icr<ncr; icr++) {
        TEvdCrystal* cr = disk->Crystal(icr);
        cr->ListOfHits()->Clear();
        cr->fEDep = 0;
        
        // crystal color will depend on whether the crystal has hits
        cr->SetLineColor(kOrange + 1);
        cr->SetFillColor(kOrange + 1);
        cr->SetTransparency(0);

      }
    }
    
    int nhits = fCaloHitBlock->NHits();
    for (int i=0; i<nhits; i++) {
      TCaloHit* hit = fCaloHitBlock->Hit(i);
      int idisk = hit->Disk();
      int icr   = hit->Cid() % 674; // 674 crystals per disk
      TEvdCrystal* cr = calo->Disk(idisk)->Crystal(icr);

      cr->AddHit(hit);
      
      cr->SetLineColor(kRed + 1);
      cr->SetFillColor(kRed + 1);
    }
    
    return rc;
  }

//-----------------------------------------------------------------------------
  void TEvdCaloVisNode::Draw(Option_t* Opt) {
    fAssembly->Draw(Opt);
  }
}
