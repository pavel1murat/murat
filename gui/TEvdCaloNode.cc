#


#include "murat/gui/TEvdManager.hh"
#include "murat/gui/TEvdCaloNode.hh"
#include "murat/gui/TEvdCalorimeter.hh"


// ClassImp(TEvdCaloNode)

namespace murat {

//-----------------------------------------------------------------------------
  TEvdCaloNode::TEvdCaloNode(const char* Name) {
    
    auto vm = TEvdManager::Instance();  // has to be initialized at this point
    
    TEvdCalorimeter* calo = (TEvdCalorimeter*) vm->FindSubdetector("CALO");

    // prepare flat view of the two disks
    fAssembly = new TGeoVolumeAssembly("calo_disks");
    fAssembly->AddNode(calo->fDisk[0],1, new TGeoTranslation(-800, 0, 0));
    fAssembly->AddNode(calo->fDisk[1],1, new TGeoTranslation( 800, 0, 0));
  }

//-----------------------------------------------------------------------------
  TEvdCaloNode::~TEvdCaloNode() {
  }

//-----------------------------------------------------------------------------
  int TEvdCaloNode::InitEvent() {
    int rc(0);

    auto vm = TEvdManager::Instance();  // has to be initialized at this point
    
    TEvdCalorimeter* calo = (TEvdCalorimeter*) vm->FindSubdetector("CALO");

    rc = calo->InitEvent();
        
    return rc;
  }

//-----------------------------------------------------------------------------
  bool TEvdCaloNode::Initialized() {
    bool initialized(false);

    auto vm = TEvdManager::Instance();
    
    
    if ((fCaloHitBlock->EventNumber ()  == f_EventNumber) and
        (fCaloHitBlock->SubrunNumber() == f_SubrunNumber) and
        (fCaloHitBlock->RunNumber   () == f_RunNumber   )     ) {
      initialized = true;
    }
    return initialized;
  }

//-----------------------------------------------------------------------------
  void TEvdCaloNode::Draw(Option_t* Opt) {
    fAssembly->Draw(Opt);
  }
}
