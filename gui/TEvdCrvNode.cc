#

#include "murat/gui/TEvdManager.hh"
#include "murat/gui/TEvdCrvNode.hh"
#include "murat/gui/TEvdCrv.hh"


// ClassImp(TEvdCrvNode)

namespace murat {

//-----------------------------------------------------------------------------
  TEvdCrvNode::TEvdCrvNode(const char* Name) {
    
    auto vm = TEvdManager::Instance();  // has to be initialized at this point
    
    TEvdCrv* crv = (TEvdCrv*) vm->FindSubdetector("CRV");

    // prepare flat view of the two disks
    fAssembly = new TGeoVolumeAssembly("crv");
    int ns = crv->GetNSectors();
    
    for (int i=0; i<ns; i++) {
      // need to figure the transformation so all bars woudl be visible
      fAssembly->AddNode(crv->GetSector(i),1, new TGeoTranslation(-800, 0, 0));
    }
  }

//-----------------------------------------------------------------------------
  TEvdCrvNode::~TEvdCrvNode() {
  }

//-----------------------------------------------------------------------------
  int TEvdCrvNode::InitEvent() {
    int rc(0);

    auto vm = TEvdManager::Instance();  // has to be initialized at this point
    
    TEvdCrv* crv = (TEvdCrv*) vm->FindSubdetector("CRV");

    rc = crv->InitEvent();
        
    return rc;
  }

//-----------------------------------------------------------------------------
  bool TEvdCrvNode::Initialized() {
    bool initialized(false);

    // auto vm = TEvdManager::Instance();
    
    if ((fCrvpBlock->EventNumber () == f_EventNumber ) and
        (fCrvpBlock->SubrunNumber() == f_SubrunNumber) and
        (fCrvpBlock->RunNumber   () == f_RunNumber   )     ) {
      initialized = true;
    }
    return initialized;
  }

//-----------------------------------------------------------------------------
  void TEvdCrvNode::Draw(Option_t* Opt) {
    fAssembly->Draw(Opt);
  }
}
