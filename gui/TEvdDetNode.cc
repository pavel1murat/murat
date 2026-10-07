///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////


#include "murat/gui/TEvdManager.hh"
#include "murat/gui/TEvdDetNode.hh"

#include "murat/gui/TEvdCalorimeter.hh"
#include "murat/gui/TEvdTracker.hh"
#include "murat/gui/TEvdCrvSector.hh"
#include "murat/gui/TEvdCrvModule.hh"
#include "murat/gui/TEvdCrvLayer.hh"
#include "murat/gui/TEvdCrv.hh"

#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/obj/TCaloHitBlock.hh"

// ClassImp(TEvdDetNode)

namespace murat {

//-----------------------------------------------------------------------------
  TEvdDetNode::TEvdDetNode(const char* Name) : TEvdVisNode() {
    
    auto vm = TEvdManager::Instance();  // has to be initialized at this point
    
    // build the geometry - that can be done in standalone mode
    fAssembly = new TGeoVolumeAssembly("mu2e");

    auto calo = vm->FindSubdetector("CALO");
    fAssembly->AddNode(calo,1, new TGeoTranslation(0, 0, 0));

    auto crv = vm->FindSubdetector("CRV");
    fAssembly->AddNode(crv,1, new TGeoTranslation(0, 0, 0));

    auto tracker = vm->FindSubdetector("TRACKER");
    fAssembly->AddNode(tracker,1, new TGeoTranslation(0, 0, 0));
  }

//-----------------------------------------------------------------------------
  TEvdDetNode::~TEvdDetNode() {
  }

//-----------------------------------------------------------------------------
// here it is assumed that all data blocks have been set
//-----------------------------------------------------------------------------
  int TEvdDetNode::InitCalorimeter() {
    int rc(0);
    auto vm = TEvdManager::Instance();  // has to be initialized at this point
    
    TEvdCalorimeter* calo = (TEvdCalorimeter*) vm->FindSubdetector("CALO");

                                        // reset colors
    for (int i=0; i<2; i++) {
      TEvdDisk* disk = calo->Disk(i);
      int ncr = disk->NCrystals();
      for (int icr=0; icr<ncr; icr++) {
        TEvdCrystal* crystal = disk->Crystal(icr);
        crystal->Clear();
      }
    }

    TCaloHitBlock* chb = (TCaloHitBlock*) vm->GetDataBlock("CaloHitBlock");
    int nhits = chb->NHits();
    for (int i=0; i<nhits; i++) {
      TCaloHit* hit = chb->Hit(i);
      int idisk = hit->Disk();
      int icr   = hit->Cid() % 674; // 674 crystals per disk
      TEvdCrystal* crystal = calo->Disk(idisk)->Crystal(icr);
      crystal->AddHit(hit);
    }

    return rc;
  }
  
//-----------------------------------------------------------------------------
// here it is assumed that all data blocks have been set
//-----------------------------------------------------------------------------
  int TEvdDetNode::InitCrv() {
    int rc(0);
    auto vm = TEvdManager::Instance();  // has to be initialized at this point

    TEvdCrv* crv = (TEvdCrv*) vm->FindSubdetector("CRV");

    int nsectors = crv->GetNSectors();

    for (int i=0; i<nsectors; i++) {
      TEvdCrvSector* sector = crv->GetSector(i);
      sector->Clear();
    }

    TCrvPulseBlock* crvp = (TCrvPulseBlock*) vm->GetDataBlock("CrvpBlock");
    
    int npulses = crvp->NPulses();
    for (int i=0; i<npulses; i++) {
      TCrvRecoPulse* pulse = crvp->Pulse(i);
      int sbid = pulse->Sbid();

      TEvdCrvSector*  sector(nullptr);
      TEvdCrvModule*  module(nullptr);
      TEvdCrvLayer*   layer (nullptr);
      TEvdCrvCounter* counter(nullptr);

      int found = crv->GetCounterLocation(sbid,sector,module,layer,counter);

      if (found) {
        counter->AddRecoPulse(pulse);
      
        layer->IncrementNHits ();
        module->IncrementNHits();
        sector->IncrementNHits();
      }
      else {
        printf("ERROR: scintillation counter %i not found\n",sbid);
      }
    }

    return rc;
  }

//-----------------------------------------------------------------------------
// here it is assumed that all data blocks have been set
//-----------------------------------------------------------------------------
  int TEvdDetNode::InitTracker() {
    int rc(0);
    return rc;
  }

//-----------------------------------------------------------------------------
// here it is assumed that all data blocks have been set
//-----------------------------------------------------------------------------
  int TEvdDetNode::InitTracks() {
    int rc(0);
    
    auto vm = TEvdManager::Instance();  // has to be initialized at this point

    TStrTrackBlock* tb = (TStrTrackBlock*) vm->GetDataBlock("TrackBlock");
    int ntracks = tb->NTracks();
    for (int i=0; i<ntracks;i++) {
      
    }
    return rc;
  }

//-----------------------------------------------------------------------------
// here it is assumed that all data blocks have been set
//-----------------------------------------------------------------------------
  int TEvdDetNode::InitEvent() {
    int rc(0);

    auto vm = TEvdManager::Instance();

    // all subdetectors have to be initialized
    // list of tracks too - make list of tracks a separate node ?
    
    int nsd = vm->GetNSubdetectors();
    for (int i=0; i<nsd; i++) {
      TEvdSubdetector* sd = vm->GetSubdetector(i);
      if (sd->InitEvent() != 0) {
        rc = -1;
      }
    }

    // InitTracks();

    return rc;
  }

//-----------------------------------------------------------------------------
  bool TEvdDetNode::Initialized() {
    bool initialized(true);

    auto vm = TEvdManager::Instance();

    // all subdetectors have to be initialized
    // list of tracks too - make list of tracks a separate node ?
    
    int nsd = vm->GetNSubdetectors();
    for (int i=0; i<nsd; i++) {
      TEvdSubdetector* sd = vm->GetSubdetector(i);
      if (not sd->Initialized()) {
        initialized = false;
        break;
      }
    }
    
    return initialized;
  }

//-----------------------------------------------------------------------------
  void TEvdDetNode::Draw(Option_t* Opt) {
    fAssembly->Draw(Opt);
  }
}
