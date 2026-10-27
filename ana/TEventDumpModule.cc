//////////////////////////////////////////////////////////////////////////////
// use of tmp:
//
// use of debug bits: bits 0-2 are reserved
//  0  : all events
//  1  : passed events
//  2  : rejected events
//  3  : N(Calc  disk 0 ) > 0 and N(calc disk1 > 0)
//  4  : N(Calc  disk 0 ) > 0 and N(calc disk1 > 0) and a track with N>=20 hits
//  5  : used for something
//  6  : fNCalh > 0
///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <format>

#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TEnv.h"
#include "TSystem.h"
#include "TGLViewer.h"

#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/alg/TStntuple.hh"
#include "Stntuple/geom/TDisk.hh"
#include "Stntuple/val/stntuple_val_functions.hh"
//------------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
#include "murat/gui/TEvdCaloView.hh"
#include "murat/gui/TEvdCaloNode.hh"
#include "murat/gui/TEvdDetNode.hh"
#include "murat/gui/TEvdManager.hh"

#include "murat/ana/TEventDumpModule.hh"


ClassImp(murat::TEventDumpModule)

namespace murat {
//-----------------------------------------------------------------------------
TEventDumpModule::TEventDumpModule(const char* name, const char* title):
  TStnModule(name,title)
{
  TH1::AddDirectory(0);
//-----------------------------------------------------------------------------
// detector numerology
//-----------------------------------------------------------------------------
  fTpm            = TrkPanelMap::Instance();
  fCaloChannelMap = TCaloChannelMap::Instance();
  
  fFrRef          = nullptr;
  fHist           = new Hist_t;
}

//-----------------------------------------------------------------------------
TEventDumpModule::~TEventDumpModule() {
}

//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TEventDumpModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
  RegisterDataBlock("CaloDigiBlock"    ,"TCaloDigiBlock"      ,&fCaloDigiBlock    );
  RegisterDataBlock("CaloHitBlock"     ,"TCaloHitBlock"       ,&fCaloHitBlock     );
  RegisterDataBlock("CaloRecoDigiBlock","TCaloRecoDigiBlock"  ,&fCaloRecoDigiBlock);
  RegisterDataBlock("CaloClusterBlock" ,"TStnClusterBlock"    ,&fCaloClusterBlock );
  
  RegisterDataBlock("CrvdBlock"        ,"TCrvDigiBlock"       ,&fCrvdBlock        );
  RegisterDataBlock("CrvpBlock"        ,"TCrvPulseBlock"      ,&fCrvpBlock        );
  RegisterDataBlock("CrvcBlock"        ,"TCrvClusterBlock"    ,&fCrvcBlock        );

  RegisterDataBlock("StrawDigiBlock"   ,"TStrawDigiBlock"     ,&fSdBlock          );
  RegisterDataBlock("StrawHitBlock"    ,"TStrawHitBlock"      ,&fShBlock          );
  RegisterDataBlock("ComboHitBlock"    ,"TComboHitBlock"      ,&fChBlock          );
  RegisterDataBlock("TimeClusterBlock" ,"TStnTimeClusterBlock",&fTcBlock          );
  RegisterDataBlock("TrackBlock"       ,"TStrTrackBlock"      ,&fTrackBlock       );
 
  return 0;
}


//_____________________________________________________________________________
int TEventDumpModule::BeginRun() {
  fRunNumber = GetHeaderBlock()->RunNumber();
  TStntuple::Init(fRunNumber);
  fTpm->Init(fRunNumber);
//-----------------------------------------------------------------------------
// book histograms at begin run, there will be a flag telling whether we need
// to re-book them for every new run
// don't do that by default
//-----------------------------------------------------------------------------
  // relative positions of subdetectors wrt each other
  int print_level(0);
  fDetPos.GetFromRunDb(fRunNumber,print_level);
  
  //  BookHistograms(fHist,fFolder);
  return 0;
}


//-----------------------------------------------------------------------------
// 2026-09-09
//-----------------------------------------------------------------------------
int TEventDumpModule::Event(int IEntry) {

  //  TLorentzVector        mom;

  fSdBlock->GetEntry(IEntry);
  fShBlock->GetEntry(IEntry);
  fChBlock->GetEntry(IEntry);
  fTcBlock->GetEntry(IEntry);
  fTrackBlock->GetEntry(IEntry);
  
  fCrvpBlock->GetEntry(IEntry);
  fCrvcBlock->GetEntry(IEntry);
  fCrvdBlock->GetEntry(IEntry);
  
  fCaloDigiBlock->GetEntry(IEntry);
  fCaloHitBlock->GetEntry(IEntry);
  fCaloRecoDigiBlock->GetEntry(IEntry);
  fCaloClusterBlock->GetEntry(IEntry);
//-----------------------------------------------------------------------------
// assume electron in the first particle, otherwise the logic will need to 
// be changed
//-----------------------------------------------------------------------------
  fNCh           = fChBlock->NHits();
  fNSh           = fShBlock->NHits();
  fNSd           = fSdBlock->NDigis();
  fNTc           = fTcBlock->NTimeClusters();
  fNTrk          = fTrackBlock->NTracks();

  fNCrvc         = fCrvcBlock->NClusters();
  fNCrvd         = fCrvdBlock->NDigis();
  fNCrvp         = fCrvpBlock->NPulses();
   
  fNCald         = fCaloDigiBlock->NDigis();
  fNCalrd        = fCaloRecoDigiBlock->NDigis();
  fNCalh         = fCaloHitBlock->NHits();
  fNCaloClusters = fCaloClusterBlock->NClusters();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TEventDumpModule::Debug() {

  if (GetDebugBit(3) == 1) {
    if ((fNCcDisk[0] > 0) and (fNCcDisk[1] > 0)) {
      GetHeaderBlock()->Print(Form("NCcDisk[0]:%2i fNCcDisk[1]:%2i",
                                   fNCcDisk[0],fNCcDisk[1]));
      fCaloClusterBlock->Print();
    }
  }

  if (GetDebugBit(4) == 1) {
    if ((fNCcDisk[0] > 0) and (fNCcDisk[1] > 0)) {
      if (fNTrk >= 1) {
        int nh = fTrackBlock->Track(0)->NHits();
        if (nh >= 20) {
          GetHeaderBlock()->Print(Form("NCcDisk[0]:%2i fNCcDisk[1]:%2i N_trk_hits = %3d",
                                       fNCcDisk[0],fNCcDisk[1],nh));
          fCaloClusterBlock->Print();
        }
      }
    }
  }

  if (GetDebugBit(5) == 1) {
    if (fMaxHitsPerCrystal > 1) {
      int n2plus = 0;
      int n_hit_crystals = fHitCrystals.size();
      for (int i=0; i<n_hit_crystals; i++) {
        crystal_t* cr = fHitCrystals.at(i);
        // figure the number of good hits in the crystal
        int nh = cr->NHits();
        int n_good_hits = 0;
        for (int ih=0; ih<nh; ih++) {
          TCaloHit* calh = cr->Hit(ih);
          if ((calh->fNSipms == 2) and (calh->fEDep > 10.)) {
            n_good_hits += 1;
          }
        }
        if (n_good_hits >= 2) {
          n2plus += 1;
          GetHeaderBlock()->Print(Form("CID:%4i n_good_hits:%2i fNtrk:%d",
                                       cr->Cid(),n_good_hits,fNTrk));
        }
      }
      if (n2plus > 0) {
          fCaloHitBlock->Print();
      }
    }
  }

  // something in the calorimeter
  if (GetDebugBit(6) == 1) {
    if (fNCalh > 0) {
      GetHeaderBlock()->Print(Form("bit_006: fNCalh:%d",fNCalh));
    }
  }

}

//_____________________________________________________________________________
int TEventDumpModule::EndJob() {
  return 0;
}


//-----------------------------------------------------------------------------
int TEventDumpModule::PrintTracks() {
  std::cout << std::format(" i      T0         Z0       Nx     Ny    Nz    Chi2D     Xc[0]      Yc[0]      Xc[1]      Yc[1]    DxCal[0]    DyCal[0]    DxCal[1]   DyCal[1]\n");
  for (int i=0; i<fNTrk; i++) {
    TStrTrack*   trk = fTrackBlock->Track(i);
    trk_param_t* tp  = &fListOfTrkParam.at(i);
    std::cout << std::format("{:2d} {:10.3f} {:10.3f} {:6.3f} {:6.3f} {:6.3f} {:7.2f}",
                             i,trk->fT0,trk->fZ0,trk->fNx,trk->fNy,trk->fNz,trk->fChi2/trk->fNDof);
    
    std::cout << std::format(" {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f}\n",
                             tp->x_disk [0],tp->y_disk [0],tp->x_disk[1],tp->y_disk[1],
                             tp->dx_calc[0],tp->dy_calc[0],tp->dx_calc[1],tp->dy_calc[1]);
  }
  return 0;
}
//-----------------------------------------------------------------------------
  int TEventDumpModule::InitEvd() {
    int rc(0);

    auto vm = TEvdManager::Instance();
    vm->InitGeometry(); // if not initialized
    
    vm->AddDataBlock("HeaderBlock"      ,GetHeaderBlock());
    vm->AddDataBlock("CaloHitBlock"     ,fCaloHitBlock);
    vm->AddDataBlock("CaloRecoDigiBlock",fCaloRecoDigiBlock);
    vm->AddDataBlock("CaloClusterBlock" ,fCaloClusterBlock);
    
    vm->AddDataBlock("CrvpBlock"        ,fCrvpBlock);
    vm->AddDataBlock("CrvcBlock"        ,fCrvcBlock);

    vm->AddDataBlock("TrackBlock"       ,fTrackBlock);
//-----------------------------------------------------------------------------
// calorimeter view
//-----------------------------------------------------------------------------
    TEvdCaloNode* node = new TEvdCaloNode("calo_node");
    TEvdCaloView*    view = new TEvdCaloView("calo_view",TGLViewer::kCameraOrthoXOY);
    view->AddNode(node);
    vm->AddView(view);
//-----------------------------------------------------------------------------
// 3D detector view
//-----------------------------------------------------------------------------
    TEvdDetNode* det_node = new TEvdDetNode("det_node");
    TEvdView*    det_view = new TEvdView   ("det_view",TGLViewer::kCameraPerspXOY);
    det_view->AddNode(det_node);
    vm->AddView(det_view);

    return rc;
  }
  
//-----------------------------------------------------------------------------
  int TEventDumpModule::DisplayEvent() {
    int rc(0);
    auto vm = TEvdManager::Instance();

    vm->DisplayEvent();

    return rc;
  }
  
}
