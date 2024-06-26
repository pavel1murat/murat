//////////////////////////////////////////////////////////////////////////////
// use of TStnTrack::fTmp:
//
// Tmp(0) : corrected momentum at the tracker front - not yet
// 
// use of debug bits: 
//  0  : one line per track
//  1  : passed events
//  2  : rejected events
// 
//  5  : events with N(tracks) > 1
//
// 3 different ID : 
///////////////////////////////////////////////////////////////////////////////
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TEnv.h"
#include "TSystem.h"
#include "TDirectory.h"

#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/alg/TStntuple.hh"
#include "Stntuple/obj/TStnNode.hh"
#include "Stntuple/val/stntuple_val_functions.hh"

#include "murat/ana/TFilterModule.hh"

ClassImp(murat::TFilterModule)

namespace murat {
//-----------------------------------------------------------------------------
TFilterModule::TFilterModule(const char* name, const char* title): 
  TAnaModule(name,title) { 
  fTrackBlockName  = "TrackBlock";
  fMinNTracks      = 1;
//-----------------------------------------------------------------------------
// initialize Track ID
// 0: SetC  1: TrkQual>0.1 2:TrkQual>0.4
// what about number of hits ? - 3: no cuts on the number of hits
//-----------------------------------------------------------------------------
  fNID             = 1;

  fTrackID[0]      = new TStnTrackID();

  int mask0 = TStnTrackID::kNActiveBit | TStnTrackID::kChi2DofBit | TStnTrackID::kMomErrBit;
  fTrackID[0]->SetMaxChi2Dof(3.0 );
  fTrackID[0]->SetMinNActive(25  );
  fTrackID[0]->SetMaxMomErr (0.25);

  mask0    |= TStnTrackID::kTanDipBit;
  fTrackID[0]->SetMinTanDip (0.6);
  fTrackID[0]->SetMaxTanDip (1.0);

  mask0    |= TStnTrackID::kD0Bit;
  fTrackID[0]->SetMinD0 (-100);         // mm
  fTrackID[0]->SetMaxD0 ( 100);

  mask0    |= TStnTrackID::kT0Bit;
  fTrackID[0]->SetMinT0 ( 250);         // ns
  fTrackID[0]->SetMaxT0 (1700);

  fTrackID[0]->SetUseMask(mask0);

  for (int i=0; i<20; i++) {
    TrackPar_t* tp  = fTrackPar+i;
    tp->fTrackID[0] = TAnaModule::fTrackID[0];
  }
}

//-----------------------------------------------------------------------------
TFilterModule::~TFilterModule() {
  for (int i=0; i<fNID; i++) delete fTrackID[i];
}


//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TFilterModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
  RegisterDataBlock(fTrackBlockName.Data()        ,"TStnTrackBlock"     ,&fTrackBlock  );
  // RegisterDataBlock("HelixBlock"      ,"TStnHelixBlock"      ,&fHelixBlock      );
  // RegisterDataBlock("TimeClusterBlock","TStnTimeClusterBlock",&fTimeClusterBlock);

  RegisterDataBlock("ClusterBlock" ,"TStnClusterBlock" ,&fClusterBlock );
  // RegisterDataBlock("CalDataBlock" ,"TCalDataBlock"    ,&fCalDataBlock );
  // RegisterDataBlock("GenpBlock"    ,"TGenpBlock"       ,&fGenpBlock    );
  // RegisterDataBlock("SimpBlock"    ,"TSimpBlock"       ,&fSimpBlock    );
  // RegisterDataBlock("SpmcBlockVDet","TStepPointMCBlock",&fVDetBlock    );
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();

  return 0;
}

//-----------------------------------------------------------------------------
void TFilterModule::BookEventHistograms(HistBase_t* Hist, const char* Folder) {
  EventHist_t* hist = (EventHist_t*) Hist;
  HBook1F(hist->fNTracks[0] ,"ntrk_0"     ,Form("%s: Number of Reconstructed Tracks[0]"  ,Folder),100,0,100,Folder);
}

//_____________________________________________________________________________
void TFilterModule::BookHistograms() {

  //  char name [200];
  //  char title[200];

  TFolder* fol;
  TFolder* hist_folder;
  char     folder_name[200];

  DeleteHistograms();
  hist_folder = (TFolder*) GetFolder()->FindObject("Hist");
//-----------------------------------------------------------------------------
// book event histograms
//-----------------------------------------------------------------------------
  int book_event_histset[kNEventHistSets];
  for (int i=0; i<kNEventHistSets; i++) book_event_histset[i] = 0;

  book_event_histset[ 0] = 1;		// all events
  book_event_histset[ 1] = 1;	        // events with a reconstructed track
 
  for (int i=0; i<kNEventHistSets; i++) {
    if (book_event_histset[i] != 0) {
      sprintf(folder_name,"evt_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fEvent[i] = new EventHist_t;
      BookEventHistograms(fHist.fEvent[i],Form("Hist/%s",folder_name));
    }
  }
}

//_____________________________________________________________________________
void TFilterModule::FillEventHistograms(HistBase_t* Hist, EventPar_t* Evp) {

  EventHist_t* hist = (EventHist_t*) Hist;

  hist->fNTracks[0]->Fill(Evp->fNTracks[0]);
}


//_____________________________________________________________________________
void TFilterModule::FillHistograms() {

//-----------------------------------------------------------------------------
// event histograms
//-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0],&fEvtPar);
  if (fEvtPar.fNTracks[0]> 0) FillEventHistograms(fHist.fEvent[1],&fEvtPar);
}


//-----------------------------------------------------------------------------
// 2014-04-30: it looks that reading the straw hits takes a lot of time - 
//              turn off by default by commenting it out
//-----------------------------------------------------------------------------
int TFilterModule::Event(int ientry) {

  //  TDiskCalorimeter::GeomData_t disk_geom;

  fTrackBlock->GetEntry(ientry);
  fClusterBlock->GetEntry(ientry);
  // fTrackSeedBlock->GetEntry(ientry);
  // fHelixBlock->GetEntry(ientry);
  // fTimeClusterBlock->GetEntry(ientry);
  //  fTrackStrawHitBlock->GetEntry(ientry);
  // fStrawHitBlock->GetEntry(ientry);
  // fCalDataBlock->GetEntry(ientry);
  // fGenpBlock->GetEntry(ientry);
  // fSimpBlock->GetEntry(ientry);
  // fVDetBlock->GetEntry(ientry);
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// may want to revisit the definition of fSimp in the future
//-----------------------------------------------------------------------------
  fSimPar.fParticle = NULL;
  fSimPar.fTFront   = NULL;
  fSimPar.fTMid     = NULL;
  fSimPar.fTBack    = NULL;
//-----------------------------------------------------------------------------
// process virtual detectors - for fSimp need parameters at tracker entrance
// use the first fit
//-----------------------------------------------------------------------------
  fSimPar.fTFront   = NULL;
  fSimPar.fTMid     = NULL;
  fSimPar.fTBack    = NULL;

  int ntrk = fTrackBlock->NTracks();

  InitTrackPar(fTrackBlock,fClusterBlock,fTrackPar,&fSimPar);

  fEvtPar.fNTracks[0] = ntrk;

  int n_good_trk = 0;
  for (int itrk=0; itrk<ntrk; itrk++) {
    // TStnTrack*   trk = fTrackBlock->Track(itrk);
    TrackPar_t*  tp  = fTrackPar+itrk;

    int id_word = tp->fIDWord[0];
    if (id_word == 0) n_good_trk++;
  }

  fEventPassedSelections = (n_good_trk > 0);

  SetPassed(fEventPassedSelections);

  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TFilterModule::Debug() {

  // TStnTrack* trk;
  // TrackPar_t* tp;
  // int ntrk = fTrackBlock->NTracks();

//   for (int itrk=0; itrk<ntrk; itrk++) {
//     trk = fTrackBlock->Track(itrk);
//     tp  = &fTrackPar[itrk];
// //-----------------------------------------------------------------------------
// // bit 0: all tracks
// //-----------------------------------------------------------------------------
//     if (GetDebugBit(0) == 1) {
//       GetHeaderBlock()->Print(Form("bit_000: All p = %10.3f",trk->Momentum()->P()));
//     }
// //-----------------------------------------------------------------------------
// // bit 1: IDWord =0 0 tracks
// //-----------------------------------------------------------------------------
//     if (GetDebugBit(1) == 1) {
//       if (tp->fIDWord[fBestID] == 0) {
// 	GetHeaderBlock()->Print(Form("bit_001: IDWord=0 p = %10.3f",trk->Momentum()->P()));
//       }
//     }
//   }
//-----------------------------------------------------------------------------
// bit 5: events with N(tracks) > 1
//-----------------------------------------------------------------------------
  if (GetDebugBit(5) == 1) {
    int ntrk = fTrackBlock->NTracks();
    if (ntrk >= fMinNTracks) {
      GetHeaderBlock()->Print(Form("NTracks = %5i",ntrk));
    }
  }
}

//_____________________________________________________________________________
int TFilterModule::EndJob() {
  printf("----- end job: ---- %s\n",GetName());
  return 0;
}

}
