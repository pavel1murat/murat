//////////////////////////////////////////////////////////////////////////////
// use of tmp:
//
// Tmp(0) : nax seg
// Tmp(1) : nst seg
// 
// use of debug bits: bits 0-2 are reserved
//  0  : all events
//  1  : passed events
//  2  : rejected events
// 
///////////////////////////////////////////////////////////////////////////////
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TEnv.h"
#include "TSystem.h"

#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/obj/TTrackStrawHitBlock.hh"
#include "Stntuple/alg/TStntuple.hh"
#include "Stntuple/val/stntuple_val_functions.hh"
//------------------------------------------------------------------------------
// Mu2e offline includes
//-----------------------------------------------------------------------------
// #include "CalorimeterGeom/inc/HexMap.hh"

#include "ana/TTrackPidAnaModule.hh"

ClassImp(TTrackPidAnaModule)
//-----------------------------------------------------------------------------
TTrackPidAnaModule::TTrackPidAnaModule(const char* name, const char* title):
  TStnModule(name,title)
{
  for (int i=0; i<8; i++) {
    fTrackBlock   [i] = NULL;
    fTrackHitBlock[i] = NULL;
    fPidBlock     [i] = NULL;
  }
//-----------------------------------------------------------------------------
// track quality selection
//-----------------------------------------------------------------------------
  fNID             = 1;
  for (int i=0; i<fNID; i++) {
    fTrackID[i]    = new TStnTrackID();
  }

  fTrackID[0]->SetMaxMomErr (100);
  fTrackID[0]->SetMaxT0Err  (100);
  fTrackID[0]->SetMinFitCons(-1.);
  fTrackID[0]->SetMinTrkQual(0.4);
  fTrackID[0]->SetMinNActive(-1 );

  fBestID     = 0;			// best: DaveTrkQual > 0.4, no NActive cut
//-----------------------------------------------------------------------------
// track quality selection
//-----------------------------------------------------------------------------
  fLogLH   = new TEmuLogLH();
}

//-----------------------------------------------------------------------------
TTrackPidAnaModule::~TTrackPidAnaModule() {
}


//-----------------------------------------------------------------------------
void TTrackPidAnaModule::BookPidHistograms(PidHist_t* Hist, const char* Folder) {
//   char name [200];
//   char title[200];

  HBook1F(Hist->fLHEDedx   ,"lhe_dedx",Form("%s: LHEDedx"        ,Folder), 150,  -300. ,0. ,Folder);
  HBook1F(Hist->fLHMDedx   ,"lhm_dedx",Form("%s: LHMDedx"        ,Folder), 150,  -300. ,0. ,Folder);
  HBook1F(Hist->fLHRDedx   ,"lhr_dedx",Form("%s: LHRDedx"        ,Folder), 100,  -10.  , 10. ,Folder);
  HBook1F(Hist->fDrdsVadim ,"drds"    ,Form("%s: dr/ds vadim",Folder), 200,  -0.001 ,0.001 ,Folder);
  HBook1F(Hist->fDxdsVadim ,"dxds"    ,Form("%s: dx/ds vadim",Folder), 200,  -10   ,10   ,Folder);
  HBook1F(Hist->fSumAvik   ,"avik_sum",Form("%s: sum Avik"   ,Folder), 200,    0   ,0.2 ,Folder);
  HBook1F(Hist->fMeanAvik  ,"mean_sum",Form("%s: Avik sum/N"       ,Folder), 200,    0   ,0.2 ,Folder);

  HBook1F(Hist->fSq2Avik   ,"sq2"     ,Form("%s: sq2"         ,Folder), 200,    0   ,10. ,Folder);
  HBook1F(Hist->fMq2Avik   ,"mq2"     ,Form("%s: Avik sq2/Nmatchedall",Folder), 200,    0   ,10 ,Folder);

  HBook1F(Hist->fDrdsOs    ,"drds_os" ,Form("%s: dr/ds OS",Folder), 200,  -0.002 ,0.002 ,Folder);
  HBook1F(Hist->fDxdsOs    ,"dxds_os" ,Form("%s: dx/ds OS",Folder), 200,  -100. , 100. ,Folder);
  
  HBook1F(Hist->fDrdsSs    ,"drds_ss" ,Form("%s: dr/ds SS",Folder), 200,  -0.002 ,0.002 ,Folder);
  HBook1F(Hist->fDxdsSs    ,"dxds_ss" ,Form("%s: dx/ds SS",Folder), 200,  -100. ,100. ,Folder);

  HBook1F(Hist->fNUsedSsH  ,"nu_ss_h" ,Form("%s: N(used) SS H",Folder), 100, 0,100 ,Folder);
  HBook1F(Hist->fNUsedOsH  ,"nu_os_h" ,Form("%s: N(used) OS H",Folder), 100, 0,100 ,Folder);
  HBook1F(Hist->fNUsedOsD  ,"nu_os_d" ,Form("%s: N(used) OS D",Folder), 100, 0,100 ,Folder);

  HBook1F(Hist->fSumAvikOs ,"avik_os" ,Form("%s: Avik OS",Folder), 200, 0,100. ,Folder);

}

//-----------------------------------------------------------------------------
void TTrackPidAnaModule::BookEventHistograms(EventHist_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];


  for (int i=0; i<8; i++) {
    HBook1F(Hist->fNTracks[i],Form("ntrk_%i",i) ,Form("%s: Number of Reconstructed Tracks[%i]"  ,Folder,i),5,0,5,Folder);
  }
}

//-----------------------------------------------------------------------------
void TTrackPidAnaModule::BookTrackHistograms(TrackHist_t* Hist, const char* Folder) {
//   char name [200];
//   char title[200];

  HBook1F(Hist->fP       ,"p"        ,Form("%s: Track P(Z1)"       ,Folder), 400,  90  ,110. ,Folder);
}

//_____________________________________________________________________________
void TTrackPidAnaModule::BookHistograms() {

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

  for (int i=0; i<kNEventHistSets; i++) {
    if (book_event_histset[i] != 0) {
      sprintf(folder_name,"evt_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fEvent[i] = new EventHist_t;
      BookEventHistograms(fHist.fEvent[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book track histograms
//-----------------------------------------------------------------------------
  int book_track_histset[kNTrackHistSets];
  for (int i=0; i<kNTrackHistSets; i++) book_track_histset[i] = 0;

  for (int ih=0; ih<8; ih++) {
    int loc = 100*ih;
    book_track_histset[loc   ] = 1;		// all tracks 
    book_track_histset[loc+ 1] = 1;		// BEST_ID tracks
  }

  for (int i=0; i<kNTrackHistSets; i++) {
    if (book_track_histset[i] != 0) {
      sprintf(folder_name,"trk_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fTrack[i] = new TrackHist_t;
      BookTrackHistograms(fHist.fTrack[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book track PID histograms
//-----------------------------------------------------------------------------
  int book_pid_histset[kNPidHistSets];
  for (int i=0; i<kNPidHistSets; i++) book_pid_histset[i] = 0;

  for (int ih=0; ih<8; ih++) {
    int loc = 100*ih;
    book_pid_histset[loc   ] = 1;		// all tracks 
    book_pid_histset[loc+ 1] = 1;		// BEST_ID tracks
  }

  for (int i=0; i<kNPidHistSets; i++) {
    if (book_pid_histset[i] != 0) {
      sprintf(folder_name,"pid_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fPid[i] = new PidHist_t;
      BookPidHistograms(fHist.fPid[i],Form("Hist/%s",folder_name));
    }
  }

}

//-----------------------------------------------------------------------------
// need MC truth branch
//-----------------------------------------------------------------------------
void TTrackPidAnaModule::FillEventHistograms(EventHist_t* Hist) {
  // double            p;
  // TLorentzVector    mom;

  for (int i=0; i<8; i++) {
    Hist->fNTracks[i]->Fill(fNTracks[i]);
  }
}

//-----------------------------------------------------------------------------
// for DIO : ultimately, one would need to renormalize the distribution
//-----------------------------------------------------------------------------
void TTrackPidAnaModule::FillTrackHistograms(TrackHist_t* Hist, TStnTrack* Track) {

  //  int              itrk;
  //  TrackPidPar_t*   tp;
					// pointer to local track parameters
  //  itrk = Track->Number();
  //  tp   = fTrackPidPar+itrk;

  Hist->fP->Fill (Track->fP);
}

//-----------------------------------------------------------------------------
void TTrackPidAnaModule::FillPidHistograms(PidHist_t* Hist, TStnPid* Pid) {

  double lhr_dedx;

  //  int              itrk;
  //  TrackPidPar_t*   tp;
					// pointer to local track parameters
  //  itrk = Track->Number();
  //  tp   = fTrackPidPar+itrk;

  Hist->fLHEDedx->Fill (Pid->fLogDedxProbEle);
  Hist->fLHMDedx->Fill (Pid->fLogDedxProbMuo);
  
  lhr_dedx = Pid->fLogDedxProbEle-Pid->fLogDedxProbMuo;

  Hist->fLHRDedx->Fill(lhr_dedx);

  Hist->fDrdsVadim->Fill(Pid->fDrdsVadim);

  float dxds = Pid->fDrdsVadim/Pid->fDrdsVadimErr;

  Hist->fDxdsVadim->Fill(dxds);
					// these are sums over doublets
  Hist->fSumAvik->Fill(Pid->fSumAvik);
  Hist->fMeanAvik->Fill(Pid->fSumAvik/Pid->fNMatched);

  Hist->fSq2Avik->Fill(Pid->fSq2Avik);
  Hist->fMq2Avik->Fill(Pid->fSq2Avik/Pid->fNMatchedAll);

  Hist->fDrdsOs->Fill(Pid->fDrdsOs);

  double xdrds_os = Pid->fDrdsOs/Pid->fDrdsOsErr;

  Hist->fDxdsOs->Fill(xdrds_os);
  Hist->fDrdsSs->Fill(Pid->fDrdsSs);

  double xdrds_ss = Pid->fDrdsSs/Pid->fDrdsSsErr;

  Hist->fDxdsSs->Fill(xdrds_ss);
  Hist->fNUsedSsH->Fill(Pid->fNUsedSsH);
//-----------------------------------------------------------------------------
// osds: 'opposite side doublet slopes'
//-----------------------------------------------------------------------------
  Hist->fSumAvikOs->Fill(Pid->fSumAvikOs);
  Hist->fNUsedOsH->Fill(Pid->fNUsedOsH);

  Hist->fNUsedOsD->Fill(Pid->fNUsedOsD);
}


//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TTrackPidAnaModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks
// follow Stntuple/fcl/templates.fcl and assume that the first, DEM, 
// track data block is called just "TrackBlock"
//-----------------------------------------------------------------------------
  RegisterDataBlock("TrackBlockDem"   ,"TStnTrackBlock"     ,&fTrackBlock[kDem]);
  RegisterDataBlock("TrackBlockDmp"   ,"TStnTrackBlock"     ,&fTrackBlock[kDmp]);
  RegisterDataBlock("TrackBlockUmm"   ,"TStnTrackBlock"     ,&fTrackBlock[kUmm]);
  RegisterDataBlock("TrackBlockUmp"   ,"TStnTrackBlock"     ,&fTrackBlock[kUmp]);

  RegisterDataBlock("TrackHitBlockDem","TTrackStrawHitBlock",&fTrackHitBlock[kDem]);
  RegisterDataBlock("TrackHitBlockDmp","TTrackStrawHitBlock",&fTrackHitBlock[kDmp]);
  RegisterDataBlock("TrackHitBlockUmm","TTrackStrawHitBlock",&fTrackHitBlock[kUmm]);
  RegisterDataBlock("TrackHitBlockUmp","TTrackStrawHitBlock",&fTrackHitBlock[kUmp]);

  RegisterDataBlock("PidBlockDem"     ,"TStnPidBlock"       ,&fPidBlock[kDem]     );
  RegisterDataBlock("PidBlockDmp"     ,"TStnPidBlock"       ,&fPidBlock[kDmp]     );
  RegisterDataBlock("PidBlockUmm"     ,"TStnPidBlock"       ,&fPidBlock[kUmm]     );
  RegisterDataBlock("PidBlockUmp"     ,"TStnPidBlock"       ,&fPidBlock[kUmp]     );
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();
// //-----------------------------------------------------------------------------
// // initialize likelihood histograms
// //-----------------------------------------------------------------------------
//   fTrackID->SetMinT0(fMinT0);
//-----------------------------------------------------------------------------
// initialize likelihood histograms
// TRK 19: "Set C" plus reconstructed and matched cluster
//-----------------------------------------------------------------------------
  const char   *pid_version;
  pid_version = gEnv->GetValue("mu2e.PidVersion","_none_");
  fLogLH->Init(pid_version);

  return 0;
}


//_____________________________________________________________________________
int TTrackPidAnaModule::BeginRun() {
  int rn = GetHeaderBlock()->RunNumber();
  TStntuple::Init(rn);
  return 0;
}


//_____________________________________________________________________________
void TTrackPidAnaModule::FillHistograms() {

//-----------------------------------------------------------------------------
// event histograms
//-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0]);
//-----------------------------------------------------------------------------
// track histograms - just test
//-----------------------------------------------------------------------------
  TStnTrack*      trk;
  TStnPid*        pid;

  for (int ihyp=0; ihyp<8; ihyp++) {
    TStnTrackBlock* tb = fTrackBlock[ihyp];
    TStnPidBlock*   pb = fPidBlock  [ihyp];
    if (tb) {
      int loc = 100*ihyp;
      int nt  = tb->NTracks();
      for (int it=0; it<nt; it++ ) {
	trk = tb->Track(it);
	pid = pb->Pid(it);
//-----------------------------------------------------------------------------
// SET LOC+0: all tracks
//-----------------------------------------------------------------------------
	FillTrackHistograms(fHist.fTrack[loc],trk);
	FillPidHistograms(fHist.fPid    [loc],pid);

	if (trk->fIDWord == 0) {
//-----------------------------------------------------------------------------
// SET LOC+1: BEST_ID tracks
//-----------------------------------------------------------------------------
	  FillTrackHistograms(fHist.fTrack[loc+1],trk);
	  FillPidHistograms(fHist.fPid    [loc+1],pid);
	}
      }
    }
  }

}


//-----------------------------------------------------------------------------
int TTrackPidAnaModule::FindTrack(TAnaPart* Part, int Index) {
  TStnTrack           *t1, *t2;
  TTrackStrawHitData  *h1, *h2;
  int                  nh1, nh2, ncommon;

  Part->fTrack[Index] = NULL;

  if (fTrackBlock[Index] == NULL) return 0;

  t1  = Part->fTrack[0];
  nh1 = t1->NHits();

  int nt2 = fTrackBlock[Index]->NTracks();
  for (int it2=0; it2<nt2; it2++) {
    t2  = fTrackBlock[Index]->Track(it2);
    nh2 = t2->NHits();
					// check hit pattern
    ncommon = 0;
    for (int ih1=0; ih1<t1->NHits(); ih1++) {
      h1  = fTrackHitBlock[0]->Hit(0,ih1);
      for (int ih2=0; ih2<t2->NHits(); ih2++) {
	h2  = fTrackHitBlock[Index]->Hit(it2,ih2);

	if (h1->Index() == h2->Index()) {
					// common hit
	  ncommon += 1;
	}
      }
    }
//-----------------------------------------------------------------------------
// require more than 50% of hit overlap
//-----------------------------------------------------------------------------
    if (ncommon > (nh1+nh2)/4.) {
					// the same track
      Part->fTrack[Index] = t2;
      break;
    }
  }

  return 0;
}

//-----------------------------------------------------------------------------
// for each potential CE candidate (DEM track) make an analysis object and run PID 
//-----------------------------------------------------------------------------
int TTrackPidAnaModule::MakeParticles() {
  TStnTrack           *t1;
  int ndem = fTrackBlock[kDem]->NTracks();

  for (int it1=0; it1<ndem; it1++) {
    t1                      = fTrackBlock[kDem]->Track(it1);
    fPart[it1].fTrack[kDem] = t1;
//-----------------------------------------------------------------------------
// loop over all other blocks and initialize other allowed hypotheses
//-----------------------------------------------------------------------------
    for (int i=1; i<8; i++) {
      FindTrack(&fPart[it1],i);
    }
  }
//-----------------------------------------------------------------------------
// particle ID 
//-----------------------------------------------------------------------------
  for (int it1=0; it1<ndem; it1++) {
    t1                    = fTrackBlock[kDem]->Track(it1);
    if ((t1->Ep() > 0) && (t1->Ep() < 1.15)) {
//-----------------------------------------------------------------------------
// a track has an associated cluster - calculate track-calo likelihoods
// for an upstream hypotheses the expected energy distribution can't be predicted
// In this case, set log(probability) to a positive value and do not consider the 
// energy part of the likelihood when comparing the two hypotheses
//-----------------------------------------------------------------------------
      // #####
      
    }
    else {
//-----------------------------------------------------------------------------
// a track doesn't have a cluster - use tracker-only particle ID 
//-----------------------------------------------------------------------------
    }
    
  }

  return 0;
}

//-----------------------------------------------------------------------------
// 2014-04-30: it looks that reading the straw hits takes a lot of time - 
//              turn off by default by commenting it out
//-----------------------------------------------------------------------------
int TTrackPidAnaModule::Event(int ientry) {

  //  double                xs, p;
  //  TEmuLogLH::PidData_t  dat;
  TStnTrack*            track;
  int                   id_word;
  TLorentzVector        mom;

  //  TDiskCalorimeter::GeomData_t disk_geom;

  for (int i=0; i<8; i++) {
    if (fTrackBlock[i] != NULL) {
      fTrackBlock   [i]->GetEntry(ientry);
      fTrackHitBlock[i]->GetEntry(ientry);
      fPidBlock     [i]->GetEntry(ientry);

      fNTracks[i] = fTrackBlock[i]->NTracks();
    }
  }
//-----------------------------------------------------------------------------
// make analysis object
//-----------------------------------------------------------------------------
  MakeParticles();
//-----------------------------------------------------------------------------
// assume electron in the first particle, otherwise the logic will need to be changed
//-----------------------------------------------------------------------------
  int ndem = fTrackBlock[kDem]->NTracks();

  for (int i=0; i<ndem; i++) {
					// assume less 20 tracks
    //    tp             = fTrackPidPar+itrk;

    track          = fPart[i].fTrack[kDem];
    id_word        = fTrackID[fBestID]->IDWord(track);
    track->fIDWord = id_word;
  }

  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TTrackPidAnaModule::Debug() {

  TStnTrack*     trk;
  //  TrackPidPar_t* tp;
  int ntrk = fTrackBlock[kDem]->NTracks();

  for (int itrk=0; itrk<ntrk; itrk++) {
    trk = fTrackBlock[kDem]->Track(itrk);
    //    tp  = fTrackPidPar+itrk;
//-----------------------------------------------------------------------------
// bit 3: Set C tracks with large DX : 70mm < |DX| < 90mm
//-----------------------------------------------------------------------------
    if (GetDebugBit(3) == 1) {
      if (trk->fIDWord == 0) {
	// TStnTrack::InterData_t*    vr = trk->fVMaxEp; // residuals
	// if ((vr && (fabs(vr->fDx) > 70) && (fabs(vr->fDx) < 90))) {
	//   GetHeaderBlock()->Print(Form("large DX: %f",vr->fDx));
	// }
      }
    }
  }
}

//_____________________________________________________________________________
int TTrackPidAnaModule::EndJob() {
  printf("----- end job: ---- %s\n",GetName());
  return 0;
}

//_____________________________________________________________________________
void TTrackPidAnaModule::Test001() {
}

