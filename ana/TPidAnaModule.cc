///////////////////////////////////////////////////////////////////////////////
// use of tmp:
//
// Tmp(0) : 
// Tmp(1) : 
// ...
// 
// use of debug bits: bits 0-2 are reserved
// 0  : all events
// 1  : passed events
// 2  : rejected events
// 
// 3  : UNUSED
// 4  : UNUSED
// 5  : UNUSED
///////////////////////////////////////////////////////////////////////////////
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TEnv.h"
#include "TSystem.h"

#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/alg/TStntuple.hh"
#include "Stntuple/val/stntuple_val_functions.hh"
//------------------------------------------------------------------------------
#include "ana/TPidAnaModule.hh"

ClassImp(murat::TPidAnaModule)

namespace murat {
//-----------------------------------------------------------------------------
TPidAnaModule::TPidAnaModule(const char* name, const char* title):
  TAnaModule(name,title)
{
  fNID = 1;
  for (int i=0; i<fNID; i++) fTrackID[i] = new TStnTrackID();

  fTrackBlockName[0] = "TrackBlockPar";
  fTrackBlockName[1] = "TrackBlockDar";
}

//-----------------------------------------------------------------------------
TPidAnaModule::~TPidAnaModule() {
  for (int i=0; i<fNID; i++) delete fTrackID[i];
}


//-----------------------------------------------------------------------------
void TPidAnaModule::BookPidHistograms(PidHist_t* Hist, const char* Folder) {
  //     char name [200];
  //     char title[200];
  //-----------------------------------------------------------------------------
  //  
  //-----------------------------------------------------------------------------

  HBook1F(Hist->fTrkNumber ,"trk"        ,Form("%s: Trk Number" ,Folder), 100, 0, 100,Folder);
  HBook1F(Hist->fNMatched     ,"nmatched"    ,Form("%s: N(Matched)"     ,Folder), 100, 0, 100,Folder);
  HBook1F(Hist->fNMatchedAll  ,"nmatched_all",Form("%s: N(matched all"  ,Folder), 100, 0, 100,Folder);

  HBook1F(Hist->fNUsedOsH  ,"nused_os_h",Form("%s: N(used OS H)" ,Folder), 100, 0, 100,Folder);

  HBook1F(Hist->fNUsedSsH  ,"nused_ss_h",Form("%s: N(used SS H)" ,Folder), 100, 0, 100,Folder);

  HBook1F(Hist->fNUsedOsD  ,"nused_os_d",Form("%s: N(used OS D)" ,Folder), 100, 0, 100,Folder);

  HBook1F(Hist->fLogDedxProbEle,"log_dedx_prob_ele",Form("%s: Log DedxProb(Ele)" ,Folder), 200, -400,   0,Folder);
  HBook1F(Hist->fLogDedxProbMuo,"log_dedx_prob_muo",Form("%s: Log DedxProb(Muo)" ,Folder), 200, -400,   0,Folder);
  HBook1F(Hist->fLhrDedx,       "lhr_dedx"         ,Form("%s: LHR Dedx(ELE/Muo)" ,Folder), 200, -10 , 10 ,Folder);

  HBook1F(Hist->fDrdsVadim  ,"drds_vadim" ,Form("%s: DrDs Vadim" ,Folder), 200, -2.0e-3, 2.0e-3,Folder);
  HBook1F(Hist->fXdrdsVadim ,"xdrds_vadim",Form("%s: XDrDs Vadim" ,Folder), 400, -20, 20,Folder);

  HBook1F(Hist->fSumAvik  ,"sum_avik" ,Form("%s: Sum Avik" ,Folder), 500, 0, 5,Folder);

  HBook1F(Hist->fDrdsOs  ,"drds_os" ,Form("%s: DrDs  OS",Folder), 200, -2.0e-3, 2.0e-3,Folder);
  HBook1F(Hist->fXdrdsOs ,"xdrds_os",Form("%s: XdrDs OS",Folder), 400, -20, 20,Folder);

  HBook1F(Hist->fSumAvikOs  ,"sum_avik_os" ,Form("%s: Sum Avik OS" ,Folder), 500, 0, 10,Folder);
  //  HBook1F(Hist->fSumAvikOsMuo  ,"sum_avik_os_muo" ,Form("%s: Sum Avik OS Muo " ,Folder), 500, 0, 10,Folder);

  HBook1F(Hist->fDrdsSs  ,"drds_ss" ,Form("%s: DrDs SS ",Folder), 200, -2.0e-3, 2.0e-3,Folder);
  HBook1F(Hist->fXdrdsSs ,"xdrds_ss",Form("%s: XdrDs SS",Folder), 400, -20, 20,Folder);
//   HBook1F(Hist->fDrdsSsMuo  ,"drds_ss_muo" ,Form("%s: DrDs SS Muo ",Folder), 200, -2.0e-3, 2.0e-3,Folder);
//   HBook1F(Hist->fXdrdsSsMuo ,"xdrds_ss_muo",Form("%s: XdrDs SS Muo",Folder), 400, -20, 20,Folder);

}


//_____________________________________________________________________________
void TPidAnaModule::BookHistograms() {

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
  book_event_histset[ 1] = 0;		// 

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
// book PID histograms
//-----------------------------------------------------------------------------
  int book_pid_histset[kNPidHistSets];
  for (int i=0; i<kNPidHistSets; i++) book_pid_histset[i] = 0;

  book_pid_histset[0] = 1;		// all tracks
  book_pid_histset[1] = 1;		// "Set C" clusters

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
void TPidAnaModule::FillPidHistograms(PidHist_t* Hist, TStnPid* Pid) {

  Hist->fTrkNumber->Fill(Pid->fTrkNumber);

  Hist->fNMatched->Fill   (Pid->fNMatched);
  Hist->fNMatchedAll->Fill(Pid->fNMatchedAll);

  Hist->fNUsedOsH->Fill(Pid->fNUsedOsH);

  Hist->fNUsedOsD->Fill(Pid->fNUsedOsD);

  Hist->fNUsedSsH->Fill (Pid->fNUsedSsH);

  Hist->fLogDedxProbEle->Fill (Pid->fLogDedxProbEle);
  Hist->fLogDedxProbMuo->Fill (Pid->fLogDedxProbMuo);
  Hist->fLhrDedx->Fill        (Pid->fLogDedxProbEle-Pid->fLogDedxProbMuo);

  Hist->fDrdsVadim ->Fill(Pid->fDrdsVadim);
  Hist->fXdrdsVadim->Fill(Pid->fDrdsVadim/Pid->fDrdsVadimErr);

  Hist->fSumAvik->Fill(Pid->fSumAvik);

  Hist->fDrdsOs->Fill(Pid->fDrdsOs);
  Hist->fXdrdsOs->Fill(Pid->fDrdsOs/Pid->fDrdsOsErr);

  Hist->fSumAvikOs->Fill(Pid->fSumAvikOs);

  Hist->fDrdsSs->Fill(Pid->fDrdsSs);
  Hist->fXdrdsSs->Fill(Pid->fDrdsSs/Pid->fDrdsSsErr);

}

//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TPidAnaModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks for downstream electrons and muons
//-----------------------------------------------------------------------------
  RegisterDataBlock("PidBlockDem"     ,"TStnPidBlock"    , &fPidDataBlock[0]);
  RegisterDataBlock("PidBlockDmm"     ,"TStnPidBlock"    , &fPidDataBlock[1]);

  RegisterDataBlock(fTrackBlockName[0], "TStnTrackBlock" , &fTrackBlock  [0]);
  RegisterDataBlock(fTrackBlockName[1], "TStnTrackBlock" , &fTrackBlock  [1]);

  RegisterDataBlock("ClusterBlock"    ,"TStnClusterBlock", &fClusterBlock   );
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();

  return 0;
}


//_____________________________________________________________________________
void TPidAnaModule::FillHistograms() {

//-----------------------------------------------------------------------------
// event histograms
//
// EVT_0: all events
//-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0],&fEvtPar);
//-----------------------------------------------------------------------------
// 2020-10-04: to be implemented properly
//-----------------------------------------------------------------------------
//  int            nh;
  TStnPid*     pid;
  TrackPar_t*  tp;

  for (int i=0; i<fNTracks[0]; i++) {
    tp  = fTrackPar[0]+i;
    pid = fPidDataBlock[0]->Pid(i);
    if (pid) {
      FillPidHistograms(fHist.fPid[0],pid);
      if (tp->fIDWord == 0) {
	FillPidHistograms(fHist.fPid[1],pid);
      }
    }
    else {
      printf("TPidAnaModule::FillHistograms ERROR: pid[%i] = nullptr\n",i);
    }
  }
}



//_____________________________________________________________________________
int TPidAnaModule::BeginRun() {
  int rn = GetHeaderBlock()->RunNumber();
  TStntuple::Init(rn);
  return 0;
}


//_____________________________________________________________________________
int TPidAnaModule::Event(int ientry) {

  //  double                p;
  //  TLorentzVector        mom;

  fPidDataBlock[0]->GetEntry(ientry);
  fPidDataBlock[1]->GetEntry(ientry);

  fTrackBlock  [0]->GetEntry(ientry);
  fTrackBlock  [1]->GetEntry(ientry);

  fClusterBlock->GetEntry(ientry);
//-----------------------------------------------------------------------------
// event parameters - standardized
//-----------------------------------------------------------------------------
  fEvtPar.fNCrvClusters     = -1;
  fEvtPar.fNCrvPulses       = -1;
  fEvtPar.fNCrvCoincidences = -1;
  fEvtPar.fNGenp            = 0;  // fGenpBlock->NParticles();
  fEvtPar.fParticle         = nullptr;
//-----------------------------------------------------------------------------
// assume electron in the first particle, otherwise the logic will need to 
// be changed
// if there are several hits, use the first one
//-----------------------------------------------------------------------------
  fSimPar.fParticle = nullptr; // fSimpBlock->Particle(0);
  fSimPar.fTFront   = NULL;
  fSimPar.fTMid     = NULL;
  fSimPar.fTBack    = NULL;
  fSimPar.fGenp     = NULL;
//  fNGenp      = fGenpBlock->NParticles();
//  TStnPid* hit;

  fNTracks[0] = fTrackBlock[0]->NTracks();
  fNTracks[1] = fTrackBlock[1]->NTracks();

  for (int i=0; i<2; i++) {
    for (int itrk=0; itrk<fNTracks[i]; itrk++) {
      TrackPar_t* tp  = fTrackPar[i]+itrk;
      tp->fTrackID[0] = TAnaModule::fTrackID_BOX;
      tp->fTrackID[1] = TAnaModule::fTrackID_MVA;
      tp->fLogLH      = TAnaModule::fLogLH;
    }

    InitTrackPar(fTrackBlock[i],fClusterBlock,fTrackPar[i],&fSimPar);
  }

  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TPidAnaModule::Debug() {

//-----------------------------------------------------------------------------
// bit 4: events with NHitsTF > 1
//-----------------------------------------------------------------------------
  if (GetDebugBit(4) == 1) {
//     if (fNHitsTF > 1) {
//       GetHeaderBlock()->Print(Form("NHits(TF) = %5i",fNHitsTF));
//     }
  }
}

//_____________________________________________________________________________
int TPidAnaModule::EndJob() {
  printf("----- end job: ---- %s\n",GetName());
  return 0;
}

//_____________________________________________________________________________
void TPidAnaModule::Test001() {
}

}
