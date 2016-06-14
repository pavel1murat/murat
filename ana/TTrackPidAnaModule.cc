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
#include "Stntuple/alg/TStntuple.hh"
#include "Stntuple/obj/TDisk.hh"
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
  fMinT0 = 700; 

  fTrackID = new TStnTrackID();
  fLogLH   = new TEmuLogLH();
}

//-----------------------------------------------------------------------------
TTrackPidAnaModule::~TTrackPidAnaModule() {
}


//-----------------------------------------------------------------------------
void TTrackPidAnaModule::BookTrackHistograms(TrackHist_t* Hist, const char* Folder) {
//   char name [200];
//   char title[200];

  HBook1F(Hist->fP       ,"p"        ,Form("%s: Track P(Z1)"       ,Folder), 400,  90  ,110. ,Folder);
}

//-----------------------------------------------------------------------------
void TTrackPidAnaModule::BookEventHistograms(EventHist_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];

  HBook1F(Hist->fNTracks[0],"ntrk_0"     ,Form("%s: Number of Reconstructed Tracks"  ,Folder),5,0,5,Folder);
  HBook1F(Hist->fNTracks[1],"ntrk_1"     ,Form("%s: Number of Set C         Tracks"  ,Folder),5,0,5,Folder);

  HBook1F(Hist->fNDem,"ndem"   ,Form("%s: NTRK Dem"                        ,Folder),5,0,5,Folder);
  HBook1F(Hist->fNDmm,"ndmm"   ,Form("%s: NTRK Dmm"                        ,Folder),5,0,5,Folder);

  HBook2F(Hist->fNDmmVsNDem,"ndmm_vs_ndem"     ,Form("%s: N(Dmm) vs N(Dem)"          ,Folder),5,0,5,5,0,5,Folder);
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

  book_track_histset[  0] = 1;		// all tracks e-
  book_track_histset[  1] = 1;		// all tracks e- passing Set C cuts 
  
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
  int book_track_pid_histset[kNTrackPidHistSets];
  for (int i=0; i<kNTrackPidHistSets; i++) book_track_pid_histset[i] = 0;

  book_track_pid_histset[0] = 1;		// all tracks
  book_track_pid_histset[1] = 1;		// Set C tracks

  for (int i=0; i<kNTrackPidHistSets; i++) {
    if (book_track_pid_histset[i] != 0) {
      sprintf(folder_name,"pid_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fTrackPid[i] = new TrackPidHist_t;
      BookTrackPidHistograms(fHist.fTrackPid[i],Form("Hist/%s",folder_name));
    }
  }

}

//-----------------------------------------------------------------------------
// need MC truth branch
//-----------------------------------------------------------------------------
void TTrackPidAnaModule::FillEventHistograms(EventHist_t* Hist) {
  // double            p;
  // TLorentzVector    mom;

  Hist->fNTracks[0]->Fill  (fNTracks[0]);
  Hist->fNTracks[1]->Fill  (fNTracks[1]);

  Hist->fNDem->Fill(fNTracksDem);
  Hist->fNDmm->Fill(fNTracksDmm);

  Hist->fNDmmVsNDem->Fill(fNTracksDem,fNTracksDmm);
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
void TTrackPidAnaModule::FillTrackPidHistograms(TrackPidHist_t* Hist, TStnPid* Pid) {

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
void TTrackPidAnaModule::BookTrackPidHistograms(TrackPidHist_t* Hist, const char* Folder) {
//   char name [200];
//   char title[200];

  HBook1F(Hist->fLHEDedx     ,"lhe_dedx",Form("%s: LHEDedx"        ,Folder), 150,  -300. ,0. ,Folder);
  HBook1F(Hist->fLHMDedx     ,"lhm_dedx",Form("%s: LHMDedx"        ,Folder), 150,  -300. ,0. ,Folder);
  HBook1F(Hist->fLHRDedx     ,"lhr_dedx",Form("%s: LHRDedx"        ,Folder), 100,  -10.  , 10. ,Folder);
  HBook1F(Hist->fDrdsVadim   ,"drds",Form("%s: dr/ds vadim",Folder), 200,  -0.001 ,0.001 ,Folder);
  HBook1F(Hist->fDxdsVadim   ,"dxds",Form("%s: dx/ds vadim",Folder), 200,  -10   ,10   ,Folder);
  HBook1F(Hist->fSumAvik     ,"avik",Form("%s: sum Avik"   ,Folder), 200,    0   ,10 ,Folder);
  HBook1F(Hist->fMeanAvik    ,"mean",Form("%s: Avik sum/N"       ,Folder), 200,    0   ,0.2 ,Folder);

  HBook1F(Hist->fSq2Avik     ,"sq2",Form("%s: sq2"         ,Folder), 200,    0   ,10. ,Folder);
  HBook1F(Hist->fMq2Avik     ,"mq2",Form("%s: Avik sq2/Nmatchedall",Folder), 200,    0   ,0.2 ,Folder);

  HBook1F(Hist->fDrdsOs,"drds_os",Form("%s: dr/ds OS",Folder), 200,  -0.002 ,0.002 ,Folder);
  HBook1F(Hist->fDxdsOs,"dxds_os",Form("%s: dx/ds OS",Folder), 200,  -100. , 100. ,Folder);
  
  HBook1F(Hist->fDrdsSs,"drds_ss",Form("%s: dr/ds SS",Folder), 200,  -0.002 ,0.002 ,Folder);
  HBook1F(Hist->fDxdsSs,"dxds_ss",Form("%s: dx/ds SS",Folder), 200,  -100. ,100. ,Folder);

  HBook1F(Hist->fNUsedSsH ,"nu_ss_h" ,Form("%s: N(used) SS H",Folder), 100, 0,100 ,Folder);
  HBook1F(Hist->fNUsedOsH ,"nu_os_h" ,Form("%s: N(used) OS H",Folder), 100, 0,100 ,Folder);
  HBook1F(Hist->fNUsedOsD ,"nu_os_d" ,Form("%s: N(used) OS D",Folder), 100, 0,100 ,Folder);

  HBook1F(Hist->fSumAvikOs,"avik_os",Form("%s: Avik OS",Folder), 200, 0,100. ,Folder);

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
  RegisterDataBlock("TrackBlock"   ,"TStnTrackBlock",&fTrackBlockDem);
  RegisterDataBlock("TrackBlockDmm","TStnTrackBlock",&fTrackBlockDmm);
  RegisterDataBlock("PidBlock"     ,"TStnPidBlock"  ,&fTrackPidBlock);
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();
//-----------------------------------------------------------------------------
// initialize likelihood histograms
//-----------------------------------------------------------------------------
  fTrackID->SetMinT0(fMinT0);
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

  //  double       cos_th (-2.),  cl_e(-1.);
  //  int          disk_id(-1), alg_mask, nsh, nactive;
  //  float        pfront, ce_pitch, reco_pitch, fcons, t0, sigt, sigp, p; 
  TStnPid      *pid;
//-----------------------------------------------------------------------------
// event histograms
//-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0]);
//-----------------------------------------------------------------------------
// track histograms, fill them only for the downstream e- hypothesis
//-----------------------------------------------------------------------------
  TStnTrack*      trk;
  //  TrackPidPar_t*  tp;

  for (int i=0; i<fNTracks[0]; ++i ) {
    trk = fTrackBlockDem->Track(i);
    //    tp  = fTrackPidPar+i;

    FillTrackHistograms(fHist.fTrack[0],trk);

    if (trk->fIDWord == 0) {
					// track passes selection "C" 

      FillTrackHistograms(fHist.fTrack[1],trk);
    }
  }
//-----------------------------------------------------------------------------
// particle ID histograms
//-----------------------------------------------------------------------------
  for (int i=0; i<fNPid; ++i ) {
    trk = fTrackBlockDem->Track(i);
    pid = fTrackPidBlock->Pid(i);
    FillTrackPidHistograms(fHist.fTrackPid[0],pid);

    if ((trk->fIDWord == 0) && (trk->fP > 100.)) {

					// track with P > 103 MeV/c passes selection "C" 

      FillTrackPidHistograms(fHist.fTrackPid[1],pid);

      if (GetDebugBit(4) != 0) {
	GetHeaderBlock()->Print(Form("-----bit004: p = %12.5f T0: %10.3f",trk->fP,trk->fT0));
	printf("NMatched, NMatchedAll, sum, sq2  %3i %3i %12.5f %12.5f \n", 
	       pid->fNMatched,
	       pid->fNMatchedAll,
	       pid->fSumAvik,
	       pid->fSq2Avik);
      }
    }
  }
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

  fTrackBlockDem->GetEntry(ientry);
  fTrackBlockDmm->GetEntry(ientry);
  fTrackPidBlock->GetEntry(ientry);
//-----------------------------------------------------------------------------
// assume electron in the first particle, otherwise the logic will need to 
// be changed
//-----------------------------------------------------------------------------
  fNTracks[0] = fTrackBlockDem->NTracks();
  fNTracksDem = fTrackBlockDem->NTracks();
  fNTracksDmm = fTrackBlockDmm->NTracks();
  fNPid       = fTrackPidBlock->NTracks();

  fNTracks[1] = 0;

  if (fNTracks[0] == 0) fTrack = 0;
  else                  fTrack = fTrackBlockDem->Track(0);

  int ntrk = fNTracksDem;
  //  int alg_mask;

  //  TrackPidPar_t*   tp;

  for (int itrk=0; itrk<ntrk; itrk++) {
					// assume less 20 tracks
    //    tp             = fTrackPidPar+itrk;

    track          = fTrackBlockDem->Track(itrk);
    id_word        = fTrackID->IDWord(track);
    track->fIDWord = id_word;
    if (id_word == 0) {
      fNTracks[1] += 1;
    }
  }

  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TTrackPidAnaModule::Debug() {

  TStnTrack*     trk;
  //  TrackPidPar_t* tp;
  int ntrk = fTrackBlockDem->NTracks();

  for (int itrk=0; itrk<ntrk; itrk++) {
    trk = fTrackBlockDem->Track(itrk);
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

