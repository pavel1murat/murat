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

  double lhr_dedx, lhr_sum(-1.e10), lhr_sq2(-1.e10), lhr_osds(-1.e-10), lhr_osds_norm(-1.e-10) ;

  //  int              itrk;
  //  TrackPidPar_t*   tp;
					// pointer to local track parameters
  //  itrk = Track->Number();
  //  tp   = fTrackPidPar+itrk;

  Hist->fLHEDedx->Fill (Pid->fLogDedxProbEle);
  Hist->fLHMDedx->Fill (Pid->fLogDedxProbMuo);
  
  lhr_dedx = Pid->fLogDedxProbEle-Pid->fLogDedxProbMuo;

  Hist->fLHRDedx->Fill(lhr_dedx);

  Hist->fDrdsVadimEle->Fill(Pid->fDrdsVadimEle);
  Hist->fDrdsVadimMuo->Fill(Pid->fDrdsVadimMuo);

  float dxds_ele = Pid->fDrdsVadimEle/Pid->fDrdsVadimEleErr;
  float dxds_muo = Pid->fDrdsVadimMuo/Pid->fDrdsVadimMuoErr;

  Hist->fDxdsVadimEle->Fill(dxds_ele);
  Hist->fDxdsVadimMuo->Fill(dxds_muo);
					// these are sums over doublets
  Hist->fSumAvikEle->Fill(Pid->fSumAvikEle);
  Hist->fSumAvikMuo->Fill(Pid->fSumAvikMuo);
  Hist->fMeanAvikEle->Fill(Pid->fSumAvikEle/Pid->fNMatched);
  Hist->fMeanAvikMuo->Fill(Pid->fSumAvikMuo/Pid->fNMatched);

  if (Pid->fSumAvikEle > 0) lhr_sum = log(Pid->fSumAvikMuo/Pid->fSumAvikEle);

  Hist->fLHRSumAvik->Fill(lhr_sum);

  Hist->fSq2AvikEle->Fill(Pid->fSq2AvikEle);
  Hist->fSq2AvikMuo->Fill(Pid->fSq2AvikMuo);
  Hist->fMq2AvikEle->Fill(Pid->fSq2AvikEle/Pid->fNMatchedAll);
  Hist->fMq2AvikMuo->Fill(Pid->fSq2AvikMuo/Pid->fNMatchedAll);

  lhr_sq2 = -1.e10;
  if (Pid->fSq2AvikEle > 0) lhr_sq2 = log(Pid->fSq2AvikMuo/Pid->fSq2AvikEle);

  Hist->fLHRSq2Avik->Fill(lhr_sq2);

  Hist->fDrdsOsEle->Fill(Pid->fDrdsOsEle);
  Hist->fDrdsOsMuo->Fill(Pid->fDrdsOsMuo);

  double xdrds_os_ele = Pid->fDrdsOsEle/Pid->fDrdsOsEleErr;
  double xdrds_os_muo = Pid->fDrdsOsMuo/Pid->fDrdsOsMuoErr;

  Hist->fDxdsOsEle->Fill(xdrds_os_ele);
  Hist->fDxdsOsMuo->Fill(xdrds_os_muo);

  Hist->fDrdsSsEle->Fill(Pid->fDrdsSsEle);
  Hist->fDrdsSsMuo->Fill(Pid->fDrdsSsMuo);

  double xdrds_ss_ele = Pid->fDrdsSsEle/Pid->fDrdsSsEleErr;
  double xdrds_ss_muo = Pid->fDrdsSsMuo/Pid->fDrdsSsMuoErr;

  Hist->fDxdsSsEle->Fill(xdrds_ss_ele);
  Hist->fDxdsSsMuo->Fill(xdrds_ss_muo);

  Hist->fNUsedSsEleH->Fill(Pid->fNUsedSsEleH);
  Hist->fNUsedSsMuoH->Fill(Pid->fNUsedSsMuoH);
//-----------------------------------------------------------------------------
// osds: 'opposite side doublet slopes'
//-----------------------------------------------------------------------------
  Hist->fSumAvikOsEle->Fill(Pid->fSumAvikOsEle);
  Hist->fSumAvikOsMuo->Fill(Pid->fSumAvikOsMuo);

  Hist->fNUsedOsEleH->Fill(Pid->fNUsedOsEleH);
  Hist->fNUsedOsMuoH->Fill(Pid->fNUsedOsMuoH);

  Hist->fNUsedOsEleD->Fill(Pid->fNUsedOsEleD);
  Hist->fNUsedOsMuoD->Fill(Pid->fNUsedOsMuoD);

  if (Pid->fSumAvikOsEle > 0) {
    lhr_osds      = log(Pid->fSumAvikOsMuo/Pid->fSumAvikOsEle);
    lhr_osds_norm = log((Pid->fSumAvikOsMuo/Pid->fNUsedOsMuoD)/(Pid->fSumAvikOsEle/Pid->fNUsedOsEleD));
  }
  Hist->fLHROsDs->Fill(lhr_osds);
  Hist->fLHROsDsNorm->Fill(lhr_osds_norm);
//-----------------------------------------------------------------------------
// now ony need to combine the time and dedx likelihoods - later
//-----------------------------------------------------------------------------
  double lhr_sq2_dedx = lhr_dedx + lhr_sq2;
  Hist->fLHRSq2Dedx->Fill(lhr_sq2_dedx);

  double lhr_osdsn_dedx = lhr_dedx + lhr_osds_norm;
  Hist->fLHROsdsnDedx->Fill(lhr_osdsn_dedx);

  Hist->fLHRDedxVsSq2->Fill(lhr_sq2,lhr_dedx);

}


//-----------------------------------------------------------------------------
void TTrackPidAnaModule::BookTrackPidHistograms(TrackPidHist_t* Hist, const char* Folder) {
//   char name [200];
//   char title[200];

  HBook1F(Hist->fLHEDedx     ,"lhe_dedx",Form("%s: LHEDedx"        ,Folder), 150,  -300. ,0. ,Folder);
  HBook1F(Hist->fLHMDedx     ,"lhm_dedx",Form("%s: LHMDedx"        ,Folder), 150,  -300. ,0. ,Folder);
  HBook1F(Hist->fLHRDedx     ,"lhr_dedx",Form("%s: LHRDedx"        ,Folder), 100,  -10.  , 10. ,Folder);
  HBook1F(Hist->fDrdsVadimEle,"drds_ele",Form("%s: dr/ds vadim ele",Folder), 200,  -0.001 ,0.001 ,Folder);
  HBook1F(Hist->fDrdsVadimMuo,"drds_muo",Form("%s: dr/ds vadim muo",Folder), 200,  -0.001 ,0.001 ,Folder);
  HBook1F(Hist->fDxdsVadimEle,"dxds_ele",Form("%s: dx/ds vadim ele",Folder), 200,  -10   ,10   ,Folder);
  HBook1F(Hist->fDxdsVadimMuo,"dxds_muo",Form("%s: dx/ds vadim muo",Folder), 200,  -10   ,10.  ,Folder);
  HBook1F(Hist->fSumAvikEle  ,"avik_ele",Form("%s: Avik ele"       ,Folder), 200,    0   ,10 ,Folder);
  HBook1F(Hist->fSumAvikMuo  ,"avik_muo",Form("%s: Avik muo"       ,Folder), 200,    0   ,10 ,Folder);
  HBook1F(Hist->fMeanAvikEle  ,"mean_ele",Form("%s: Avik sum/N ele"       ,Folder), 200,    0   ,0.2 ,Folder);
  HBook1F(Hist->fMeanAvikMuo  ,"mean_muo",Form("%s: Avik sum/N muo"       ,Folder), 200,    0   ,0.2 ,Folder);
  HBook1F(Hist->fLHRSumAvik  ,"lhr_sum" ,Form("%s: lhr sum"        ,Folder), 200,  -10   ,10. ,Folder);

  HBook1F(Hist->fSq2AvikEle  ,"sq2_ele",Form("%s: sq2 ele"         ,Folder), 200,    0   ,10. ,Folder);
  HBook1F(Hist->fSq2AvikMuo  ,"sq2_muo",Form("%s: sq2 muo"         ,Folder), 200,    0   ,10. ,Folder);
  HBook1F(Hist->fMq2AvikEle  ,"mq2_ele",Form("%s: Avik sq2/Nmatchedall ELE",Folder), 200,    0   ,0.2 ,Folder);
  HBook1F(Hist->fMq2AvikMuo  ,"mq2_muo",Form("%s: Avik sq2/Nmatchedall MUO",Folder), 200,    0   ,0.2 ,Folder);
  HBook1F(Hist->fLHRSq2Avik  ,"lhr_sq2",Form("%s: lhr sq2"         ,Folder), 200,  -10   ,10. ,Folder);

  HBook1F(Hist->fDrdsOsEle,"drds_os_ele",Form("%s: dr/ds OS ele",Folder), 200,  -0.002 ,0.002 ,Folder);
  HBook1F(Hist->fDrdsOsMuo,"drds_os_muo",Form("%s: dr/ds OS muo",Folder), 200,  -0.002 ,0.002 ,Folder);
  HBook1F(Hist->fDxdsOsEle,"dxds_os_ele",Form("%s: dx/ds OS ele",Folder), 200,  -100. , 100. ,Folder);
  HBook1F(Hist->fDxdsOsMuo,"dxds_os_muo",Form("%s: dx/ds OS muo",Folder), 200,  -100. , 100. ,Folder);
  
  HBook1F(Hist->fDrdsSsEle,"drds_ss_ele",Form("%s: dr/ds SS ele",Folder), 200,  -0.002 ,0.002 ,Folder);
  HBook1F(Hist->fDrdsSsMuo,"drds_ss_muo",Form("%s: dr/ds SS muo",Folder), 200,  -0.002 ,0.002 ,Folder);
  HBook1F(Hist->fDxdsSsEle,"dxds_ss_ele",Form("%s: dx/ds SS ele",Folder), 200,  -100. ,100. ,Folder);
  HBook1F(Hist->fDxdsSsMuo,"dxds_ss_muo",Form("%s: dx/ds SS muo",Folder), 200,  -100. ,100. ,Folder);

  HBook1F(Hist->fNUsedSsEleH ,"nu_ss_eleh" ,Form("%s: N(used) SS EleH",Folder), 100, 0,100 ,Folder);
  HBook1F(Hist->fNUsedSsMuoH ,"nu_ss_muoh" ,Form("%s: N(used) SS MuoH",Folder), 100, 0,100 ,Folder);

  HBook1F(Hist->fSumAvikOsEle,"avik_os_ele",Form("%s: Avik OS Ele",Folder), 200, 0,100. ,Folder);
  HBook1F(Hist->fSumAvikOsMuo,"avik_os_muo",Form("%s: Avik OS Muo",Folder), 200, 0,100. ,Folder);

  HBook1F(Hist->fNUsedOsEleH ,"nu_os_eleh" ,Form("%s: N(used) OS EleH",Folder), 100, 0,100 ,Folder);
  HBook1F(Hist->fNUsedOsMuoH ,"nu_os_muoh" ,Form("%s: N(used) OS MuoH",Folder), 100, 0,100 ,Folder);

  HBook1F(Hist->fNUsedOsEleD ,"nu_os_eled" ,Form("%s: N(used) OS EleD",Folder), 100, 0,100 ,Folder);
  HBook1F(Hist->fNUsedOsMuoD ,"nu_os_muod" ,Form("%s: N(used) OS MuoD",Folder), 100, 0,100 ,Folder);

  HBook1F(Hist->fLHROsDs     ,"lhr_osds"   ,Form("%s: LHR OS doublets",Folder), 200, -10,10 ,Folder);
  HBook1F(Hist->fLHROsDsNorm ,"lhr_osds_n" ,Form("%s: LHR OS doub N"  ,Folder), 200, -10,10 ,Folder);

  HBook1F(Hist->fLHRSq2Dedx  ,"lhr_sq2_dedx",Form("%s: LHR DeDx+Sq2",Folder), 400, -20,20 ,Folder);
  HBook1F(Hist->fLHROsdsnDedx,"lhr_osdsn_dedx",Form("%s: LHR OsdsN+Dedx",Folder), 400, -20,20 ,Folder);

  HBook2F(Hist->fLHRDedxVsSq2,"lhr_dedx_vs_sq2",Form("%s: LHR OsdsN+Dedx",Folder),200,-10,10,200,-10,10,Folder);
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
	printf("NMatched, NMatchedAll, sum_ele, sum_muo, sq2_ele, sq2_muo: %3i %3i %12.5f %12.5f %12.5f %12.5f \n", 
	       pid->fNMatched,
	       pid->fNMatchedAll,
	       pid->fSumAvikEle,
	       pid->fSumAvikMuo,
	       pid->fSq2AvikEle,
	       pid->fSq2AvikMuo);
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

