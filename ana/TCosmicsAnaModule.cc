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
//  3  : print NTracks[0-3] for all events
//  4  : mu+ passed all selections
//  5  : tracks passing all ID but TanDip with TanDip > 1
//  6  : print TanDip for all tracks
//
//  ReportedErrors: 
//  ---------------
//  0  : E/p discrepancy
///////////////////////////////////////////////////////////////////////////////
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TEnv.h"
#include "TSystem.h"

#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/alg/TStntuple.hh"
#include "Stntuple/geom/TDisk.hh"
#include "Stntuple/val/stntuple_val_functions.hh"
#include "Offline/DataProducts/inc/VirtualDetectorId.hh"
//------------------------------------------------------------------------------
// Mu2e offline includes
//-----------------------------------------------------------------------------
// #include "CalorimeterGeom/inc/HexMap.hh"

#include "ana/TCosmicsAnaModule.hh"


ClassImp(murat::TCosmicsAnaModule)

namespace murat {
//-----------------------------------------------------------------------------
TCosmicsAnaModule::TCosmicsAnaModule(const char* name, const char* title):
  TAnaModule(name,title),
  fTrackBlockName("TrackBlock")
{
  fPtMin  = 1.;
  fTrackNumber.Set(100);

  for (int i=0; i<kMaxNErrors; i++) {
    fError[i].fNReports   =  0;
    fError[i].fMaxNReports = 10;
  }

  fDiskCalorimeter = new TDiskCalorimeter();
  fCalorimeterType = 2;
					// cut on track quality only
  fNID        = 2;
  
  fBestID     = 0;
  fTrackID[0] = fTrackID_BOX;
  fTrackID[1] = fTrackID_MVA;
}

//-----------------------------------------------------------------------------
TCosmicsAnaModule::~TCosmicsAnaModule() {
}


//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TCosmicsAnaModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks
// note that the first one has name'TrackBlock'
//-----------------------------------------------------------------------------
  RegisterDataBlock(fTrackBlockName.Data(),"TStnTrackBlock"   ,&fTrackBlockDem  );
  RegisterDataBlock("TrackBlockDmm"       ,"TStnTrackBlock"   ,&fTrackBlockDmm  );

  RegisterDataBlock("TrackBlockDep" ,"TStnTrackBlock"   ,&fTrackBlockDep  );
  RegisterDataBlock("TrackBlockDmp" ,"TStnTrackBlock"   ,&fTrackBlockDmp  );

  RegisterDataBlock("TrackBlockUem" ,"TStnTrackBlock"   ,&fTrackBlockUem  );
  RegisterDataBlock("TrackBlockUmm" ,"TStnTrackBlock"   ,&fTrackBlockUmm  );

  RegisterDataBlock("TrackBlockUep" ,"TStnTrackBlock"   ,&fTrackBlockUep  );
  RegisterDataBlock("TrackBlockUmp" ,"TStnTrackBlock"   ,&fTrackBlockUmp  );

  RegisterDataBlock("ClusterBlock"  ,"TStnClusterBlock" ,&fClusterBlock);
  RegisterDataBlock("CalDataBlock"  ,"TCalDataBlock"    ,&fCalDataBlock);
  RegisterDataBlock("StrawHitBlock" ,"TStrawHitBlock"   ,&fStrawHitBlock);
  RegisterDataBlock("GenpBlock"     ,"TGenpBlock"       ,&fGenpBlock);
  RegisterDataBlock("SimpBlock"     ,"TSimpBlock"       ,&fSimpBlock);
  RegisterDataBlock("VDetBlock"     ,"TStepPointMCBlock",&fVDetBlock);

  RegisterDataBlock("CrvClusterBlock"     ,"TCrvClusterBlock"   , &fCrvClusterBlock);
  RegisterDataBlock("CrvPulseBlock"       ,"TCrvPulseBlock"     , &fCrvPulseBlock  );
//-----------------------------------------------------------------------------
// cache multiple track block pointers for convenience
//-----------------------------------------------------------------------------
  fTrackBlock[0] = fTrackBlockDem;
  fTrackBlock[1] = fTrackBlockDmm;
  fTrackBlock[2] = fTrackBlockDep;
  fTrackBlock[3] = fTrackBlockDmp;
  fTrackBlock[4] = fTrackBlockUem;
  fTrackBlock[5] = fTrackBlockUmm;
  fTrackBlock[6] = fTrackBlockUep;
  fTrackBlock[7] = fTrackBlockUmp;
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();
//-----------------------------------------------------------------------------
// initialize likelihood
//-----------------------------------------------------------------------------
  // const char   *pid_version;
  // pid_version = gEnv->GetValue("mu2e.PidVersion","_none_");
  //  fLogLH->Init(pid_version);

  return 0;
}

//-----------------------------------------------------------------------------
int TCosmicsAnaModule::BeginRun() {
  int rn = GetHeaderBlock()->RunNumber();
  TStntuple::Init(rn);
  return 0;
}

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
void TCosmicsAnaModule::BookCaloHistograms(HistBase_t* HistBase, const char* Folder) {
  //     char name [200];
  //     char title[200];

  CaloHist_t* Hist = (CaloHist_t*) HistBase;

  HBook1F(Hist->fDiskID ,"disk_id",Form("%s: Disk ID"       ,Folder), 10, 0,  10,Folder);

  for (int i=0; i<kNDisks; i++) {
    HBook1F(Hist->fEnergy  [i],Form("energy_%i",i),Form("%s: Hit Energy[%i]",Folder,i),200, 0, 100,Folder);
    HBook1F(Hist->fTime    [i],Form("time_%i"  ,i),Form("%s: Hit time  [%i]",Folder,i),200, 0,2000,Folder);
    HBook1F(Hist->fNHits   [i],Form("nhits_%i" ,i),Form("%s: NHits     [%i]",Folder,i), 50, 0,  50,Folder);
    HBook1F(Hist->fRadius  [i],Form("r_%i"     ,i),Form("%s: Radius    [%i]",Folder,i),100, 0,1000,Folder);
    HBook1F(Hist->fRadiusWE[i],Form("rwe_%i"   ,i),Form("%s: RadiusWE  [%i]",Folder,i),100, 0,1000,Folder);

    HBook1F(Hist->fE700    [i],Form("e700_%i",i),Form("%s: Hit Energy[%i] (T > 700ns)",Folder,i),200, 0, 100,Folder);
    HBook1F(Hist->fT700    [i],Form("t700_%i",i),Form("%s: Hit time  [%i] (T > 700ns)",Folder,i),200, 0,2000,Folder);
    HBook1F(Hist->fN700    [i],Form("n700_%i",i),Form("%s: NHits     [%i] (T > 700ns)",Folder,i), 50, 0,  50,Folder);

    HBook1F(Hist->fR700  [i],Form("r700_%i"  ,i),Form("%s: Radius (T>700) [%i]",Folder,i),100, 0,1000,Folder);
    HBook1F(Hist->fRWE700[i],Form("rwe700_%i",i),Form("%s: Radius*E(T>700)[%i]",Folder,i),100, 0,1000,Folder);
  }
}

//_____________________________________________________________________________
void TCosmicsAnaModule::BookHistograms() {

  //  char name [200];
  //  char title[200];

  TFolder* fol;
  TFolder* hist_folder;
  char     folder_name[200];

  DeleteHistograms();
  hist_folder = (TFolder*) GetFolder()->FindObject("Hist");

//-----------------------------------------------------------------------------
// book crystal histograms
//-----------------------------------------------------------------------------
  HBook1F(fHist.fCrystalR[0],"rc_0"     ,Form("disk [0] crystal radius"),100,0,1000,"Hist");
  HBook1F(fHist.fCrystalR[1],"rc_1"     ,Form("disk [1] crystal radius"),100,0,1000,"Hist");
//-----------------------------------------------------------------------------
// book event histograms
//-----------------------------------------------------------------------------
  int book_event_histset[kNEventHistSets];
  for (int i=0; i<kNEventHistSets; i++) book_event_histset[i] = 0;

  book_event_histset[ 0] = 1;		// all events
  book_event_histset[ 1] = 1;	        // events with a reconstructed track
  book_event_histset[ 2] = 1;	        // events without reconstructed tracks
  book_event_histset[ 3] = 1;	        // events with a reconstructed cluster
  book_event_histset[ 4] = 1;	        // events without reconstructed clusters
  book_event_histset[ 5] = 1;	        // events with DEM BESTID=0 
  book_event_histset[ 6] = 1;	        // events with DEM BESTID=0 NTRKU=0
  book_event_histset[ 7] = 1;	        // events with E(cluster) > 60 MeV
  book_event_histset[ 8] = 1;	        // events with the highest energy cluster on the 1st disk
  book_event_histset[ 9] = 1;	        // events with the highest energy cluster on the 2nd disk
  book_event_histset[10] = 0;	        // 
  book_event_histset[11] = 1;	        // events with N(straw hits by MC part) > 20
  book_event_histset[12] = 1;	        // evt_11 + (pfront > 100)
  book_event_histset[13] = 1;	        // DEM reconstructed 
  book_event_histset[14] = 1;	        // DEM IDWORD=0
  book_event_histset[15] = 1;	        // TRKPATREC or CALPATREC present 
  book_event_histset[16] = 0;	        // 
  book_event_histset[17] = 0;	        // 
  book_event_histset[18] = 0;	        // 
  book_event_histset[19] = 0;	        // 
  book_event_histset[20] = 0;	        // 
  book_event_histset[21] = 0;	        // 
  book_event_histset[22] = 0;	        // 
  book_event_histset[23] = 0;	        // 
					// TrkPatRec tracks
  book_event_histset[24] = 1;	        // events with at least one reco track
  book_event_histset[25] = 1;	        // 
  book_event_histset[26] = 1;	        // 
  book_event_histset[27] = 1;	        // 
  book_event_histset[28] = 1;	        // 

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
// book simp histograms
//-----------------------------------------------------------------------------
  int book_simp_histset[kNSimpHistSets];
  for (int i=0; i<kNSimpHistSets; i++) book_simp_histset[i] = 0;

  book_simp_histset[ 0] = 1;		// all events

  for (int i=0; i<kNSimpHistSets; i++) {
    if (book_simp_histset[i] != 0) {
      sprintf(folder_name,"sim_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fSimp[i] = new SimpHist_t;
      BookSimpHistograms(fHist.fSimp[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book track histograms
//-----------------------------------------------------------------------------
  int book_track_histset[kNTrackHistSets];
  for (int i=0; i<kNTrackHistSets; i++) book_track_histset[i] = 0;

  book_track_histset[  0] = 1;		// all tracks e-
  book_track_histset[  1] = 1;		// all tracks e- IDWORD=0
  book_track_histset[  2] = 1;		// good e- + N(upstream = 0)
  book_track_histset[  3] = 1;		// trk_2 + (nclusters>0
  book_track_histset[  4] = 1;		// trk_2 + (nclusters=0)
  book_track_histset[  5] = 1;		// trk_2 + E/P cut
  book_track_histset[  6] = 1;		// trk_2 + E/P + (chi2_tcm < 100)
  book_track_histset[  7] = 1;		// trk_2 + E/P + (chi2_tcm < 100) + (llhr > 0)
  book_track_histset[  8] = 1;		// trk_7 + e-
  book_track_histset[  9] = 1;		// trk_7 + mu-
  book_track_histset[ 10] = 1;		// trk_7 + e+
  book_track_histset[ 11] = 1;		// trk_7 + mu+
  book_track_histset[ 12] = 0;		// tracks with fcons < 1.e-4
  book_track_histset[ 13] = 0;		// "Set C" tracks with 100 <= P < 110 
  book_track_histset[ 14] = 0;		// [13] + no upstream electrons
  book_track_histset[ 15] = 0;		// [14] + calorimeter preselection
  book_track_histset[ 16] = 0;		// [15] + LLHR > 1.5
  book_track_histset[ 17] = 0;		// [16] + electrons
  book_track_histset[ 18] = 0;		// [16] + muons
  book_track_histset[ 19] = 0;		// Set C tracks with E/P > 0

  book_track_histset[ 20] = 0;		// 

  book_track_histset[ 21] = 1;		// e- all selections but no PID
  book_track_histset[ 22] = 1;		// mu-
  book_track_histset[ 23] = 1;		// e+
  book_track_histset[ 24] = 1;		// mu+

  book_track_histset[ 31] = 1;		// e- after all selections
  book_track_histset[ 32] = 1;		// mu-
  book_track_histset[ 33] = 1;		// e+
  book_track_histset[ 34] = 1;		// mu+

  book_track_histset[100] = 1.;         // DMM
  book_track_histset[101] = 1.;         // DMM

  book_track_histset[200] = 1.;		// DEP
  book_track_histset[201] = 1.;		// DEP

  book_track_histset[300] = 1.;		// DMP
  book_track_histset[301] = 1.;		// DMP

  book_track_histset[400] = 1.;		// UEM
  book_track_histset[401] = 1.;		// UEM

  book_track_histset[500] = 1.;		// UMM
  book_track_histset[501] = 1.;		// UMM

  book_track_histset[600] = 1.;		// UEP
  book_track_histset[601] = 1.;		// UEP

  book_track_histset[700] = 1.;		// UMP
  book_track_histset[701] = 1.;		// UMP

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
// book track ID histograms
//-----------------------------------------------------------------------------
  int book_track_id_histset[kNTrackIDHistSets];
  for (int i=0; i<kNTrackIDHistSets; i++) book_track_id_histset[i] = 0;

  book_track_id_histset[  0] = 1;          // all tracks on input

  for (int i=0; i<kNTrackIDHistSets; i++) {
    if (book_track_id_histset[i] != 0) {
      sprintf(folder_name,"tid_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fTrackID[i] = new TStnTrackID::Hist_t;
      BookTrackIDHistograms(fHist.fTrackID[i],Form("Hist/%s",folder_name));
    }
  }

//-----------------------------------------------------------------------------
// book cluster histograms
//-----------------------------------------------------------------------------
  int book_cluster_histset[kNClusterHistSets];
  for (int i=0; i<kNClusterHistSets; i++) book_cluster_histset[i] = 0;

  book_cluster_histset[0] = 1;		// all clusters
  book_cluster_histset[1] = 1;		// clusters in events with the reconstructed e-
  book_cluster_histset[2] = 1;		// clusters in events with the track passing SetC cuts
  book_cluster_histset[3] = 1;		// clusters in events w/track passing SetC cuts and |dt|<2.5ns 
  book_cluster_histset[4] = 1;		// clusters > 10 MeV
  book_cluster_histset[5] = 1;		// clusters > 60 MeV
  book_cluster_histset[6] = 1;		// clusters disk#0
  book_cluster_histset[7] = 1;		// clusters disk#1

  for (int i=0; i<kNClusterHistSets; i++) {
    if (book_cluster_histset[i] != 0) {
      sprintf(folder_name,"cls_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fCluster[i] = new ClusterHist_t;
      BookClusterHistograms(fHist.fCluster[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book calorimeter histograms
//-----------------------------------------------------------------------------
  int book_calo_histset[kNCaloHistSets];
  for (int i=0; i<kNCaloHistSets; i++) book_calo_histset[i] = 0;

  book_calo_histset[0] = 1;		// all crystals
  book_calo_histset[1] = 1;		// all crystals, e > 0
  book_calo_histset[2] = 1;		// all crystals, e > 0.1
  book_calo_histset[3] = 1;		// all crystals, e > 1.0

  for (int i=0; i<kNCaloHistSets; i++) {
    if (book_calo_histset[i] != 0) {
      sprintf(folder_name,"cal_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fCalo[i] = new CaloHist_t;
      BookCaloHistograms(fHist.fCalo[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book Genp histograms
//-----------------------------------------------------------------------------
  int book_genp_histset[kNGenpHistSets];
  for (int i=0; i<kNGenpHistSets; i++) book_genp_histset[i] = 0;

  book_genp_histset[0] = 1;		// all particles
//   book_genp_histset[1] = 1;		// all crystals, e > 0
//   book_genp_histset[2] = 1;		// all crystals, e > 0.1
//   book_genp_histset[3] = 1;		// all crystals, e > 1.0

  for (int i=0; i<kNGenpHistSets; i++) {
    if (book_genp_histset[i] != 0) {
      sprintf(folder_name,"gen_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fGenp[i] = new GenpHist_t;
      BookGenpHistograms(fHist.fGenp[i],Form("Hist/%s",folder_name));
    }
  }

}

//-----------------------------------------------------------------------------
void TCosmicsAnaModule::FillCaloHistograms(HistBase_t* HistBase, TStnCrystal* Cr) {

  int                    nhits;
  float                  t, e, r, e700, n700;
  TCalHitData*           hit;

  CaloHist_t* Hist = (CaloHist_t*) HistBase;
					// determine crystal coordinates
  TDisk* disk = Cr->Disk();

  int idisk = disk->SectionID();
  // time needs to be defiend
  //    t  = Cr->Time();
  e     = Cr->Energy();
  r     = Cr->Radius();
  nhits = Cr->NHits();

  Hist->fDiskID->Fill(idisk);

  Hist->fEnergy  [idisk]->Fill(e);
  Hist->fNHits   [idisk]->Fill(nhits);
  //    Hist->fTime    [idisk]->Fill(t);
  Hist->fRadius  [idisk]->Fill(r);
  Hist->fRadiusWE[idisk]->Fill(r,e);
    
  e700 = 0;
  n700 = 0;
  for (int i=0; i<nhits; i++) {
    hit  = Cr->CalHitData(i);
    t   = hit->Time();
    Hist->fTime[idisk]->Fill(t);
    if (t > 700.) {
      n700 += 1;
      e700 += hit->Energy();
      Hist->fT700[idisk]->Fill(t);
    }
  }

  Hist->fE700   [idisk]->Fill(e700);
  Hist->fN700   [idisk]->Fill(n700);

  if (n700 > 0) {
    Hist->fR700  [idisk]->Fill(r);
    Hist->fRWE700[idisk]->Fill(r,e700);
  }
}

//_____________________________________________________________________________
void TCosmicsAnaModule::FillHistograms() {

  double       cl_e(-1.);
  int          disk_id(-1), nsh;
  float        pfront; 
  TStnCluster  *cl0;

  //  cos_th = fEle->momentum().pz()/fEle->momentum().vect().mag();

  if (fNClusters > 0) {
    cl0     = fClusterBlock->Cluster(0);
    cl_e    = cl0->Energy();
    disk_id = cl0->DiskID();
  }
  //-----------------------------------------------------------------------------
  // event histograms
  //-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0],&fEvtPar);

  if (fNTracks[0]> 0) FillEventHistograms(fHist.fEvent[1],&fEvtPar);
  else                FillEventHistograms(fHist.fEvent[2],&fEvtPar);

  if (fNClusters > 0) FillEventHistograms(fHist.fEvent[3],&fEvtPar);
  else                FillEventHistograms(fHist.fEvent[4],&fEvtPar);

  if (fNGoodTracks[0] > 0) { 
    FillEventHistograms(fHist.fEvent[5],&fEvtPar); 

    if ((fNTrkUNeg == 0) && (fNTrkDPos == 0)) {
      FillEventHistograms(fHist.fEvent[6],&fEvtPar); 
    }
  }

  if (cl_e > 60.) {
    FillEventHistograms(fHist.fEvent[7],&fEvtPar); 
    if (GetDebugBit(34)) {
      if (fNTracks[0] <= 0) {
	GetHeaderBlock()->Print(Form(" bit:034 cl_e = %10.3f",cl_e));
      }
    }
  }

  if      (disk_id == 0) FillEventHistograms(fHist.fEvent[8],&fEvtPar);
  else if (disk_id == 1) FillEventHistograms(fHist.fEvent[9],&fEvtPar);
//-----------------------------------------------------------------------------
// Dave's ladder for all tracks
// 1. N(straw hits) > 20
//-----------------------------------------------------------------------------
  if (fSimPar.fParticle) {
    nsh    = fSimPar.fParticle->NStrawHits();
    pfront = fSimPar.fParticle->fMomTrackerFront;
  }
  else {
    nsh    = -1;
    pfront = -1.e6;
  }
  
  if (nsh >= 20) {
    FillEventHistograms(fHist.fEvent[11],&fEvtPar);
    if (pfront > 100.) {
      FillEventHistograms(fHist.fEvent[12],&fEvtPar);
      
      if (fNTracks[0] > 0) {
//-----------------------------------------------------------------------------
// DEM track is reconstructed
//-----------------------------------------------------------------------------
	FillEventHistograms(fHist.fEvent[13],&fEvtPar);

	TrackPar_t* tp = &fTrackPar[0][0];
	
	if (tp->fIDWord[fBestID] == 0) {
	  FillEventHistograms(fHist.fEvent[14],&fEvtPar);
	}
	
      }
    }
  }
//-----------------------------------------------------------------------------
// Simp histograms
//-----------------------------------------------------------------------------
  if (fSimPar.fParticle) {
    FillSimpHistograms(fHist.fSimp[0],fSimPar.fParticle);
  }
//-----------------------------------------------------------------------------
// track histograms, fill them only for the downstream e- hypothesis
//-----------------------------------------------------------------------------
  TStnTrack*   trk;
  TrackPar_t*  tp;

  for (int i=0; i<fNTracks[0]; ++i ) {
    trk = fTrackBlockDem->Track(i);
    tp  = fTrackPar[0]+i;

    FillTrackHistograms(fHist.fTrack[0],trk,tp,&fSimPar);

    fTrackID[fBestID]->FillHistograms(fHist.fTrackID[0],trk);

    if (tp->fIDWord[fBestID] == 0) {
					// GOOD track: IDWORD=0

      FillTrackHistograms(fHist.fTrack[1],trk,tp,&fSimPar);

      if ((fNTrkUNeg == 0) && (fNTrkDPos == 0)) {
//-----------------------------------------------------------------------------
// no other track which couldn't be mistaken with the DEM
//----------------------------------------------------------------------------- 
	FillTrackHistograms(fHist.fTrack[2],trk,tp,&fSimPar);

	if   (fNClusters > 0) FillTrackHistograms(fHist.fTrack[3],trk,tp,&fSimPar);
	else                  FillTrackHistograms(fHist.fTrack[4],trk,tp,&fSimPar);
//-----------------------------------------------------------------------------
// TRK 5 : add E/P cut
//-----------------------------------------------------------------------------
	if ((tp->fEp > 0.4) && ( tp->fEp < 1.15)) {
	  FillTrackHistograms(fHist.fTrack[5],trk,tp,&fSimPar);

	  if (tp->fChi2Tcm < 100) {
	    FillTrackHistograms(fHist.fTrack[6],trk,tp,&fSimPar);

	    double llhr_cal = trk->LogLHRCal();
	    
	    if      (trk->fPdgCode ==  11) FillTrackHistograms(fHist.fTrack[21],trk,tp,&fSimPar); // e-
	    else if (trk->fPdgCode ==  13) {
	      FillTrackHistograms(fHist.fTrack[22],trk,tp,&fSimPar); // mu-
	      if (tp->fDt > 0) {
		if (GetDebugBit(6)) {
		  GetHeaderBlock()->Print(Form("bit006: mu- with DT = %10.3f",tp->fDt));
		}
	      }
	    }
	    else if (trk->fPdgCode == -11) FillTrackHistograms(fHist.fTrack[23],trk,tp,&fSimPar); // e+
	    else if (trk->fPdgCode == -13) FillTrackHistograms(fHist.fTrack[24],trk,tp,&fSimPar); // mu+

	    if (llhr_cal > 2.) {
	      FillTrackHistograms(fHist.fTrack[7],trk,tp,&fSimPar);

	      if      (trk->fPdgCode ==  11) FillTrackHistograms(fHist.fTrack[31],trk,tp,&fSimPar); // e-
	      else if (trk->fPdgCode ==  13) FillTrackHistograms(fHist.fTrack[32],trk,tp,&fSimPar); // mu-
	      else if (trk->fPdgCode == -11) FillTrackHistograms(fHist.fTrack[33],trk,tp,&fSimPar); // e+
	      else if (trk->fPdgCode == -13) {
		FillTrackHistograms(fHist.fTrack[34],trk,tp,&fSimPar); // mu+
		if (GetDebugBit(4)) {
		  GetHeaderBlock()->Print(Form("bit004: mu+ passed all selections"));
		}
	      }
	    }
	  }
	}
      }
    }
  }
//-----------------------------------------------------------------------------
// distributions for tracks from different blocks
//-----------------------------------------------------------------------------
  for (int ib=1; ib<kNTrackBlocks; ib++) {
    int loc = 100*ib;
    for (int it=0; it<fNTracks[ib]; it++) {
      trk = fTrackBlock[ib]->Track(it);
      tp  = &fTrackPar[ib][it];

      FillTrackHistograms(fHist.fTrack[loc+0],trk,tp,&fSimPar);
      if (tp->fIDWord[fBestID] == 0) {
	FillTrackHistograms(fHist.fTrack[loc+1],trk,tp,&fSimPar);
      }
    }
  }
//-----------------------------------------------------------------------------
// cluster histograms
//-----------------------------------------------------------------------------
  TStnCluster* cl;
  int id;
  for (int i=0; i<fNClusters; ++i ) {
    cl = fClusterBlock->Cluster(i);
    id = cl->DiskID();
    FillClusterHistograms(fHist.fCluster[0],cl);

    if (fNTracks[0]        >  0 ) FillClusterHistograms(fHist.fCluster[1],cl);
    if (fNGoodTracks[0]    >  0 ) FillClusterHistograms(fHist.fCluster[2],cl);
    if (fNMatchedTracks[0] >  0 ) FillClusterHistograms(fHist.fCluster[3],cl);
    if (cl->Energy()       > 10.) FillClusterHistograms(fHist.fCluster[4],cl);
    if (cl->Energy()       > 60.) FillClusterHistograms(fHist.fCluster[5],cl);

    if      (id == 0         ) FillClusterHistograms(fHist.fCluster[6],cl);
    else if (id == 1         ) FillClusterHistograms(fHist.fCluster[7],cl);
  }
//-----------------------------------------------------------------------------
// calorimeter histograms
//-----------------------------------------------------------------------------
  TDisk*         disk;
  TStnCrystal*   cr;

  if (fCalorimeterType == 2) {
    int nd = fDiskCalorimeter->NDisks();

    for (int i=0; i<nd; i++) {
      disk = fDiskCalorimeter->Disk(i);
      for (int ic=0; ic<disk->NCrystals(); ic++) {
	cr = disk->Crystal(ic);
	FillCaloHistograms(fHist.fCalo[0],cr);

	if (cr->Energy() > 0  ) FillCaloHistograms(fHist.fCalo[1],cr);
	if (cr->Energy() > 0.1) FillCaloHistograms(fHist.fCalo[2],cr);
	if (cr->Energy() > 1.0) FillCaloHistograms(fHist.fCalo[3],cr);
      }
    }
  }
//-----------------------------------------------------------------------------
// radial distributions for crystals
//-----------------------------------------------------------------------------
  static int first_entry(1);

  if (first_entry == 1) {
    if (fCalorimeterType == 2) {
      int nd = fDiskCalorimeter->NDisks();
	
      for (int i=0; i<nd; i++) {
	disk = fDiskCalorimeter->Disk(i);
	for (int ic=0; ic<disk->NCrystals(); ic++) {
	  cr = disk->Crystal(ic);

	  fHist.fCrystalR[i]->Fill(cr->Radius());
	}
      }
    }
  }

//-----------------------------------------------------------------------------
// fill GENP histograms
// GEN_0: all particles
//-----------------------------------------------------------------------------
  TGenParticle* genp;
  for (int i=0; i<fNGenp; i++) {
    genp = fGenpBlock->Particle(i);
    FillGenpHistograms(fHist.fGenp[0],genp);
  }
  first_entry = 0;
}



//-----------------------------------------------------------------------------
// 2014-04-30: it looks that reading the straw hits takes a lot of time - 
//              turn off by default by commenting it out
//-----------------------------------------------------------------------------
int TCosmicsAnaModule::Event(int ientry) {

  //  double                p;
  TLorentzVector        mom;

  TDiskCalorimeter::GeomData_t disk_geom;

  for (int i=0; i< kNTrackBlocks; i++) fTrackBlock[i]  ->GetEntry(ientry);

  fClusterBlock->GetEntry(ientry);
  fCalDataBlock->GetEntry(ientry);
  fGenpBlock->GetEntry(ientry);
  fSimpBlock->GetEntry(ientry);
  fVDetBlock->GetEntry(ientry);

  fEvtPar.fNCrvClusters     = fCrvClusterBlock->NClusters();
  fEvtPar.fNCrvPulses       = fCrvPulseBlock->NPulses();
  fEvtPar.fNCrvCoincidences = fCrvPulseBlock->NCoincidences();
//-----------------------------------------------------------------------------
// assume electron in the first particle, otherwise the logic will need to 
// be changed
//-----------------------------------------------------------------------------
  // fNGenp            = fGenpBlock->NParticles();
  fEvtPar.fNSimp    = fSimpBlock->NParticles();

  for (int i=fEvtPar.fNSimp-1; i>=0; i--) {
    TSimParticle* simp = fSimpBlock->Particle(i);
    int pdg_code       = simp->PDGCode();
    int generator_id   = simp->GeneratorID();         // MC process code

    if ((pdg_code == fPDGCode) and (generator_id == fMCProcessCode)) {
      fEvtPar.fSimp  = simp;
      fEvtPar.fPartE = simp->StartMom()->Energy();
    }
  }
					// may want to revisit the definition of fSimp

  fSimPar.fParticle = fEvtPar.fSimp;
//-----------------------------------------------------------------------------
// process virtual detectors - for fSimp need parameters at tracker entrance
//-----------------------------------------------------------------------------
  fSimPar.fTFront   = NULL;
  fSimPar.fTMid     = NULL;
  fSimPar.fTBack    = NULL;

  int nvdhits = fVDetBlock->NStepPoints();
  for (int i=0; i<nvdhits; i++) {
    TStepPointMC* vdhit = fVDetBlock->StepPointMC(i);
    if (vdhit->PDGCode() == fSimPar.fParticle->fPdgCode) {
      if ((vdhit->VolumeID() == 13) || (vdhit->VolumeID() == 14)) {
	fSimPar.fTFront = vdhit;
      }
      else if ((vdhit->VolumeID() == 11) || (vdhit->VolumeID() == 12)) {
	fSimPar.fTMid = vdhit;
      }
      else if (vdhit->VolumeID() == mu2e::VirtualDetectorId::TT_Back) {
	fSimPar.fTBack = vdhit;
      }
    }
  }
					// this is a kludge, to be removed at the next 
					// ntupling 

  if (fDiskCalorimeter->Initialized() == 0) {
    disk_geom.fNDisks = fCalDataBlock->NDisks();

    for (int i=0; i<disk_geom.fNDisks; i++) {
      //      disk_geom.fNCrystals[i] = fCalDataBlock->fNCrystals[i];
      disk_geom.fRMin[i]      = fCalDataBlock->fRMin[i];
      disk_geom.fRMax[i]      = fCalDataBlock->fRMax[i];
      disk_geom.fZ0  [i]      = fCalDataBlock->fZ0  [i];
    }

    disk_geom.fHexSize          = fCalDataBlock->CrystalSize()*2;
    // kludge , so far
    disk_geom.fMinFraction      = 1.; // fCalDataBlock->MinFraction();
    disk_geom.fWrapperThickness = fCalDataBlock->WrapperThickness();
    disk_geom.fShellThickness   = fCalDataBlock->ShellThickness();

    fDiskCalorimeter->Init(&disk_geom);
  }

  fNCalHits   = fCalDataBlock->NHits();
  fNStrawHits = fStrawHitBlock->NHits();

  fDiskCalorimeter->InitEvent(fCalDataBlock);

  fNHyp       = -1;
  fBestHyp[0] = -1;
  fBestHyp[1] = -1;
//-----------------------------------------------------------------------------
// initialize additional track parameters
//-----------------------------------------------------------------------------
  for (int i=0; i<kNTrackBlocks; i++) {
    fNTracks[i] = fTrackBlock[i]->NTracks();
    for (int it=0; it<fNTracks[i]; it++) {
      TrackPar_t* tp  = fTrackPar[i]+it;
      tp->fTrackID[0] = TAnaModule::fTrackID_BOX;
      tp->fTrackID[1] = TAnaModule::fTrackID_MVA;
      //      tp->fLogLH      = TAnaModule::fLogLH;
    }
    InitTrackPar(fTrackBlock[i],fClusterBlock,fTrackPar[i],&fSimPar);
  }

  fNTrkDNeg     = fNTracks[kDem]+fNTracks[kDmm];
  fNTrkDPos     = fNTracks[kDep]+fNTracks[kDmp];
  fNTrkUNeg     = fNTracks[kUem]+fNTracks[kUmm];
  fNTrkUPos     = fNTracks[kUep]+fNTracks[kUmp];

  fNTrkUpstream = fNTrkUNeg+fNTrkUPos;

  if (fNTracks[0] == 0) fTrack = 0;
  else                  fTrack = fTrackBlockDem->Track(0);

  fNClusters = fClusterBlock->NClusters();
  if (fNClusters == 0) fCluster = 0;
  else                 fCluster = fClusterBlock->Cluster(0);

  //  fDiskCalorimeter->InitEvent(fCalDataBlock);

  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TCosmicsAnaModule::Debug() {

//   TStnTrack* trk;
//   TrackPar_t* tp;

  char s[1000];

//-----------------------------------------------------------------------------
// bit 3: Set C tracks with large DX : 70mm < |DX| < 90mm
//-----------------------------------------------------------------------------
  if (GetDebugBit(3) == 1) {
    sprintf(s,"NTracks[0,1,2,3]: %3i %3i %3i %3i",
	    fNTracks[0],fNTracks[1],fNTracks[2],fNTracks[3]);
	  
    GetHeaderBlock()->Print(Form(":bit003: %s",s));
  }

  int ntrk = fTrackBlockDem->NTracks();

  for (int itrk=0; itrk<ntrk; itrk++) {
    TStnTrack* trk  = fTrackBlockDem->Track(itrk);
     TrackPar_t* tp = &fTrackPar[0][itrk];

     if (GetDebugBit(5) && ((tp->fIDWord[fBestID] & ~TStnTrackID::kTanDipBit) == 0)) {
       GetHeaderBlock()->Print(Form("trk->TanDip = %10.3f",trk->TanDip()));
     }

     if (GetDebugBit(6)) {
       GetHeaderBlock()->Print(Form("trk->PDG: %7i trk->TanDip = %10.3f tp->fLogLHDedm: %10.5f",
				    trk->fPdgCode,trk->TanDip(),tp->fLogLHDedm));
     }
  }

}

//_____________________________________________________________________________
int TCosmicsAnaModule::EndJob() {
  printf("----- end job: ---- %s\n",GetName());
  return 0;
}

//_____________________________________________________________________________
void TCosmicsAnaModule::Test001() {

  // mu2e::HexMap* hmap      = new mu2e::HexMap();

  // mu2e::HexLK hex_index(0,0);

  // for (int i=0; i<40; i++) {
  //   hex_index = hmap->lk(i);
  //   printf(" i,l,k = %5i %5i %5i\n",i,hex_index._l,hex_index._k);
  // }
}

}
