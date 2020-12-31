//////////////////////////////////////////////////////////////////////////////
// use of TStnTrack::fTmp:
//
// Tmp(0) : corrected momentum at the tracker front - not yet
// 
// use of debug bits: bits 0-2 are reserved
//  0  : all events
//  1  : passed events
//  2  : rejected events
// 
//  3  : events with set C tracks and 70mm < |dx|  < 90 mm
//  4  : events with DpF > 1 MeV : obviously, misreconstructed ones
//  5  : events with N(tracks) > 1
//  6  : events trk_41 with 0.8< E/P < 1.1 - tracks missed by CalPatRec
//  7  : events (muo) with LogLHRCal >   20
//  8  : events (ele) with LogLHRCal < - 20
//  9  : events (muo) with 0.42 < E/P < 0.46
// 10  : events (muo) with Set C track with ECL > 80 MeV
// 28  : Set C DEM tracks with E/P > 1.1
// 29  : TRK_19 (Set C DEM tracks with a cluster) and LLHR(cal) < 0
// 31  : EVT_6 events with ce_costh > 0.8 
// 32  : TRK_1 events with chi2tcm > 100. 
// 33  : DU < -80mm - study edge effects
// 34  : EVT_7: events with E_CL > 60 and no CalPatRec tracks 
// 35  : TRK_1: events with P > 106 MeV/c - misreconstruction
// 36  : TRK_23 events with P < 80: odd misidentified muons - turned out to be DIO electrons
// 37  : TRK_26 LLHR_CAL > 5
// 38  : EVT_7: events with E_CL > 60 and no CalPatRec tracks and TrkPatRec
// 39  : trk_1: events with |SIN_TC| > 0.6
// 40  : EVT_7: events with E_CL > 60 and no tracks at all
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
#include "Stntuple/geom/TDisk.hh"
#include "Stntuple/obj/TStnNode.hh"
#include "Stntuple/val/stntuple_val_functions.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"

//------------------------------------------------------------------------------
// Mu2e offline includes
//-----------------------------------------------------------------------------

#include "murat/ana/TTrackAnaModule.hh"

ClassImp(murat::TTrackAnaModule)

namespace murat {
//-----------------------------------------------------------------------------
TTrackAnaModule::TTrackAnaModule(const char* name, const char* title): TAnaModule(name,title)
{
  fTrackBlockName         = "TrackBlockPar";
  fTrackStrawHitBlockName = "TrackStrawHitBlock";
  fTrackNumber.Set(100);

  fDiskCalorimeter = new TDiskCalorimeter();
  fCalorimeterType = 2;

  fMinT0           = 700; 
  fMinETrig        = 50.;               // MeV
					// track-cluster matching timing cut
  fMinDtTcm        = -5.;
  fMaxDtTcm        =  8.;
//-----------------------------------------------------------------------------
// initialize Track ID
// 0: SetC  1: TrkQual>0.1 2:TrkQual>0.4
// what about number of hits ? - 3: no cuts on the number of hits
//-----------------------------------------------------------------------------
  fNID             = 2;
  fTrackID[0]      = TAnaModule::fTrackID_BOX;
  fTrackID[1]      = TAnaModule::fTrackID_MVA;

  fBestID          = 1;
//-----------------------------------------------------------------------------
// MC truth: define which MC particle to consider as signal
// 2:conversionGun, 28:StoppedParticleReactionGun - see 
//-----------------------------------------------------------------------------
  fPdgCode       = 11;
  fGeneratorCode =  2;
  fDirection     =  1;

  fApplyCorrections = 0;
}

//-----------------------------------------------------------------------------
TTrackAnaModule::~TTrackAnaModule() {
  for (int i=0; i<fNID; i++) delete fTrackID[i];
}


//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TTrackAnaModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
  RegisterDataBlock(fTrackBlockName.Data()        ,"TStnTrackBlock"     ,&fTrackBlock  );
  RegisterDataBlock(fTrackStrawHitBlockName.Data(),"TTrackStrawHitBlock",&fTrackStrawHitBlock);

  RegisterDataBlock("ClusterBlock"        ,"TStnClusterBlock" ,&fClusterBlock);
  RegisterDataBlock("CalDataBlock"        ,"TCalDataBlock"    ,&fCalDataBlock);
  RegisterDataBlock("StrawDataBlock"      ,"TStrawDataBlock"  ,&fStrawDataBlock);
  RegisterDataBlock("GenpBlock"           ,"TGenpBlock"       ,&fGenpBlock);
  RegisterDataBlock("SimpBlock"           ,"TSimpBlock"       ,&fSimpBlock);
  RegisterDataBlock("SpmcBlockVDet"       ,"TStepPointMCBlock",&fVDetBlock);
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();

  return 0;
}


//_____________________________________________________________________________
int TTrackAnaModule::BeginRun() {
  int rn = GetHeaderBlock()->RunNumber();
  TStntuple::Init(rn);
  return 0;
}


//-----------------------------------------------------------------------------
void TTrackAnaModule::BookCaloHistograms(CaloHist_t* Hist, const char* Folder) {
  //     char name [200];
  //     char title[200];
  //-----------------------------------------------------------------------------
  //  
  //-----------------------------------------------------------------------------
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

//-----------------------------------------------------------------------------
void TTrackAnaModule::BookClusterHistograms(ClusterHist_t* Hist, const char* Folder) {
//   char name [200];
//   char title[200];

  HBook1F(Hist->fDiskID ,"disk_id",Form("%s: Disk ID"       ,Folder), 10, 0,  10,Folder);
  HBook1F(Hist->fEnergy ,"energy" ,Form("%s: Cluster Energy",Folder),500, 0, 250,Folder);
  HBook1F(Hist->fT0     ,"t0"     ,Form("%s: cluster T0"    ,Folder),200, 0,2000,Folder);
  HBook1F(Hist->fRow    ,"row"    ,Form("%s: cluster Row"   ,Folder),200, 0, 200,Folder);
  HBook1F(Hist->fCol    ,"col"    ,Form("%s: cluster column",Folder),200, 0, 200,Folder);
  HBook1F(Hist->fX      ,"x"      ,Form("%s: cluster X"     ,Folder),200, -5000,5000,Folder);
  HBook1F(Hist->fY      ,"y"      ,Form("%s: cluster Y"     ,Folder),200,-1000,1000,Folder);
  HBook1F(Hist->fZ      ,"z"      ,Form("%s: cluster Z"     ,Folder),200, 11500,13500,Folder);
  HBook1F(Hist->fR      ,"r"      ,Form("%s: cluster Radius",Folder),100, 0,  1000,Folder);
  HBook1F(Hist->fYMean  ,"ymean"  ,Form("%s: cluster YMean" ,Folder),400,-200,200,Folder);
  HBook1F(Hist->fZMean  ,"zmean"  ,Form("%s: cluster ZMean" ,Folder),400,-200,200,Folder);
  HBook1F(Hist->fSigY   ,"sigy"   ,Form("%s: cluster SigY"  ,Folder),100, 0,100,Folder);
  HBook1F(Hist->fSigZ   ,"sigz"   ,Form("%s: cluster SigZ"  ,Folder),100, 0,100,Folder);
  HBook1F(Hist->fSigR   ,"sigr"   ,Form("%s: cluster SigR"  ,Folder),100, 0,100,Folder);
  HBook1F(Hist->fNCr0   ,"ncr0"   ,Form("%s: cluster NCR[0]",Folder),100, 0,100,Folder);
  HBook1F(Hist->fNCr1   ,"ncr1"   ,Form("%s: cluster NCR[1]",Folder),100, 0,100,Folder);
  HBook1F(Hist->fFrE1   ,"fre1"   ,Form("%s: E1/Etot"       ,Folder),200, 0,  1,Folder);
  HBook1F(Hist->fFrE2   ,"fre2"   ,Form("%s: (E1+E2)/Etot"  ,Folder),200, 0,  1,Folder);
  HBook1F(Hist->fSigE1  ,"sige1"   ,Form("%s: SigmaE/Etot"  ,Folder),200, 0, 10,Folder);
  HBook1F(Hist->fSigE2  ,"sige2"   ,Form("%s: SigmaE/Emean" ,Folder),200, 0, 10,Folder);
}


//_____________________________________________________________________________
void TTrackAnaModule::BookHistograms() {

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
  book_event_histset[ 5] = 1;	        // events w/o reconstructed tracks, |costh|<0.4
  book_event_histset[ 6] = 1;	        // events with tracks passing "Set C" cuts
  book_event_histset[ 7] = 1;	        // events with E(cluster) > 60 MeV
  book_event_histset[ 8] = 1;	        // events with the highest energy cluster on the 1st disk
  book_event_histset[ 9] = 1;	        // events with the highest energy cluster on the 2nd disk
  book_event_histset[10] = 1;	        // events with SetC tracks 103.5 < p < 105.

  book_event_histset[11] = 1;	        // Track ID ladder : 11-20
  book_event_histset[12] = 1;	        // 
  book_event_histset[13] = 1;	        // 
  book_event_histset[14] = 1;	        // 
  book_event_histset[15] = 1;	        // 
  book_event_histset[16] = 1;	        // 
  book_event_histset[17] = 1;	        // 
  book_event_histset[18] = 1;	        // 
  book_event_histset[19] = 1;	        // 
  book_event_histset[20] = 1;	        // 

  book_event_histset[21] = 0;	        // 
  book_event_histset[22] = 0;	        // 
  book_event_histset[23] = 0;	        // 

					// TrkPatRec tracks
  book_event_histset[24] = 1;	        // events with at least one reco track
  book_event_histset[25] = 1;	        // 
  book_event_histset[26] = 1;	        // 
  book_event_histset[27] = 1;	        // 
  book_event_histset[28] = 1;	        // 

  book_event_histset[31] = 1;		// EFFICIENCY - Dave's definition
  book_event_histset[32] = 1;
  book_event_histset[33] = 1;
  

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
  book_track_histset[  1] = 1;		// all tracks e- passing Set C cuts 
  book_track_histset[  2] = 1;		// all tracks e- passing Set C cuts, events with clusters 
  book_track_histset[  3] = 1;		// all tracks e- passing Set C cuts, events w/o  clusters
  book_track_histset[  4] = 1;		// all tracks e- passing Set C cuts, events with clusters, no closest
  book_track_histset[  5] = 1;		// all tracks e- passing Set C cuts, |dt| < 2.5ns
  book_track_histset[  6] = 1;		// all tracks e-  E/P > 0.4
  book_track_histset[  7] = 1;		// all tracks e- passing Set C cuts, E/P > 0.4
  book_track_histset[  8] = 1;		// Set C tracks e- , |xslope| < 3.
  book_track_histset[  9] = 1;		// all  tracks in the event when there is no EM clusters E > 60 MeV
  book_track_histset[ 10] = 1;		// Set C tracks in the event when there is no EM clusters E > 60 MeV
  book_track_histset[ 11] = 1;		// all tracks with P > 103.5
  book_track_histset[ 12] = 1;		// tracks with fcons < 1.e-4
  book_track_histset[ 13] = 1;		// SetC tracks with 100 <= P < 110 
  book_track_histset[ 14] = 1;		// tracks with fcons < 1.e-2
  book_track_histset[ 15] = 1;		// all tracks intersecting the 1st disk
  book_track_histset[ 16] = 1;		// all tracks intersecting the 2nd disk
  book_track_histset[ 17] = 1;		// all tracks with no calorimeter intersections
  book_track_histset[ 18] = 1;		// Set C tracks, T0 > 700
  book_track_histset[ 19] = 1;		// Set C tracks with E/P > 0

  book_track_histset[ 20] = 1;		// BEST_ID tracks intersecting the 1st disk 
  book_track_histset[ 21] = 1;		// BEST_ID tracks intersecting the 2nd disk 
  book_track_histset[ 22] = 1;		// BEST_ID tracks with E/P > 0 and chi2(match) < 100
  book_track_histset[ 23] = 1;		// BEST_ID tracks with E/P > 0 and chi2(match) < 100 and LLHR(cal) > 0 (interesting for muons)
  book_track_histset[ 24] = 1;		// BEST_ID tracks with E/P > 0 and chi2(match) < 100 and LLHR(cal) < 0 (interesting for electrons)

  book_track_histset[ 25] = 1;		// BEST_ID tracks, 100<p<110, 0<E/P<1.15, chi2tcm<100, -5<DT<8, 
  book_track_histset[ 26] = 1;		// [25] + LLHR_CAL > 0 - interesting for muons
  book_track_histset[ 27] = 1;		// [25] + LLHR_CAL < 0 - interesting for electrons
  book_track_histset[ 28] = 1;		// Set C tracks, E/P > 1.1
  book_track_histset[ 29] = 1;		// Set C tracks 100 < P < 110 , 0 < E/P < 1.15 *precursor for TRK_32*

  book_track_histset[ 30] = 1;		// tracks with Nhits >= 25
  book_track_histset[ 31] = 1;		// tracks with Nhits >= 25 and chi/Ndof < 3
  book_track_histset[ 32] = 1;		// Set C tracks 100 < P < 110 , 0 < E/P < 1.15, chi2tcm<100

  book_track_histset[ 33] = 1;		// DaveTrkQual tracks Q > 0.4
  book_track_histset[ 34] = 1;		// DaveTrkQual tracks Q > 0.1

  book_track_histset[ 40] = 1;		// all tracks, alg_mask = 1
  book_track_histset[ 41] = 1;		// Set "C" tracks, alg_mask = 1
  book_track_histset[ 42] = 1;		// Set "C" tracks, alg_mask = 1, T > 700

  book_track_histset[ 50] = 1;		// all tracks, alg_mask = 2
  book_track_histset[ 51] = 1;		// Set "C" tracks, alg_mask = 2
  book_track_histset[ 52] = 1;		// Set "C" tracks, alg_mask = 2, T > 700

  book_track_histset[ 60] = 1;		// all tracks, alg_mask = 3
  book_track_histset[ 61] = 1;		// Set "C" tracks, alg_mask = 3
  book_track_histset[ 62] = 1;		// Set "C" tracks, alg_mask = 3, T > 700
  book_track_histset[ 63] = 1;          // best ID, best = 0 (TrkPatRec)
  book_track_histset[ 64] = 1;          // best ID, best = 1 (CalPatRec)
  book_track_histset[ 65] = 1;          // tp->fDtClZ0 > -0.6
  book_track_histset[ 66] = 0;
  book_track_histset[ 67] = 0;
  book_track_histset[ 68] = 0;
  book_track_histset[ 69] = 0;
  book_track_histset[ 70] = 0;

  book_track_histset[ 71] = 1;		// Set "C" tracks   103.5 < P < 105 : interesting for DIO

  book_track_histset[ 81] = 1;
  book_track_histset[ 82] = 1;
  book_track_histset[ 83] = 1;
  book_track_histset[ 84] = 1;
  book_track_histset[ 85] = 1;
  book_track_histset[ 86] = 1;
  

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
// need MC truth branch
//-----------------------------------------------------------------------------
// void TTrackAnaModule::FillEventHistograms(EventHist_t* Hist) {
//   double            cos_th, dio_wt, xv(-1.e6), yv(-1.e6), rv(-1.e6), zv(-1.e6), p;
//   double            e, m, r;
//   TLorentzVector    mom(1.,0.,0.,0);

//   if (fParticle) { 
//     fParticle->Momentum(mom);
//     xv = fParticle->Vx()+3904.;
//     yv = fParticle->Vy();
//     rv = sqrt(xv*xv+yv*yv);
//     zv = fParticle->Vz();
//   }

//   p      = mom.P();
//   cos_th = mom.Pz()/p;
//   dio_wt = TStntuple::DioWeightAl(p);

//   Hist->fLumWt->Fill(fLumWt);
//   Hist->fRv->Fill(rv);
//   Hist->fZv->Fill(zv);
//   Hist->fEleMom->Fill(p);
//   Hist->fDioMom->Fill(p,dio_wt);
//   Hist->fEleCosTh->Fill(cos_th);

//   Hist->fNClusters->Fill(fNClusters);
//   Hist->fNTracks->Fill  (fNTracks[0]);

//   int nsh = GetHeaderBlock()->fNStrawHits;
//   Hist->fNStrawHits[0]->Fill(nsh);
//   Hist->fNStrawHits[1]->Fill(nsh);

//   double emax   = -1;
//   double t0_cls = -1;
//   double dt     = 9999.;

//   TStnCluster* cluster(0);
//   if (fNClusters > 0) cluster = fClusterBlock->Cluster(0);

//   TStnTrack* track(0);
//   if (fNTracks > 0) track = fTrackBlock->Track(0);

//   if (cluster) {
//     emax   = cluster->Energy();
//     t0_cls = cluster->Time();
//   }

//   double t0_trk = -1;
//   int    alg_mask(-1);
//   if (track) {
//     t0_trk   = track->fT0;
//     alg_mask = track->AlgMask();
//   }

//   if (track && cluster) {
//     dt = t0_cls-t0_trk;
//   }

//   Hist->fAlgMask->Fill(alg_mask);
//   Hist->fDtClT->Fill(dt);
//   Hist->fEMax->Fill(emax);

//   TStrawHitData*  sh;
//   int n_good_hits = 0;
//   int nstraw_hits = fStrawDataBlock->NHits();
//   for (int i=0; i<nstraw_hits; i++ ) {
//     sh = fStrawDataBlock->Hit(i);
//     dt = t0_cls-sh->Time() + 15;
//     Hist->fDtClS->Fill(dt);
//     Hist->fSHTime->Fill(sh->Time());

//     if (fabs(dt+15.)< 50) n_good_hits += 1;
//   }

//   Hist->fNGoodSH->Fill(n_good_hits);

//   Hist->fNHyp->Fill(fNHyp);
//   Hist->fBestHyp[0]->Fill(fBestHyp[0]);
//   Hist->fBestHyp[1]->Fill(fBestHyp[1]);
//   Hist->fNGenp->Fill(fNGenp);
// //-----------------------------------------------------------------------------
// // crystals - count crystals with E > 1MeV
// //-----------------------------------------------------------------------------
//   int          ndisks, n_hit_crystals[kNDisks], n_hit_crystals_tot;
//   double       etot[kNDisks];
//   //  TCalHitData* cch;

//   //  int n_cch_1mev = 0;

//   //  int   nhits_vane[2][kNDisks], nhits_row [2][20], nhits_col[2][50];
//   //  int   crystal_id, vane_id, local_id, vane_row, vane_col;
//   int   bin, hit_id, idisk, nhits;
//   int   nhits_r[kNDisks][100], n_hit_crystals_r[kNDisks][100];

//   if (fCalorimeterType == 2) {
// //-----------------------------------------------------------------------------
// // disk calorimeter
// //-----------------------------------------------------------------------------

//     TCalHitData* hit;

//     ndisks = fDiskCalorimeter->NDisks();

//     if (ndisks > 0) {

//       for (int id=0; id<kNDisks; id++) {
// 	n_hit_crystals[id] = 0;
// 	etot[id]           = 0;
	
// 	for (int ib=0; ib<100; ib++) {
// 	  nhits_r         [id][ib] = 0;
// 	  n_hit_crystals_r[id][ib] = 0;
// 	}
//       }

//       nhits = fCalDataBlock->NHits();

//       for (int i=0; i< nhits; i++) {
// 	hit    = fCalDataBlock->CalHitData(i);
	
// 	hit_id = hit->ID();
// 	idisk  = fDiskCalorimeter->DiskNumber(hit_id);
// 	r      = fDiskCalorimeter->CrystalRadius(hit_id);
// 	e      = hit->Energy(); 
      
// 	etot          [idisk] += e;
// 	n_hit_crystals[idisk] += 1;
	
// 	Hist->fECrVsR[idisk]->Fill(r,e);
// 	Hist->fNCrVsR[idisk]->Fill(r,1);
	
// 	bin  = (int) (r/10.);
	
// 	nhits_r         [idisk][bin] += 1;
// //-----------------------------------------------------------------------------
// // this is not correct, one needs to check whether this crystal has been hit,
// // for the moment, to get going, ignore that
// //-----------------------------------------------------------------------------
// 	n_hit_crystals_r[idisk][bin] += 1;
//       }
//     }

//     n_hit_crystals_tot = 0;

//     double ecal = 0;
//     for (int id=0; id<ndisks; id++) {
//       n_hit_crystals_tot += n_hit_crystals[id];
//       ecal += etot[id];
// //-----------------------------------------------------------------------------
// // fill 'per-disk' histograms
// //-----------------------------------------------------------------------------
//       Hist->fETot[id]->Fill(etot[id]);

// //-----------------------------------------------------------------------------
// // 100 is fixed by the number of bins in the radial distributions
// //-----------------------------------------------------------------------------
//       for (int ib=0; ib<100; ib++) {
// 	r = (ib+0.5)*10.;
// 	Hist->fNCrystalHitsVsR[id]->Fill(r,nhits_r         [id][ib]);
// 	Hist->fNHitCrystalsVsR[id]->Fill(r,n_hit_crystals_r[id][ib]);
//       }
//     }

//     Hist->fNHitCrystalsTot->Fill(n_hit_crystals_tot);
//     Hist->fECal->Fill(ecal);

//     double ekin(-1.);
// //-----------------------------------------------------------------------------
// // there is an inconsistency in the SIMP block filling - in Mu2e offline 
// // the particle momentumis is kept in MeV/c, while the PDG mass  -in GeV/c^2..
// // thus the energy is screwed up... kludge around
// // assign muon mass
// //-----------------------------------------------------------------------------
//     if (fSimp) {
//       p    = fSimp->fStartMom.P();
//       m    = 105.658; // in MeV
//       ekin = sqrt(p*p+m*m)-m;
//     }
//     Hist->fECalOverEKin->Fill(ecal/ekin);
//   }

//   Hist->fInstLumi->Fill(GetHeaderBlock()->InstLum());
// }

//-----------------------------------------------------------------------------
void TTrackAnaModule::FillCaloHistograms(CaloHist_t* Hist, TStnCrystal* Cr) {

  int                    nhits;
  float                  t, e, r, e700, n700;
  TCalHitData*           hit;

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


//-----------------------------------------------------------------------------
void TTrackAnaModule::FillClusterHistograms(ClusterHist_t* Hist, TStnCluster* Cluster) {
  int   row, col;
  float  x, y, z, r;

  row = Cluster->Ix1();
  col = Cluster->Ix2();

  x   = Cluster->fX+3904.;
  y   = Cluster->fY;
  z   = Cluster->fZ;
  r   = sqrt(x*x+y*y);

  if ((row < 0) || (row > 9999)) row = -9999;
  if ((col < 0) || (col > 9999)) col = -9999;

  Hist->fDiskID->Fill(Cluster->DiskID());
  Hist->fEnergy->Fill(Cluster->Energy());
  Hist->fT0->Fill(Cluster->Time());
  Hist->fRow->Fill(row);
  Hist->fCol->Fill(col);
  Hist->fX->Fill(x);
  Hist->fY->Fill(y);
  Hist->fZ->Fill(z);
  Hist->fR->Fill(r);

  Hist->fYMean->Fill(Cluster->fYMean);
  Hist->fZMean->Fill(Cluster->fZMean);
  Hist->fSigY->Fill(Cluster->fSigY);
  Hist->fSigZ->Fill(Cluster->fSigZ);
  Hist->fSigR->Fill(Cluster->fSigR);
  Hist->fNCr0->Fill(Cluster->fNCrystals);
  Hist->fNCr1->Fill(Cluster->fNCr1);
  Hist->fFrE1->Fill(Cluster->fFrE1);
  Hist->fFrE2->Fill(Cluster->fFrE2);
  Hist->fSigE1->Fill(Cluster->fSigE1);
  Hist->fSigE2->Fill(Cluster->fSigE2);
}



//-----------------------------------------------------------------------------
// fill efficiency histograms : need 10 histogram sets
// pitch = 1./tan(dip)
//-----------------------------------------------------------------------------
void TTrackAnaModule::FillEfficiencyHistograms(TStnTrackBlock*  TrackBlock, 
					       TStnTrackID*     TrackID   , 
					       int              HistSet   ) {
  // having MC truth is a must for calculating efficiency!
  if (fSimp == nullptr) return;

  if (fSimp->NStrawHits() >= 20) {
    FillEventHistograms(fHist.fEvent[HistSet],&fEvtPar);

    if (fSimp->fMomTrackerFront > 100.) {
      FillEventHistograms(fHist.fEvent[HistSet+1],&fEvtPar);

      TVector3 vdmom;

      if (fSimPar.fTFront != NULL) {
//-----------------------------------------------------------------------------
// regular case
//-----------------------------------------------------------------------------
	vdmom.SetXYZ(fSimPar.fTFront->Mom()->X(),
		     fSimPar.fTFront->Mom()->Y(),		      
		     fSimPar.fTFront->Mom()->Z());
      }
      else {
//-----------------------------------------------------------------------------
// pathological case when an upstream simulated MC particle is reconstructed 
// as the downstream one and there is no hit... efficiency doens't make sense
// just make sure the code doesnt' crash
//-----------------------------------------------------------------------------
	vdmom.SetXYZ(fSimPar.fParticle->fStartMom.X(),
		     fSimPar.fParticle->fStartMom.Y(),
		     fSimPar.fParticle->fStartMom.Z());
      }

      float ce_pitch  = vdmom.Pt()/vdmom.Pz();
      float min_pitch = 1./TrackID->MaxTanDip();
      float max_pitch = 1./TrackID->MinTanDip();

      if ((min_pitch < ce_pitch) && (ce_pitch < max_pitch)) {
	FillEventHistograms(fHist.fEvent[HistSet+2],&fEvtPar);
	  
	if (TrackBlock->NTracks() > 0) {
	  TStnTrack* track = TrackBlock->Track(0);
	  int id_word      = TrackID->IDWord(track);

	  FillEventHistograms(fHist.fEvent[HistSet+3],&fEvtPar);
	  
	  if ((id_word & TStnTrackID::kTrkQualBit) == 0) {
	    FillEventHistograms(fHist.fEvent[HistSet+4],&fEvtPar);
	    
	    if ((id_word & TStnTrackID::kT0Bit) == 0) {
	      FillEventHistograms(fHist.fEvent[HistSet+5],&fEvtPar);
	      
	      if ((id_word & TStnTrackID::kTanDipBit) == 0) {
		FillEventHistograms(fHist.fEvent[HistSet+6],&fEvtPar);
		
		if (((id_word & TStnTrackID::kD0Bit  ) == 0) && 
		    ((id_word & TStnTrackID::kRMaxBit) == 0)    ) {
		  
		  FillEventHistograms(fHist.fEvent[HistSet+7],&fEvtPar);
		  
		  if ((id_word & TStnTrackID::kTanDipBit) == 0) {
		    FillEventHistograms(fHist.fEvent[HistSet+8],&fEvtPar);

		    if ((103.5 < track->fP) && (track->fP < 105)) {
		      FillEventHistograms(fHist.fEvent[HistSet+9],&fEvtPar);
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
}

//_____________________________________________________________________________
void TTrackAnaModule::FillHistograms() {

  double       cos_th (-2.),  cl_e(-1.);
  int          disk_id(-1), alg_mask, best_alg;
  TStnCluster  *cl0(NULL);

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

  if ((fNTracks[0] == 0) && (fabs(cos_th) < 0.4)) {
    FillEventHistograms(fHist.fEvent[5],&fEvtPar); 
  }

  if (fNGoodTracks > 0) {
    FillEventHistograms(fHist.fEvent[6],&fEvtPar); 

    TLorentzVector    mom(1.,0.,0.,0);
    
    if (fParticle) fParticle->Momentum(mom);

    double p, cos_th;

    p      = mom.P();
    cos_th = mom.Pz()/p;

    if (GetDebugBit(31) && (cos_th > 0.8)) {
      GetHeaderBlock()->Print(Form(" bit:031 cos_th = %10.3f p = %10.3f ntrk = %5i",
				   cos_th, p, fNTracks[0]));
    }
  }

  if ((cl_e > fMinETrig) && (cl0->Time() > 500)) {
    FillEventHistograms(fHist.fEvent[7],&fEvtPar); 

    if (GetDebugBit(34)) {
      if (fNCalPatRec <= 0) {
	GetHeaderBlock()->Print(Form(" bit:034 cl_e = %10.3f, cl_time = %10.3f fNCalPatRec = 0",cl_e,cl0->Time()));
      }
    }

    if (GetDebugBit(38)) {
      if ((fNCalPatRec <= 0) && (fNTracks[0] > 0)) {
	GetHeaderBlock()->Print(Form(" bit:038 cl_e = %10.3f, cl_time = %10.3f fNCalPatRec = 0",cl_e,cl0->Time()));
      }
    }

    if (GetDebugBit(40)) {
      if (fNTracks[0] <= 0) {
	GetHeaderBlock()->Print(Form(" bit:040 cl_e = %10.3f, cl_time = %10.3f fNTracks[0] = %3i",
				     cl_e,cl0->Time(),fNTracks[0]));
      }
    }
  }

  if      (disk_id == 0) FillEventHistograms(fHist.fEvent[8],&fEvtPar);
  else if (disk_id == 1) FillEventHistograms(fHist.fEvent[9],&fEvtPar);
//-----------------------------------------------------------------------------
// EVT_10: events with SetC tracks 103.5 < p < 105.0
//-----------------------------------------------------------------------------
  if (fNGoodTracks > 0) {
    TStnTrack* trk = fTrackBlock->Track(0);
    if ((103.5 < trk->fP) && (trk->fP < 105.)) {
      FillEventHistograms(fHist.fEvent[10],&fEvtPar);
    }
  }
//-----------------------------------------------------------------------------
// Simp histograms
//-----------------------------------------------------------------------------
  if (fSimp) {
    FillSimpHistograms(fHist.fSimp[0],fSimp);
  }
//-----------------------------------------------------------------------------
// track histograms, fill them only for the downstream e- hypothesis
//-----------------------------------------------------------------------------
  TStnTrack*   trk;
  TrackPar_t*  tp;

  for (int i=0; i<fNTracks[0]; ++i ) {
    trk = fTrackBlock->Track(i);
    tp  = fTrackPar+i;

    FillTrackHistograms(fHist.fTrack[0],trk,tp,&fSimPar);

    if (tp->fIDWord[fBestID] == 0) {
					// track passes selection "C" 

      FillTrackHistograms(fHist.fTrack[1],trk,tp,&fSimPar);

      if (GetDebugBit(32) && (trk->fVMinS != NULL)) {
	if ((tp->fDtClZ0 > -0.6) && (trk->fVMinS->fEnergy > 50.)) { 
	  GetHeaderBlock()->Print(Form("bit032: odd DT: = %10.3f energy: %10.3f",tp->fDtClZ0, trk->fVMinS->fEnergy));
	}
      }

      if (GetDebugBit(35) && (trk->fP > 106.)) {
	GetHeaderBlock()->Print(Form("bit035: P = %10.3lf",trk->fP));
      }

					// events with at least one  cluster
      if (fNClusters > 0) {
	FillTrackHistograms(fHist.fTrack[2],trk,tp,&fSimPar);
      }
      else {
					// events without a cluster
	FillTrackHistograms(fHist.fTrack[3],trk,tp,&fSimPar);
      }
//-----------------------------------------------------------------------------
// events with a good track, reconstructed clusters but without a match
//-----------------------------------------------------------------------------
      if ((fNClusters > 0) && (trk->NClusters() == 0)) {
	FillTrackHistograms(fHist.fTrack[4],trk,tp,&fSimPar);
      }
//-----------------------------------------------------------------------------
// TRK 5 : events with a good track, reconstructed clusters but without a match
//-----------------------------------------------------------------------------
      if ((trk->fVMaxEp) && (fabs(trk->fVMaxEp->fDt) < 2.5)) {
	FillTrackHistograms(fHist.fTrack[5],trk,tp,&fSimPar);
      }
//-----------------------------------------------------------------------------
// TRK 8: good track, |xslope| < 3
//-----------------------------------------------------------------------------
      if (fabs(trk->XSlope()) < 3.) {
	FillTrackHistograms(fHist.fTrack[8],trk,tp,&fSimPar);
      }
    }
//-----------------------------------------------------------------------------
// TRK 6 : events with a track and a cluster E/P > 0.4
// TRK 7 : events with a "Set C" track and a cluster E/P > 0.4
//-----------------------------------------------------------------------------
    if ((trk->Ep() > 0.4) && ( trk->Ep() < 1.2)) {
      FillTrackHistograms(fHist.fTrack[6],trk,tp,&fSimPar);
      if (tp->fIDWord[fBestID] == 0) {
	FillTrackHistograms(fHist.fTrack[7],trk,tp,&fSimPar);
      }
    }
//----------------------------------------------------------------------------
//TRK  9: events with track and with no EM cluster      with E < 60 MeV
//TRK 10: events with good track and with no EM cluster with E < 60 MeV
//TRK 10: events with good track and with no EM cluster with E < 60 MeV
//----------------------------------------------------------------------------
    if (cl_e > 60) {
      FillTrackHistograms(fHist.fTrack[9],trk,tp,&fSimPar);
      if (tp->fIDWord[fBestID] == 0) {
	FillTrackHistograms(fHist.fTrack[10],trk,tp,&fSimPar);
      }
    }
//-----------------------------------------------------------------------------
// TRK 11: tracks with P > 103.5 MeV
//-----------------------------------------------------------------------------
    if (trk->P() > 103.5) FillTrackHistograms(fHist.fTrack[11],trk,tp,&fSimPar);
//-----------------------------------------------------------------------------
// TRK_12: tracks with fcon < 1e-4
// TRK_13: "Set C" tracks with 100 <= P < 110 
// TRK_14: tracks with fcon < 1e-2
//-----------------------------------------------------------------------------
    if (trk->fFitCons < 1.e-4) FillTrackHistograms(fHist.fTrack[12],trk,tp,&fSimPar);

    if ((tp->fIDWord[fBestID] == 0) && (trk->P() >= 100.) && (trk->P() < 110.)) {
      FillTrackHistograms(fHist.fTrack[13],trk,tp,&fSimPar);
    }

    if (trk->fFitCons < 1.e-2) FillTrackHistograms(fHist.fTrack[14],trk,tp,&fSimPar);

    TStnTrack::InterData_t*    vt = trk->fVMinS;  // track-only
//-----------------------------------------------------------------------------
// TRK_15: all tracks which have intersection with the 1st disk
// TRK_16: all tracks which have intersection with the 2nd disk
// TRK_17: all tracks which do not have intersections with the calorimeter
//-----------------------------------------------------------------------------
    if (vt) {
      if      (vt->fID == 0) FillTrackHistograms(fHist.fTrack[15],trk,tp,&fSimPar);
      else if (vt->fID == 1) FillTrackHistograms(fHist.fTrack[16],trk,tp,&fSimPar);
    }
    else {
      FillTrackHistograms(fHist.fTrack[17],trk,tp,&fSimPar);
    }
//-----------------------------------------------------------------------------
// TRK_18: BEST_ID tracks with T0 > 700
//-----------------------------------------------------------------------------
    if ((tp->fIDWord[fBestID] == 0) && (trk->T0() > 700)) {
      FillTrackHistograms(fHist.fTrack[18],trk,tp,&fSimPar);
    }
//-----------------------------------------------------------------------------
// TRK_19: BEST_ID tracks with an associated cluster
//-----------------------------------------------------------------------------
    if ((tp->fIDWord[fBestID] == 0) && (trk->Ep() > 0)) {
      FillTrackHistograms(fHist.fTrack[19],trk,tp,&fSimPar);

      if (trk->LogLHRCal() < 0) {
	if (GetDebugBit(29)) {
	  GetHeaderBlock()->Print(Form(" bit:029 LLHR(CAL) = %10.3f ep = %10.3f dt = %10.3f chi2_tcm = %10.3f",
				       trk->LogLHRCal(), trk->Ep(), trk->Dt(), trk->fVMinS->fChi2Match));
	}
      }

      if (tp->fDu < -80.) {
	if (GetDebugBit(33)) {
	  GetHeaderBlock()->Print(Form(" bit:033 DU = %10.3f dv = %10.3f ep = %10.3f dt = %10.3f",
				       tp->fDu, tp->fDv, tp->fEp, trk->Dt()));
	}
      }
    }
//-----------------------------------------------------------------------------
// TRK_20: BEST_ID tracks intersecting DISK=0
// TRK_21: BEST_ID tracks intersecting DISK=1
//-----------------------------------------------------------------------------
    if (vt && (tp->fIDWord[fBestID] == 0)) {
      if      (vt->fID == 0) FillTrackHistograms(fHist.fTrack[20],trk,tp,&fSimPar);
      else if (vt->fID == 1) FillTrackHistograms(fHist.fTrack[21],trk,tp,&fSimPar);
    }
//-----------------------------------------------------------------------------
// TRK 22: BEST_ID tracks with an associated cluster and chi2(match) < 100
// TRK 23: BEST_ID tracks with an associated cluster and chi2(match) < 100 and LLHR(cal) > 0
//         this is interesting to see which muons are getting misidentified
//-----------------------------------------------------------------------------
    if ((tp->fIDWord[fBestID] == 0) && (tp->fEp > 0) && (tp->fChi2Tcm < 100.)) {

      FillTrackHistograms(fHist.fTrack[22],trk,tp,&fSimPar);
      if (trk->LogLHRCal() > 0) {
	FillTrackHistograms(fHist.fTrack[23],trk,tp,&fSimPar);
	
	if (GetDebugBit(36) && (trk->fP < 80.)) {
	  GetHeaderBlock()->Print(Form(" bit:036 trk p = %10.3f E/P = %10.3f", trk->fP,tp->fEp));
	}
      }
      else {
//-----------------------------------------------------------------------------
// TRK_24: BEST_ID  tracks with an associated cluster and chi2(match) < 100 and LLHR(cal) < 0
//         this set allows to see which electrons are getting misidentified
//-----------------------------------------------------------------------------
	FillTrackHistograms(fHist.fTrack[24],trk,tp,&fSimPar);
      }
//-----------------------------------------------------------------------------
// TRK_25: BEST_ID tracks, 100 < P < 110, 0 < E/p < 1.15,  |dt_corr| < 3, chi2(match) < 100
//-----------------------------------------------------------------------------
//      double dt_corr = trk->Dt(); // -1.;
      if ( (tp->fDt > fMinDtTcm) && (tp->fDt < fMaxDtTcm) && 
	   (tp->fEp <      1.15) && 
	   (trk->fP >      100.) && (trk->fP <       110)   ) {
	FillTrackHistograms(fHist.fTrack[25],trk,tp,&fSimPar);
//-----------------------------------------------------------------------------
// more details on the calorimeter-based likelihood 
// TRK_26 : TRK 25 events with LLHR_CAL > 0 ( interesting for muons)
// TRK_27 : TRK 25 events with LLHR_CAL < 0 ( interesting for electrons)
//-----------------------------------------------------------------------------
	double llhr_cal = trk->LogLHRCal();
	if (llhr_cal > 0) {
	  FillTrackHistograms(fHist.fTrack[26],trk,tp,&fSimPar);

	  if (llhr_cal > 5.) {
	    if (GetDebugBit(37)) {
	      GetHeaderBlock()->Print(Form("bit:037 trk p = %10.3f E/P = %10.3f", trk->fP,tp->fEp));
	    }
	  }
	}
	else {
	  FillTrackHistograms(fHist.fTrack[27],trk,tp,&fSimPar);
	}
      }

    }
//-----------------------------------------------------------------------------
// TRK_28 : events with a BEST_ID track and a cluster E/P > 1.1
//-----------------------------------------------------------------------------
    if (tp->fIDWord[fBestID] == 0) {
      if (tp->fEp > 1.1) {
	FillTrackHistograms(fHist.fTrack[28],trk,tp,&fSimPar);

	if (GetDebugBit(28)) {
	  GetHeaderBlock()->Print(Form(" bit:028 LLHR(CAL) = %10.3f ep = %10.3f dt = %10.3f",
				       trk->LogLHRCal(), trk->Ep(), trk->Dt()));
	}
      }
    }
//-----------------------------------------------------------------------------
// TRK_29: BEST_ID track 100 < P < 110, 0 < E/P < 1.15              : next to TRK_13
// TRK_32: BEST_ID track 100 < P < 110, 0 < E/P < 1.15, chi2tcm<100 : next to TRK_29
//-----------------------------------------------------------------------------
    if (tp->fIDWord[fBestID] == 0) {
      if ((tp->fEp > 0) && (trk->P() > 100.) && (trk->P() < 110.) && (tp->fEp < 1.15)) {
	FillTrackHistograms(fHist.fTrack[29],trk,tp,&fSimPar);

	if (tp->fChi2Tcm < 100.) {
	  FillTrackHistograms(fHist.fTrack[32],trk,tp,&fSimPar);
	}
      }
    }
//-----------------------------------------------------------------------------
// TRK_30: tracks with >= 25 hits
// TRK_31: tracks with >= 25 hits and Chi2/Ndof < 3
//-----------------------------------------------------------------------------
    if (trk->NActive() >= 25) {
      FillTrackHistograms(fHist.fTrack[30],trk,tp,&fSimPar);
      if (trk->Chi2Dof() < 3) {
	FillTrackHistograms(fHist.fTrack[31],trk,tp,&fSimPar);
      }
    }
//-----------------------------------------------------------------------------
// TRK_33: "DaveTrkQual" Q > 0.4 tracks
//-----------------------------------------------------------------------------
    if (tp->fIDWord[fBestID] == 0) {
      FillTrackHistograms(fHist.fTrack[33],trk,tp,&fSimPar);
    }
//-----------------------------------------------------------------------------
// TRK_34: "DaveTrkQual" Q > 0.1 tracks
//-----------------------------------------------------------------------------
    if (tp->fIDWord[1] == 0) {
      FillTrackHistograms(fHist.fTrack[34],trk,tp,&fSimPar);
    }
//-----------------------------------------------------------------------------
// split tracks by the algorithm mask: 1 , 2 , or 3
//-----------------------------------------------------------------------------
    alg_mask = trk->AlgMask();
    best_alg = trk->BestAlg();
    if      (alg_mask == 1) {
//-----------------------------------------------------------------------------
// TrkPatRec-only tracks
//-----------------------------------------------------------------------------
      FillTrackHistograms(fHist.fTrack[40],trk,tp,&fSimPar);
      if (tp->fIDWord[fBestID] == 0) {
	FillTrackHistograms(fHist.fTrack[41],trk,tp,&fSimPar);
	// print run/event numbers :
	if (GetDebugBit(6)) {
	  double ep = trk->Ep();
	  if ((ep > 0.8) && (ep < 1.1)) {
	    GetHeaderBlock()->Print(Form(" bit:006 trk_41: track E/P = %8.3f",trk->Ep()));
	  }
	}
	if (trk->T0() > 700.) FillTrackHistograms(fHist.fTrack[42],trk,tp,&fSimPar);
      }
    }
    else if (alg_mask == 2) {
//-----------------------------------------------------------------------------
// CalPatRec-only tracks
//-----------------------------------------------------------------------------
      FillTrackHistograms(fHist.fTrack[50],trk,tp,&fSimPar);
      if (tp->fIDWord[fBestID] == 0) {
	FillTrackHistograms(fHist.fTrack[51],trk,tp,&fSimPar);
	if (trk->T0() > 700.) FillTrackHistograms(fHist.fTrack[52],trk,tp,&fSimPar);
      }
    }
    else if (alg_mask == 3) {
//-----------------------------------------------------------------------------
// TrkPatRec+CalPatRec tracks
//-----------------------------------------------------------------------------
      FillTrackHistograms(fHist.fTrack[60],trk,tp,&fSimPar);
      if (tp->fIDWord[fBestID] == 0) {
	FillTrackHistograms(fHist.fTrack[61],trk,tp,&fSimPar);
	if (trk->T0() > 700.) FillTrackHistograms(fHist.fTrack[62],trk,tp,&fSimPar);
      }
    }
//-----------------------------------------------------------------------------
// TRK_63 : fIDWORD[bestID] == 0, TRkPatRec
// TRK_64 : fIDWORD[bestID] == 0, CalPatRec
//-----------------------------------------------------------------------------
    if (tp->fIDWord[fBestID] == 0) {
      if      (best_alg == 0) FillTrackHistograms(fHist.fTrack[63],trk,tp,&fSimPar);
      else if (best_alg == 1) FillTrackHistograms(fHist.fTrack[64],trk,tp,&fSimPar);

      if (tp->fDtClZ0 > -0.6) FillTrackHistograms(fHist.fTrack[65],trk,tp,&fSimPar);
    }

    
//-----------------------------------------------------------------------------
// TRK_71: SetC tracks  103.5 < p < 105 : DIO studies
//-----------------------------------------------------------------------------
    if ((tp->fIDWord[fBestID] == 0) && (trk->fP > 103.5) && (trk->fP < 105.)) {
      FillTrackHistograms(fHist.fTrack[71],trk,tp,&fSimPar);
    }
  }

//-----------------------------------------------------------------------------
// Track Reco efficiency, Dave style
//-----------------------------------------------------------------------------
  FillEfficiencyHistograms(fTrackBlock,TAnaModule::fTrackID[fBestID],11);
//-----------------------------------------------------------------------------
// cluster histograms 
//-----------------------------------------------------------------------------
  TStnCluster*  cl;
  int           id;
  for (int i=0; i<fNClusters; ++i ) {
    cl = fClusterBlock->Cluster(i);
    id = cl->DiskID();
    FillClusterHistograms(fHist.fCluster[0],cl);

    if (fNTracks[0]     >  0 ) FillClusterHistograms(fHist.fCluster[1],cl);
    if (fNGoodTracks    >  0 ) FillClusterHistograms(fHist.fCluster[2],cl);
    if (fNMatchedTracks >  0 ) FillClusterHistograms(fHist.fCluster[3],cl);
    if (cl->Energy()    > 10.) FillClusterHistograms(fHist.fCluster[4],cl);
    if (cl->Energy()    > 60.) FillClusterHistograms(fHist.fCluster[5],cl);

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

	if (cr->Energy() > 0) {
	  FillCaloHistograms(fHist.fCalo[1],cr);
	}
	if (cr->Energy() > 0.1) {
	  FillCaloHistograms(fHist.fCalo[2],cr);
	}
	if (cr->Energy() > 1.0) {
	  FillCaloHistograms(fHist.fCalo[3],cr);
	}
      }
    }
  }
//-----------------------------------------------------------------------------
// radial distributions for crystals
//-----------------------------------------------------------------------------
  static int first_entry(1);

  if (first_entry == 1) {
    first_entry = 0;

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
}



//-----------------------------------------------------------------------------
// 2014-04-30: it looks that reading the straw hits takes a lot of time - 
//              turn off by default by commenting it out
//-----------------------------------------------------------------------------
int TTrackAnaModule::Event(int ientry) {

  TDiskCalorimeter::GeomData_t disk_geom;

  fTrackBlock->GetEntry(ientry);
  fTrackStrawHitBlock->GetEntry(ientry);
  fClusterBlock->GetEntry(ientry);
  fStrawDataBlock->GetEntry(ientry);
  fCalDataBlock->GetEntry(ientry);
  fGenpBlock->GetEntry(ientry);
  fSimpBlock->GetEntry(ientry);
  fVDetBlock->GetEntry(ientry);
//-----------------------------------------------------------------------------
// at some point, the TVdetBlock class got obsolete, and now the virtual detector 
// hits are stored in TStepPointMCBlock 
// dont' try to read it as that would fail
//-----------------------------------------------------------------------------
  //  fVDetBlock->GetEntry(ientry);
  fEvtPar.fNCrvClusters     = -1;
  fEvtPar.fNCrvPulses       = -1;
  fEvtPar.fNCrvCoincidences = -1;
  fEvtPar.fNGenp            = fGenpBlock->NParticles();
//-----------------------------------------------------------------------------
// luminosity weight
//-----------------------------------------------------------------------------
  fLumWt = GetHeaderBlock()->LumWeight();
//-----------------------------------------------------------------------------
// look for signal particle defined by the PDG code and the generator code
//-----------------------------------------------------------------------------
  TGenParticle* genp;
  int           pdg_code, generator_code;
  
  fEvtPar.fParticle = NULL;
  for (int i=fNGenp-1; i>=0; i--) {
    genp           = fGenpBlock->Particle(i);
    pdg_code       = genp->GetPdgCode();
    generator_code = genp->GetStatusCode();
    if ((abs(pdg_code) == fPdgCode) && (generator_code == fGeneratorCode)) {
      fEvtPar.fParticle = genp;
      break;
    }
  }
//-----------------------------------------------------------------------------
// may want to revisit the definition of fSimp in the future
//-----------------------------------------------------------------------------
  fSimp             = fSimpBlock->Particle(0);
  fSimPar.fParticle = fSimp;
  fSimPar.fTFront   = NULL;
  fSimPar.fTMid     = NULL;
  fSimPar.fTBack    = NULL;
  fSimPar.fGenp     = fEvtPar.fParticle;
//-----------------------------------------------------------------------------
// process virtual detectors - for fSimp need parameters at tracker entrance
// use the first fit
//-----------------------------------------------------------------------------
  fSimPar.fTFront   = NULL;
  fSimPar.fTMid     = NULL;
  fSimPar.fTBack    = NULL;
  int nvdhits = fVDetBlock->NStepPoints();
  for (int i=0; i<nvdhits; i++) {
    TStepPointMC* vdhit = fVDetBlock->StepPointMC(i);
    if (vdhit->PDGCode() == fSimp->fPdgCode) {
      if ((vdhit->VolumeID() == 13) || (vdhit->VolumeID() == 14)) {
	if (fDirection*vdhit->Mom()->Z() > 0) {
	  if (fSimPar.fTFront == 0) fSimPar.fTFront = vdhit;
	}
      }
      else if ((vdhit->VolumeID() == 11) || (vdhit->VolumeID() == 12)) {
	if (fDirection*vdhit->Mom()->Z() > 0) {
	  if (fSimPar.fTMid == 0) fSimPar.fTMid = vdhit;
	}
      }
      else if (vdhit->VolumeID() == mu2e::VirtualDetectorId::TT_Back) {
	if (fDirection*vdhit->Mom()->Z() > 0) {
	  if (fSimPar.fTBack == NULL) fSimPar.fTBack = vdhit;
	}
      }
    }
  }

  if (fDiskCalorimeter->Initialized() == 0) {
    disk_geom.fNDisks = fCalDataBlock->NDisks();

    for (int i=0; i<disk_geom.fNDisks; i++) {
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

  fNTracks[0] = fTrackBlock->NTracks();
  if (fNTracks[0] == 0) fTrack = 0;
  else                  fTrack = fTrackBlock->Track(0);

  fNClusters  = fClusterBlock->NClusters();
  fNCalHits   = fCalDataBlock->NHits();
  fNStrawHits = GetHeaderBlock()->fNStrawHits;

  fDiskCalorimeter->InitEvent(fCalDataBlock);
//-----------------------------------------------------------------------------
// loop over tracks and calculate needed parameters
//-----------------------------------------------------------------------------
  fNHyp           = -1;
  fBestHyp[0]     = -1;
  fBestHyp[1]     = -1;

  fNGoodTracks    = 0;
  fNMatchedTracks = 0;
//-----------------------------------------------------------------------------
// determine the number of CalPatRec tracks - this assumes that the list of 
// tracks has been created by MergePatRec
//-----------------------------------------------------------------------------
  int ntrk = fTrackBlock->NTracks();
  fNCalPatRec = 0;
  for (int itrk=0; itrk<ntrk; itrk++) {
    TStnTrack*   trk = fTrackBlock->Track(itrk);
    TrackPar_t* tp   = fTrackPar+itrk;

    int alg_mask = trk->AlgMask();
    if (alg_mask & 0x2) fNCalPatRec += 1;

    tp->fTrackID[0] = TAnaModule::fTrackID_BOX;
    tp->fTrackID[1] = TAnaModule::fTrackID_MVA;
    //    tp->fLogLH      = TAnaModule::fLogLH;
  }

  InitTrackPar(fTrackBlock,fClusterBlock,fTrackPar,&fSimPar);
//-----------------------------------------------------------------------------
// init calorimeter clusters - remember, the first one not necessarily is the 
// most energetic
//-----------------------------------------------------------------------------
  fNClusters = fClusterBlock->NClusters();
  if (fNClusters == 0) fCluster = 0;
  else                 fCluster = fClusterBlock->Cluster(0);

  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TTrackAnaModule::Debug() {

  TStnTrack* trk;
  TrackPar_t* tp;
  int ntrk = fTrackBlock->NTracks();

  for (int itrk=0; itrk<ntrk; itrk++) {
    trk = fTrackBlock->Track(itrk);
    tp  = &fTrackPar[itrk];
//-----------------------------------------------------------------------------
// bit 0: all tracks
//-----------------------------------------------------------------------------
    if (GetDebugBit(0) == 1) {
	GetHeaderBlock()->Print(Form("bit_000: All p = %10.3f",
				     trk->Momentum()->P()));
    }
//-----------------------------------------------------------------------------
// bit 1: IDWord =0 0 tracks
//-----------------------------------------------------------------------------
    if (GetDebugBit(1) == 1) {
      if (tp->fIDWord[fBestID] == 0) {
	GetHeaderBlock()->Print(Form("bit_001: IDWord=0 p = %10.3f",
				     trk->Momentum()->P()));
      }
    }
//-----------------------------------------------------------------------------
// bit 3: Set C tracks with large DX : 70mm < |DX| < 90mm
//-----------------------------------------------------------------------------
    if (GetDebugBit(3) == 1) {
      if (tp->fIDWord[fBestID] == 0) {
	TStnTrack::InterData_t*    vr = trk->fVMaxEp; // residuals
	if ((vr && (fabs(vr->fDx) > 70) && (fabs(vr->fDx) < 90))) {
	  GetHeaderBlock()->Print(Form("large DX: %f",vr->fDx));
	}
      }
    }
//-----------------------------------------------------------------------------
// bit 4: tracks with DpF > 1MeV - positive tail...
//-----------------------------------------------------------------------------
    if (GetDebugBit(4) == 1) {
      if (tp->fDpF > 1.) {
	GetHeaderBlock()->Print(Form("pF pRec, fDpf = %10.3f  %10.3f  %10.3f",
				     trk->fPFront, trk->Momentum()->P(),tp->fDpF));
      }
    }
//-----------------------------------------------------------------------------
// bit 9: Set C tracks with DpF > 1MeV - positive tail...
//-----------------------------------------------------------------------------
    if (GetDebugBit(9) == 1) {
      double ep = trk->Ep();
      if (tp->fIDWord[fBestID] == 0) { 
	if (((ep > 0.42) && (ep < 0.46)) || ((ep > 0.35) && (ep < 0.39))) {
	  GetHeaderBlock()->Print(Form("bit:009 ep = %10.3f e = %10.3f p = %10.3f",
				       trk->fEp,trk->fEp*trk->fP,trk->fP));
	}
      }
    }
//-----------------------------------------------------------------------------
// bit 10: Set C tracks with Ecl > 80
//-----------------------------------------------------------------------------
    if (GetDebugBit(10) == 1) {
      double ecl = trk->ClusterE();
      if (tp->fIDWord[fBestID] == 0) { 
	if (ecl > 60) {
	  GetHeaderBlock()->Print(Form("bit:010 e = %10.3f p = %10.3f",
				       ecl,trk->fP));
	}
      }
    }
  }

//-----------------------------------------------------------------------------
// bit 5: events with N(tracks) > 1
//-----------------------------------------------------------------------------
  if (GetDebugBit(5) == 1) {
    int ntrk = fTrackBlock->NTracks();
    if (ntrk > 1) {
      GetHeaderBlock()->Print(Form("NTracks = %5i",ntrk));
    }
  }
}

//_____________________________________________________________________________
int TTrackAnaModule::EndJob() {
  printf("----- end job: ---- %s\n",GetName());
  return 0;
}

//_____________________________________________________________________________
void TTrackAnaModule::Test001() {
}

}
