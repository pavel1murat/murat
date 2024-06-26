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
//  3  : events with good tracks P > 70 and T0 > 300
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
#include "Stntuple/geom/TDisk.hh"
#include "Stntuple/obj/TStnNode.hh"
#include "Stntuple/val/stntuple_val_functions.hh"
#include "Offline/DataProducts/inc/VirtualDetectorId.hh"

//------------------------------------------------------------------------------
// Mu2e offline includes
//-----------------------------------------------------------------------------

#include "murat/ana/TPipenuAnaModule.hh"

ClassImp(murat::TPipenuAnaModule)

namespace murat {
//-----------------------------------------------------------------------------
TPipenuAnaModule::TPipenuAnaModule(const char* name, const char* title): 
  TAnaModule(name,title) { 
  fTrackBlockName         = "TrackBlock";
  fTrackStrawHitBlockName = "TrackStrawHitBlock";
  fTrackNumber.Set(100);

  //  fDiskCalorimeter = new TDiskCalorimeter();
  fCalorimeterType = 2;

  fMinT0           = 700; 
  //  fMinETrig        = 50.;               // MeV
					// track-cluster matching timing cut
  // fMinDtTcm        = -5.;
  // fMaxDtTcm        =  8.;
//-----------------------------------------------------------------------------
// initialize Track ID
// 0: SetC  1: TrkQual>0.1 2:TrkQual>0.4
// what about number of hits ? - 3: no cuts on the number of hits
//-----------------------------------------------------------------------------
  fNID             = 2;

  fTrackID[0]      = new TStnTrackID();

  int mask0 = TStnTrackID::kNActiveBit | TStnTrackID::kChi2DofBit | TStnTrackID::kMomErrBit ;
  fTrackID[0]->SetMaxChi2Dof(3.0 );
  fTrackID[0]->SetMinNActive(25  );
  fTrackID[0]->SetMaxMomErr (0.25);

  mask0    |= TStnTrackID::kTanDipBit;
  fTrackID[0]->SetMinTanDip (0.6);
  fTrackID[0]->SetMaxTanDip (1.0);

  mask0    |= TStnTrackID::kD0Bit;
  fTrackID[0]->SetMinD0 (-100);         // mm
  fTrackID[0]->SetMaxD0 ( 100);

  fTrackID[0]->SetUseMask(mask0);
//-----------------------------------------------------------------------------
// tight ID
//-----------------------------------------------------------------------------
  fTrackID[1]      = new TStnTrackID();

  int mask1 = TStnTrackID::kNActiveBit | TStnTrackID::kChi2DofBit | TStnTrackID::kMomErrBit ;
  fTrackID[1]->SetMaxChi2Dof(3.0 );
  fTrackID[1]->SetMinNActive(30  );     // 30 ... cut tighter on the number of hits
  fTrackID[1]->SetMaxMomErr (0.25);

  mask1    |= TStnTrackID::kTanDipBit;
  fTrackID[1]->SetMinTanDip (0.70);     // prev: 0.65 
  fTrackID[1]->SetMaxTanDip (1.20);

  mask1    |= TStnTrackID::kD0Bit;
  fTrackID[1]->SetMinD0 ( -50);         // mm
  fTrackID[1]->SetMaxD0 ( 100);

  fTrackID[1]->SetUseMask(mask1);

  fBestID           = 0;

  for (int i=0; i<20; i++) {
    TrackPar_t* tp  = fTrackPar+i;
    tp->fTrackID[0] = TAnaModule::fTrackID[0];
    tp->fTrackID[1] = TAnaModule::fTrackID[1];
  }

//-----------------------------------------------------------------------------
// MC truth: define which MC particle to consider as signal
// 2:conversionGun, 28:StoppedParticleReactionGun - see 
//-----------------------------------------------------------------------------
  fDirection        = 1;

  fDnMax            = 15;
//-----------------------------------------------------------------------------
// this is redefiend by the dataset catalog anyway
//-----------------------------------------------------------------------------
  fPDGCode          = -11;
  fMCProcessCode    = 181;   // mu2e::ProcessCode::mu2ePienu;
}

//-----------------------------------------------------------------------------
TPipenuAnaModule::~TPipenuAnaModule() {
  for (int i=0; i<fNID; i++) delete fTrackID[i];
}


//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TPipenuAnaModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
  RegisterDataBlock(fTrackBlockName.Data()        ,"TStnTrackBlock"     ,&fTrackBlock  );
  // RegisterDataBlock(fTrackStrawHitBlockName.Data(),"TTrackStrawHitBlock",&fTrackStrawHitBlock);

  RegisterDataBlock("KSFBlock"        ,"TStnTrackSeedBlock"  ,&fTrackSeedBlock  );
  RegisterDataBlock("HelixBlock"      ,"TStnHelixBlock"      ,&fHelixBlock      );
  RegisterDataBlock("TimeClusterBlock","TStnTimeClusterBlock",&fTimeClusterBlock);

  RegisterDataBlock("ClusterBlock" ,"TStnClusterBlock" ,&fClusterBlock );
  RegisterDataBlock("CalDataBlock" ,"TCalDataBlock"    ,&fCalDataBlock );
  // RegisterDataBlock("StrawHitBlock","TStrawHitBlock"   ,&fStrawHitBlock);
  RegisterDataBlock("GenpBlock"    ,"TGenpBlock"       ,&fGenpBlock    );
  RegisterDataBlock("SimpBlock"    ,"TSimpBlock"       ,&fSimpBlock    );
  RegisterDataBlock("SpmcBlockVDet","TStepPointMCBlock",&fVDetBlock    );
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();

  return 0;
}

//-----------------------------------------------------------------------------
void TPipenuAnaModule::BookPipenuHistograms(PipenuHist_t* Hist, const char* Folder) {
  HBook2F(Hist->fPVsT0 ,"p_vs_t0",Form("%s: track P s T0"       ,Folder), 200, 0,  2000, 100, 0, 100,Folder);
}

//-----------------------------------------------------------------------------
void TPipenuAnaModule::BookTimeClusterHistograms(TimeClusterHist_t* Hist, const char* Folder) {
  HBook1F(Hist->fNsh  ,"nsh"  ,Form("%s: N(single straw hits)",Folder), 200, 0,  200,Folder);
  HBook1F(Hist->fNch  ,"nch"  ,Form("%s: N(combo hits)"       ,Folder), 200, 0,  200,Folder);
  HBook1F(Hist->fT0   ,"t0"   ,Form("%s: T0"                  ,Folder), 200, 0, 2000,Folder);
  HBook1F(Hist->fT0Err,"t0err",Form("%s: T0Err"               ,Folder), 200, 0, 2000,Folder);
}

//_____________________________________________________________________________
void TPipenuAnaModule::BookHistograms() {

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
  // HBook1F(fHist.fCrystalR[0],"rc_0"     ,Form("disk [0] crystal radius"),100,0,1000,"Hist");
  // HBook1F(fHist.fCrystalR[1],"rc_1"     ,Form("disk [1] crystal radius"),100,0,1000,"Hist");
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
// book helix histograms
//-----------------------------------------------------------------------------
  int book_helix_histset[kNHelixHistSets];
  for (int i=0; i<kNHelixHistSets; i++) book_helix_histset[i] = 0;

  book_helix_histset[ 0] = 1;		// all events
  book_helix_histset[ 1] = 1;		// events with N(reconstructed tracks) > 0
  book_helix_histset[ 2] = 1;		// helices with a TrackSeed

  for (int i=0; i<kNHelixHistSets; i++) {
    if (book_helix_histset[i] != 0) {
      sprintf(folder_name,"hel_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fHelix[i] = new HelixHist_t;
      BookHelixHistograms(fHist.fHelix[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book simp histograms
//-----------------------------------------------------------------------------
  int book_simp_histset[kNSimpHistSets];
  for (int i=0; i<kNSimpHistSets; i++) book_simp_histset[i] = 0;

  book_simp_histset[ 0] = 1;		// all events
  book_simp_histset[ 1] = 1;		// events with a      track
  book_simp_histset[ 2] = 1;		// events with a good track

  book_simp_histset[10] = 1;		// all events                     volID=3015 (last foil)
  book_simp_histset[11] = 1;		// events with a reco track and a volID=3015 (last foil)

  book_simp_histset[31] = 1;		// muons events with a      track T>300
  book_simp_histset[32] = 1;		// muons events with a      track T>300
  book_simp_histset[33] = 1;		// muons events with a      track T>300

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

  book_track_histset[  0] = 1;		// all tracks

  book_track_histset[  1] = 1;		// fIDWord[0]=0 tracks
  book_track_histset[ 21] = 1;		// fIDWord[1]=0 tracks

  book_track_histset[ 50] = 1;          // all         tracks wt=pionSurvProb

  book_track_histset[ 51] = 1;          // IDWord[0]=0 tracks wt=pionSurvProb
  book_track_histset[ 71] = 1;          // IDWord[1]=0 tracks wt=pionSurvProb
  
  book_track_histset[100] = 1;		// all  tracks Q=1

  book_track_histset[101] = 1;		// IDWord[0]=0 tracks Q=1
  book_track_histset[102] = 1;		// IDWord[0]=0 tracks Q=1 p>60 
  book_track_histset[103] = 1;		// IDWord[0]=0 tracks Q=1 p>60 T>300 
  book_track_histset[104] = 1;		// IDWord[0]=0 tracks Q=1 p>60 T>300 DnTrackTc 
  book_track_histset[105] = 1;		// IDWord[0]=0 tracks Q=1 p>60 T>300 DnTrackTc DtTcTc
  book_track_histset[106] = 1;		// IDWord[0]=0 tracks Q=1 T>300 
  book_track_histset[107] = 1;		// IDWord[0]=0 tracks Q=1 T>300 DnTrackTc 
  book_track_histset[108] = 1;		// IDWord[0]=0 tracks Q=1 T>300 DnTrackTc DtTcTc

  book_track_histset[121] = 1;		// IDWord[1]=0 tracks Q=1
  book_track_histset[122] = 1;		// IDWord[1]=0 tracks Q=1 p>60 
  book_track_histset[123] = 1;		// IDWord[1]=0 tracks Q=1 p>60 T>300 
  book_track_histset[124] = 1;		// IDWord[1]=0 tracks Q=1 p>60 T>300 DnTrackTc 
  book_track_histset[125] = 1;		// IDWord[1]=0 tracks Q=1 p>60 T>300 DnTrackTc DtTcTc
  book_track_histset[126] = 1;		// IDWord[1]=0 tracks Q=1 T>300 
  book_track_histset[127] = 1;		// IDWord[1]=0 tracks Q=1 T>300 DnTrackTc 
  book_track_histset[128] = 1;		// IDWord[1]=0 tracks Q=1 T>300 DnTrackTc DtTcTc

  book_track_histset[150] = 1;          // all         tracks          wt=pionSurvProb Q=1

  book_track_histset[151] = 1;          // IDWord[0]=0 tracks, wt=pionSurvProb Q=1
  book_track_histset[152] = 1;          // IDWord[0]=0 tracks, wt=pionSurvProb Q=1 p>60 
  book_track_histset[153] = 1;          // IDWord[0]=0 tracks, wt=pionSurvProb Q=1 p>60 T>300
  book_track_histset[154] = 1;          // IDWord[0]=0 tracks, wt=pionSurvProb Q=1 p>60 T>300 DnTrackTC 
  book_track_histset[155] = 1;          // IDWord[0]=0 tracks, wt=pionSurvProb Q=1 p>60 T>300 DnTrackTC DtTcTc
  book_track_histset[156] = 1;          // IDWord[0]=0 tracks, wt=pionSurvProb Q=1 T>300
  book_track_histset[157] = 1;          // IDWord[0]=0 tracks, wt=pionSurvProb Q=1 T>300 DnTrackTC 
  book_track_histset[158] = 1;          // IDWord[0]=0 tracks, wt=pionSurvProb Q=1 T>300 DnTrackTC DtTcTc
  
  book_track_histset[171] = 1;          // IDWord[1]=0 tracks, wt=pionSurvProb Q=1
  book_track_histset[172] = 1;          // IDWord[1]=0 tracks, wt=pionSurvProb Q=1 p>60 
  book_track_histset[173] = 1;          // IDWord[1]=0 tracks, wt=pionSurvProb Q=1 p>60 T>300
  book_track_histset[174] = 1;          // IDWord[1]=0 tracks, wt=pionSurvProb Q=1 p>60 T>300 DnTrackTC 
  book_track_histset[175] = 1;          // IDWord[1]=0 tracks, wt=pionSurvProb Q=1 p>60 T>300 DnTrackTC DtTcTc
  book_track_histset[176] = 1;          // IDWord[1]=0 tracks, wt=pionSurvProb Q=1 T>300
  book_track_histset[177] = 1;          // IDWord[1]=0 tracks, wt=pionSurvProb Q=1 T>300 DnTrackTC 
  book_track_histset[178] = 1;          // IDWord[1]=0 tracks, wt=pionSurvProb Q=1 T>300 DnTrackTC DtTcTc
  
  book_track_histset[200] = 1;		// all  tracks Q=-1
  book_track_histset[201] = 1;		// IDWord[0]=0 tracks Q=-1

  book_track_histset[221] = 1;		// IDWord[1]=0 tracks Q=-1

  book_track_histset[250] = 1;          // all         tracks wt=pionSurvPro Q=-1
  book_track_histset[251] = 1;          // IDWord[0]=0 tracks wt=pionSurvPro Q=-1

  book_track_histset[271] = 1;		// IDWord[1]=0 tracks wt=pionSurvPro Q=-1

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
// book track histograms
//-----------------------------------------------------------------------------
  int book_track_id_histset[kNTrackIDHistSets];
  for (int i=0; i<kNTrackIDHistSets; i++) book_track_id_histset[i] = 0;

  book_track_id_histset[  0] = 1;		// all tracks

  for (int i=0; i<kNTrackIDHistSets; i++) {
    if (book_track_id_histset[i] == 0)                      continue;
    sprintf(folder_name,"tid_%i",i);
    fol = (TFolder*) hist_folder->FindObject(folder_name);
    if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
    fHist.fTrackID[i] = new TStnTrackID::Hist_t;
    BookTrackIDHistograms(fHist.fTrackID[i],Form("Hist/%s",folder_name));
  }
//-----------------------------------------------------------------------------
// book time cluster histograms
//-----------------------------------------------------------------------------
  int book_tc_histset[kNTimeClusterHistSets];
  for (int i=0; i<kNTimeClusterHistSets; i++) book_tc_histset[i] = 0;

  book_tc_histset[0] = 1;		// all clusters

  for (int i=0; i<kNTimeClusterHistSets; i++) {
    if (book_tc_histset[i] == 0)                            continue;
    sprintf(folder_name,"tc_%i",i);
    fol = (TFolder*) hist_folder->FindObject(folder_name);
    if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
    fHist.fTimeCluster[i] = new TimeClusterHist_t;
    BookTimeClusterHistograms(fHist.fTimeCluster[i],Form("Hist/%s",folder_name));
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
//-----------------------------------------------------------------------------
// book pipenu histograms
//-----------------------------------------------------------------------------
  int book_pipenu_histset[kNPipenuHistSets];
  for (int i=0; i<kNPipenuHistSets; i++) book_pipenu_histset[i] = 0;

  book_pipenu_histset[  0] = 1;		// all  tracks

  book_pipenu_histset[  1] = 1;		// IDWord[0]=0 tracks
  book_pipenu_histset[ 21] = 1;		// IDWord[1]=0 tracks

  book_pipenu_histset[ 50] = 1;		// all  tracks, wt = pion surv prob

  book_pipenu_histset[ 51] = 1;		// IDWord[0]=0 tracks, wt=pionSurvProb
  book_pipenu_histset[ 71] = 1;		// IDWord[1]=0 tracks, wt=pionSurvProb

  book_pipenu_histset[100] = 1;		// all  tracks Q>0

  book_pipenu_histset[101] = 1;		// IDWord[0]=0 tracks Q>0
  book_pipenu_histset[102] = 1;		// IDWord[0]=0 tracks Q>0 p>60
  book_pipenu_histset[103] = 1;		// IDWord[0]=0 tracks Q>0 p>60 T>300
  book_pipenu_histset[104] = 1;		// IDWord[0]=0 tracks Q>0 p>60 T>300 DnTrackTc
  book_pipenu_histset[105] = 1;		// IDWord[0]=0 tracks Q>0 p>60 T>300 DnTrackTc DtTcTc
  book_pipenu_histset[106] = 1;		// IDWord[0]=0 tracks Q>0 T>300
  book_pipenu_histset[107] = 1;		// IDWord[0]=0 tracks Q>0 T>300 DnTrackTc
  book_pipenu_histset[108] = 1;		// IDWord[0]=0 tracks Q>0 T>300 DnTrackTc DtTcTc

  book_pipenu_histset[121] = 1;		// IDWord[1]=0 tracks Q>0
  book_pipenu_histset[122] = 1;		// IDWord[1]=0 tracks Q>0 p>60
  book_pipenu_histset[123] = 1;		// IDWord[1]=0 tracks Q>0 p>60 T>300
  book_pipenu_histset[124] = 1;		// IDWord[1]=0 tracks Q>0 p>60 T>300 DnTrackTc
  book_pipenu_histset[125] = 1;		// IDWord[1]=0 tracks Q>0 p>60 T>300 DnTrackTc DtTcTc
  book_pipenu_histset[126] = 1;		// IDWord[1]=0 tracks Q>0 T>300
  book_pipenu_histset[127] = 1;		// IDWord[1]=0 tracks Q>0 T>300 DnTrackTc
  book_pipenu_histset[128] = 1;		// IDWord[1]=0 tracks Q>0 T>300 DnTrackTc DtTcTc

  book_pipenu_histset[150] = 1;		// all  tracks, wt = pionSurvProb Q>0

  book_pipenu_histset[151] = 1;		// IDWord[0]=0 tracks, wt = pionSurvProb Q>0
  book_pipenu_histset[152] = 1;		// IDWord[0]=0 tracks, wt = pionSurvProb Q>0 p>60
  book_pipenu_histset[153] = 1;		// IDWord[0]=0 tracks, wt = pionSurvProb Q>0 p>60 T>300
  book_pipenu_histset[154] = 1;		// IDWord[0]=0 tracks, wt = pionSurvProb Q>0 p>60 T>300 DnTrackTc
  book_pipenu_histset[155] = 1;		// IDWord[0]=0 tracks, wt = pionSurvProb Q>0 p>60 T>300 DnTrackTc DtTcTc
  book_pipenu_histset[156] = 1;		// IDWord[0]=0 tracks, wt = pionSurvProb Q>0 T>300
  book_pipenu_histset[157] = 1;		// IDWord[0]=0 tracks, wt = pionSurvProb Q>0 T>300 DnTrackTc
  book_pipenu_histset[158] = 1;		// IDWord[0]=0 tracks, wt = pionSurvProb Q>0 T>300 DnTrackTc DtTcTc

  book_pipenu_histset[171] = 1;		// IDWord[1]=0 tracks, wt = pionSurvProb Q>0
  book_pipenu_histset[172] = 1;		// IDWord[1]=0 tracks, wt = pionSurvProb Q>0 p>60
  book_pipenu_histset[173] = 1;		// IDWord[1]=0 tracks, wt = pionSurvProb Q>0 p>60 T>300
  book_pipenu_histset[174] = 1;		// IDWord[1]=0 tracks, wt = pionSurvProb Q>0 p>60 T>300 DnTrackTc
  book_pipenu_histset[175] = 1;		// IDWord[1]=0 tracks, wt = pionSurvProb Q>0 p>60 T>300 DnTrackTc DtTcTc
  book_pipenu_histset[176] = 1;		// IDWord[1]=0 tracks, wt = pionSurvProb Q>0 T>300
  book_pipenu_histset[177] = 1;		// IDWord[1]=0 tracks, wt = pionSurvProb Q>0 T>300 DnTrackTc
  book_pipenu_histset[178] = 1;		// IDWord[1]=0 tracks, wt = pionSurvProb Q>0 T>300 DnTrackTc DtTcTc

  book_pipenu_histset[200] = 1;		// all         tracks Q<0

  book_pipenu_histset[201] = 1;		// IDWord[0]=0 tracks Q<0
  book_pipenu_histset[221] = 1;		// IDWord[1]=0 tracks Q<0

  book_pipenu_histset[250] = 1;		// all  tracks, wt = pion surv prob Q<0

  book_pipenu_histset[251] = 1;		// IDWord[0]=0 tracks, wt = pion surv prob Q<0
  book_pipenu_histset[271] = 1;		// IDWord[1]=0 tracks, wt = pion surv prob Q<0

  for (int i=0; i<kNPipenuHistSets; i++) {
    if (book_pipenu_histset[i] == 0)                        continue;
    sprintf(folder_name,"pip_%i",i);
    fol = (TFolder*) hist_folder->FindObject(folder_name);
    if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
    fHist.fPipenu[i] = new PipenuHist_t;
    BookPipenuHistograms(fHist.fPipenu[i],Form("Hist/%s",folder_name));
  }
}


//-----------------------------------------------------------------------------
// fill efficiency histograms : need 10 histogram sets
// pitch = 1./tan(dip)
//-----------------------------------------------------------------------------
void TPipenuAnaModule::FillEfficiencyHistograms(TStnTrackBlock*  TrackBlock, 
                                                TStnTrackID*     TrackID   , 
                                                int              HistSet   ) {
  // having MC truth is a must for calculating efficiency!
  if (fEvtPar.fSimp == nullptr) return;

  if (fEvtPar.fSimp->NStrawHits() >= 20) {
    FillEventHistograms(fHist.fEvent[HistSet],&fEvtPar);

    if (fEvtPar.fSimp->fMomTrackerFront > 100.) {
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

//-----------------------------------------------------------------------------
void    TPipenuAnaModule::FillPipenuHistograms    (PipenuHist_t*  Hist, 
                                                   TStnTrack*     Trk, 
                                                   TrackPar_t*    Tp, 
                                                   SimPar_t*      SimPar,
                                                   double         Weight) {
  Hist->fPVsT0->Fill(Trk->fT0,Tp->fP,Weight);
}

//-----------------------------------------------------------------------------
void    TPipenuAnaModule::FillTimeClusterHistograms(TimeClusterHist_t* Hist  , 
                                                    TStnTimeCluster*   Tc    ,
                                                    double             Weight) {
  Hist->fNsh->Fill  (Tc->NHits()     ,Weight);
  Hist->fNch->Fill  (Tc->NComboHits(),Weight);
  Hist->fT0->Fill   (Tc->T0()        ,Weight);
  Hist->fT0Err->Fill(Tc->T0Err()     ,Weight);
}

//_____________________________________________________________________________
void TPipenuAnaModule::FillHistograms() {

  double       cos_th (-2.); //,  cl_e(-1.);
  int          disk_id(-1);
//-----------------------------------------------------------------------------
// event histograms
//-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0],&fEvtPar);

  if (fEvtPar.fNTracks[0]> 0) FillEventHistograms(fHist.fEvent[1],&fEvtPar);
  else                        FillEventHistograms(fHist.fEvent[2],&fEvtPar);

  if (fNClusters > 0) FillEventHistograms(fHist.fEvent[3],&fEvtPar);
  else                FillEventHistograms(fHist.fEvent[4],&fEvtPar);

  if ((fEvtPar.fNTracks[0] == 0) && (fabs(cos_th) < 0.4)) {
    FillEventHistograms(fHist.fEvent[5],&fEvtPar); 
  }

  if (fEvtPar.fNGoodTracks[0] > 0) {
    FillEventHistograms(fHist.fEvent[6],&fEvtPar); 

    TLorentzVector    mom(1.,0.,0.,0);
    
    if (fSimPar.fParticle) mom = *fSimPar.fParticle->StartMom();

    double p, cos_th;

    p      = mom.P();
    cos_th = mom.Pz()/p;

    if (GetDebugBit(31) && (cos_th > 0.8)) {
      GetHeaderBlock()->Print(Form(" bit:031 cos_th = %10.3f p = %10.3f ntrk = %5i",
				   cos_th, p, fEvtPar.fNTracks[0]));
    }
  }

  if      (disk_id == 0) FillEventHistograms(fHist.fEvent[8],&fEvtPar);
  else if (disk_id == 1) FillEventHistograms(fHist.fEvent[9],&fEvtPar);
//-----------------------------------------------------------------------------
// EVT_10: events with SetC tracks 103.5 < p < 105.0
//-----------------------------------------------------------------------------
  if (fEvtPar.fNGoodTracks[0] > 0) {
    TStnTrack* trk = fTrackBlock->Track(0);
    if ((103.5 < trk->fP) && (trk->fP < 105.)) {
      FillEventHistograms(fHist.fEvent[10],&fEvtPar);
    }
  }
//-----------------------------------------------------------------------------
// helix histograms
// HEL_0: all
// HEL_1: helices in events with reconstructed tracks
//-----------------------------------------------------------------------------
  int nhel = fHelixBlock->NHelices();
  for (int i=0; i<nhel; i++) {
    TStnHelix*  hel = fHelixBlock->Helix(i);
    HelixPar_t* hp  = fHelixPar+i;
    FillHelixHistograms(fHist.fHelix[0],hel,hp);
    if (fEvtPar.fNTracks[0]  >  0) FillHelixHistograms(fHist.fHelix[1],hel,hp);
    if (hel->fTrackSeedIndex >= 0) FillHelixHistograms(fHist.fHelix[2],hel,hp);
  }
//-----------------------------------------------------------------------------
// SIMP histograms
//-----------------------------------------------------------------------------
  if (fEvtPar.fSimp) {
    SimpData_t sd;
    sd.fParent  = nullptr;
    sd.fStepVD9 = nullptr;

    TSimParticle* mom = fSimpBlock->Particle(fEvtPar.fSimp->ParentID());

    FillSimpHistograms(fHist.fSimp[0],fEvtPar.fSimp,&sd);
    if (fEvtPar.fNTracks[0] > 0) {
      FillSimpHistograms(fHist.fSimp[ 1],fEvtPar.fSimp,&sd);
      FillSimpHistograms(fHist.fSimp[31],mom,&sd);
      if (fEvtPar.fNGoodTracks[0] > 0) {
        FillSimpHistograms(fHist.fSimp[ 2],fEvtPar.fSimp,&sd);
        FillSimpHistograms(fHist.fSimp[32],mom,&sd);
      }

      if (fN300 > 0) {
        FillSimpHistograms(fHist.fSimp[33],mom,&sd);
      }
    }

    if (fEvtPar.fSimp->StartVolumeIndex() == 3015) { 
//-----------------------------------------------------------------------------
// last foil
//-----------------------------------------------------------------------------
      FillSimpHistograms(fHist.fSimp[10],fEvtPar.fSimp,&sd);
      if (fEvtPar.fNTracks[0] > 0) FillSimpHistograms(fHist.fSimp[11],fEvtPar.fSimp,&sd);
    }
  }
//-----------------------------------------------------------------------------
// track histograms, fill them only for the downstream e- hypothesis
//-----------------------------------------------------------------------------
  TStnTrack*   trk;
  TrackPar_t*  tp;

  for (int i=0; i<fEvtPar.fNTracks[0]; ++i ) {
    trk       = fTrackBlock->Track(i);
    tp        = fTrackPar+i;

    fTrackID[0]->FillHistograms(fHist.fTrackID[0],trk);

    tp->fDnTrackTc = tp->fTimeCluster->NHits()-trk->NActive();
    tp->fDtTrackTc = tp->fTimeCluster->NHits()-trk->NActive();
//-----------------------------------------------------------------------------
// find closest time cluster 
//-----------------------------------------------------------------------------
    float dt_tc_min = 1.e6;   // signed !!!

    for (int itc=0; itc<fEvtPar.fNTimeClusters; itc++) {
      TStnTimeCluster* tc = fTimeClusterBlock->TimeCluster(itc);
      if (tc == tp->fTimeCluster)                 continue;
      float dt = tp->fTimeCluster->T0()-tc->T0();
      if (fabs(dt) < fabs(dt_tc_min)) {
        dt_tc_min = dt;
      }
    }

    tp->fDtTcTc = dt_tc_min;

    FillTrackHistograms (fHist.fTrack [0],trk,tp,&fSimPar);
    FillPipenuHistograms(fHist.fPipenu[0],trk,tp,&fSimPar);
    if (trk->Charge() > 0) {
      FillTrackHistograms (fHist.fTrack [100],trk,tp,&fSimPar);
      FillPipenuHistograms(fHist.fPipenu[100],trk,tp,&fSimPar);
    }
    else {
      FillTrackHistograms (fHist.fTrack [200],trk,tp,&fSimPar);
      FillPipenuHistograms(fHist.fPipenu[200],trk,tp,&fSimPar);
    }

    if (tp->fIDWord[fBestID] == 0) {
					// track passes the quality selection

      FillTrackHistograms (fHist.fTrack [  1],trk,tp,&fSimPar);
      FillPipenuHistograms(fHist.fPipenu[  1],trk,tp,&fSimPar);
      if (trk->Charge() > 0) {
        FillTrackHistograms (fHist.fTrack [101],trk,tp,&fSimPar);
        FillPipenuHistograms(fHist.fPipenu[101],trk,tp,&fSimPar);
//-----------------------------------------------------------------------------
// at this point, apply the analysis-specific cuts
// 1. require the number of track hits to be consistent with the number of hits 
// 2. require no close time clusters
// TRK_102: P>60
// TRK_103: P>60 T>300
// TRK_104: P>60 T>300 dNTrackTC < fDnMax
// TRK_105: P>60 T>300 dNTrackTC < fDnMax abs(fDtTcTc) > 100
//-----------------------------------------------------------------------------
        if (tp->fP > 60) {
          FillTrackHistograms (fHist.fTrack [102],trk,tp,&fSimPar);
          FillPipenuHistograms(fHist.fPipenu[102],trk,tp,&fSimPar);
          if (trk->T0() > 300) {
            FillTrackHistograms (fHist.fTrack [103],trk,tp,&fSimPar);
            FillPipenuHistograms(fHist.fPipenu[103],trk,tp,&fSimPar);
            if (tp->fDnTrackTc <= fDnMax) { 
              FillTrackHistograms (fHist.fTrack [104],trk,tp,&fSimPar);
              FillPipenuHistograms(fHist.fPipenu[104],trk,tp,&fSimPar);
              if (fabs(tp->fDtTcTc) > 100) {
                fEventPassedSelections = 1;
                FillTrackHistograms (fHist.fTrack [105],trk,tp,&fSimPar);
                FillPipenuHistograms(fHist.fPipenu[105],trk,tp,&fSimPar);
              }
            }
          }
        }
//-----------------------------------------------------------------------------
// the same, but w/o p>60 MeV/c cut
//-----------------------------------------------------------------------------
        if (trk->T0() > 300) {
          FillTrackHistograms (fHist.fTrack [106],trk,tp,&fSimPar);
          FillPipenuHistograms(fHist.fPipenu[106],trk,tp,&fSimPar);
          if (tp->fDnTrackTc <= fDnMax) { 
            FillTrackHistograms (fHist.fTrack [107],trk,tp,&fSimPar);
            FillPipenuHistograms(fHist.fPipenu[107],trk,tp,&fSimPar);
            if (fabs(tp->fDtTcTc) > 100) {
              fEventPassedSelections = 1;
              FillTrackHistograms (fHist.fTrack [108],trk,tp,&fSimPar);
              FillPipenuHistograms(fHist.fPipenu[108],trk,tp,&fSimPar);
            }
          }
        }
      }
      else {
        FillTrackHistograms (fHist.fTrack [201],trk,tp,&fSimPar);
        FillPipenuHistograms(fHist.fPipenu[201],trk,tp,&fSimPar);
      }
    }
//-----------------------------------------------------------------------------
// tight ID
//-----------------------------------------------------------------------------
    if (tp->fIDWord[1] == 0) {
					// track passes the quality selection

      FillTrackHistograms (fHist.fTrack [ 21],trk,tp,&fSimPar);
      FillPipenuHistograms(fHist.fPipenu[ 21],trk,tp,&fSimPar);
      if (trk->Charge() > 0) {
        FillTrackHistograms (fHist.fTrack [121],trk,tp,&fSimPar);
        FillPipenuHistograms(fHist.fPipenu[121],trk,tp,&fSimPar);
//-----------------------------------------------------------------------------
// at this point, apply the analysis-specific cuts
// 1. require the number of track hits to be consistent with the number of hits 
// 2. require no close time clusters
// TRK_102: P>60
// TRK_103: P>60 T>300
// TRK_104: P>60 T>300 dNTrackTC < fDnMax
// TRK_105: P>60 T>300 dNTrackTC < fDnMax abs(fDtTcTc) > 100
//-----------------------------------------------------------------------------
        if (tp->fP > 60) {
          FillTrackHistograms (fHist.fTrack [122],trk,tp,&fSimPar);
          FillPipenuHistograms(fHist.fPipenu[122],trk,tp,&fSimPar);
          if (trk->T0() > 300) {
            FillTrackHistograms (fHist.fTrack [123],trk,tp,&fSimPar);
            FillPipenuHistograms(fHist.fPipenu[123],trk,tp,&fSimPar);
            if (tp->fDnTrackTc <= fDnMax) { 
              FillTrackHistograms (fHist.fTrack [124],trk,tp,&fSimPar);
              FillPipenuHistograms(fHist.fPipenu[124],trk,tp,&fSimPar);
              if (fabs(tp->fDtTcTc) > 100) {
                fEventPassedSelections = 1;
                FillTrackHistograms (fHist.fTrack [125],trk,tp,&fSimPar);
                FillPipenuHistograms(fHist.fPipenu[125],trk,tp,&fSimPar);
              }
            }
          }
        }
//-----------------------------------------------------------------------------
// the same, but w/o p>60 MeV/c cut
//-----------------------------------------------------------------------------
        if (trk->T0() > 300) {
          FillTrackHistograms (fHist.fTrack [126],trk,tp,&fSimPar);
          FillPipenuHistograms(fHist.fPipenu[126],trk,tp,&fSimPar);
          if (tp->fDnTrackTc <= fDnMax) { 
            FillTrackHistograms (fHist.fTrack [127],trk,tp,&fSimPar);
            FillPipenuHistograms(fHist.fPipenu[127],trk,tp,&fSimPar);
            if (fabs(tp->fDtTcTc) > 100) {
              fEventPassedSelections = 1;
              FillTrackHistograms (fHist.fTrack [128],trk,tp,&fSimPar);
              FillPipenuHistograms(fHist.fPipenu[128],trk,tp,&fSimPar);
            }
          }
        }
      }
      else {
        FillTrackHistograms (fHist.fTrack [221],trk,tp,&fSimPar);
        FillPipenuHistograms(fHist.fPipenu[221],trk,tp,&fSimPar);
      }
    }
//-----------------------------------------------------------------------------
// track reco efficiency - filles Event histograms
//-----------------------------------------------------------------------------
    FillEfficiencyHistograms(fTrackBlock,TAnaModule::fTrackID[fBestID],11);
//-----------------------------------------------------------------------------
// TRK_101: all tracks pi+ --> e+ nu weighted with the pion survival prob
//-----------------------------------------------------------------------------
    if (fEvtPar.fPionSurvProb > 0) {
      FillTrackHistograms (fHist.fTrack [ 50],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
      FillPipenuHistograms(fHist.fPipenu[ 50],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
      if (trk->Charge() > 0) {
        FillTrackHistograms (fHist.fTrack [150],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
        FillPipenuHistograms(fHist.fPipenu[150],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
      }
      else {
        FillTrackHistograms (fHist.fTrack [250],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
        FillPipenuHistograms(fHist.fPipenu[250],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
      }
      
      if (tp->fIDWord[fBestID] == 0) {
        FillTrackHistograms (fHist.fTrack [ 51],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
        FillPipenuHistograms(fHist.fPipenu[ 51],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
        if (trk->Charge() > 0) {
          FillTrackHistograms (fHist.fTrack [151],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
          FillPipenuHistograms(fHist.fPipenu[151],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
          if (tp->fP > 60) {
            FillTrackHistograms (fHist.fTrack [152],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
            FillPipenuHistograms(fHist.fPipenu[152],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
            if (trk->T0() > 300) {
              FillTrackHistograms (fHist.fTrack [153],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
              FillPipenuHistograms(fHist.fPipenu[153],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
//-----------------------------------------------------------------------------
// at this point, apply the analysis-specific cuts
// 1. require the number of track hits to be consistent with the number of hits 
// 2. require no close time clusters
//-----------------------------------------------------------------------------
              if (tp->fDnTrackTc <= fDnMax) {
                FillTrackHistograms (fHist.fTrack [154],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
                FillPipenuHistograms(fHist.fPipenu[154],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
                if (fabs(tp->fDtTcTc) > 100) {
                  FillTrackHistograms (fHist.fTrack [155],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
                  FillPipenuHistograms(fHist.fPipenu[155],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
                }
              }
            }
          }
//-----------------------------------------------------------------------------
// the same but w/o p>60 MeV/c cut
//-----------------------------------------------------------------------------
          if (trk->T0() > 300) {
            FillTrackHistograms (fHist.fTrack [156],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
            FillPipenuHistograms(fHist.fPipenu[156],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
            if (tp->fDnTrackTc <= fDnMax) {
              FillTrackHistograms (fHist.fTrack [157],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
              FillPipenuHistograms(fHist.fPipenu[157],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
              if (fabs(tp->fDtTcTc) > 100) {
                FillTrackHistograms (fHist.fTrack [158],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
                FillPipenuHistograms(fHist.fPipenu[158],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
              }
            }
          }
        }
        else {
          FillTrackHistograms (fHist.fTrack [251],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
          FillPipenuHistograms(fHist.fPipenu[251],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
        }
      }
//-----------------------------------------------------------------------------
// tight ID
//-----------------------------------------------------------------------------
      if (tp->fIDWord[1] == 0) {
        FillTrackHistograms (fHist.fTrack [ 71],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
        FillPipenuHistograms(fHist.fPipenu[ 71],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
        if (trk->Charge() > 0) {
          FillTrackHistograms (fHist.fTrack [171],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
          FillPipenuHistograms(fHist.fPipenu[171],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
          if (tp->fP > 60) {
            FillTrackHistograms (fHist.fTrack [172],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
            FillPipenuHistograms(fHist.fPipenu[172],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
            if (trk->T0() > 300) {
              FillTrackHistograms (fHist.fTrack [173],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
              FillPipenuHistograms(fHist.fPipenu[173],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
//-----------------------------------------------------------------------------
// at this point, apply the analysis-specific cuts
// 1. require the number of track hits to be consistent with the number of hits 
// 2. require no close time clusters
//-----------------------------------------------------------------------------
              if (tp->fDnTrackTc <= fDnMax) {
                FillTrackHistograms (fHist.fTrack [174],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
                FillPipenuHistograms(fHist.fPipenu[174],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
                if (fabs(tp->fDtTcTc) > 100) {
                  FillTrackHistograms (fHist.fTrack [175],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
                  FillPipenuHistograms(fHist.fPipenu[175],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
                }
              }
            }
          }
//-----------------------------------------------------------------------------
// the same but w/o p>60 MeV/c cut
//-----------------------------------------------------------------------------
          if (trk->T0() > 300) {
            FillTrackHistograms (fHist.fTrack [176],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
            FillPipenuHistograms(fHist.fPipenu[176],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
            if (tp->fDnTrackTc <= fDnMax) {
              FillTrackHistograms (fHist.fTrack [177],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
              FillPipenuHistograms(fHist.fPipenu[177],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
              if (fabs(tp->fDtTcTc) > 100) {
                FillTrackHistograms (fHist.fTrack [178],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
                FillPipenuHistograms(fHist.fPipenu[178],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
              }
            }
          }
        }
        else {
          FillTrackHistograms (fHist.fTrack [271],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
          FillPipenuHistograms(fHist.fPipenu[271],trk,tp,&fSimPar,fEvtPar.fPionSurvProb);
        }
      }
    }
  }
//-----------------------------------------------------------------------------
// fill track ID histograms
//-----------------------------------------------------------------------------
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
int TPipenuAnaModule::Event(int ientry) {

  //  TDiskCalorimeter::GeomData_t disk_geom;

  fTrackBlock->GetEntry(ientry);
  fTrackSeedBlock->GetEntry(ientry);
  fHelixBlock->GetEntry(ientry);
  fTimeClusterBlock->GetEntry(ientry);
  //  fTrackStrawHitBlock->GetEntry(ientry);
  // fClusterBlock->GetEntry(ientry);
  // fStrawHitBlock->GetEntry(ientry);
  // fCalDataBlock->GetEntry(ientry);
  fGenpBlock->GetEntry(ientry);
  fSimpBlock->GetEntry(ientry);
  fVDetBlock->GetEntry(ientry);
//-----------------------------------------------------------------------------
// at some point, the TVdetBlock class got obsolete, and now the virtual detector 
// hits are stored in TStepPointMCBlock 
// dont' try to read it as that would fail
//-----------------------------------------------------------------------------
  fEvtPar.fDioLOWt          = 1.;
  fEvtPar.fDioLLWt          = 1.;
  fEvtPar.fNCrvClusters     = -1;
  fEvtPar.fNCrvPulses       = -1;
  fEvtPar.fNCrvCoincidences = -1;
  fEvtPar.fNStrawHits       = GetHeaderBlock()->NStrawHits();
  fEvtPar.fNComboHits       = GetHeaderBlock()->NComboHits();
  fEvtPar.fNGenp            = fGenpBlock->NParticles();
  fEvtPar.fNSimp            = fSimpBlock->NParticles();
  fEvtPar.fSimp             = nullptr;
  fEvtPar.fPartE            = -1.;
  fEvtPar.fPionSurvProb     = 1.;
  fEvtPar.fNHelices         = fHelixBlock->NHelices();
  fEvtPar.fNTimeClusters    = fTimeClusterBlock->NTimeClusters();
//-----------------------------------------------------------------------------
// luminosity weight
//-----------------------------------------------------------------------------
  fLumWt = GetHeaderBlock()->LumWeight();
//-----------------------------------------------------------------------------
// figure the MC truth
// pi+ --> e+ nu case : determine the event weight
//-----------------------------------------------------------------------------
  for (int i=fEvtPar.fNSimp-1; i>=0; i--) {
    TSimParticle* simp = fSimpBlock->Particle(i);
    int pdg_code       = simp->PDGCode();
    int generator_id   = simp->GeneratorID();         // MC process code

    if ((pdg_code == fPDGCode) and (generator_id == fMCProcessCode)) {
      fEvtPar.fSimp  = simp;
      fEvtPar.fPartE = simp->StartMom()->Energy();
    }

    if ((abs(pdg_code) == 211) && (simp->GeneratorID() == 56)) {
//-----------------------------------------------------------------------------
// found the pion, survival probability
//-----------------------------------------------------------------------------
      fEvtPar.fPionSurvProb  = exp(-simp->fEndProperTime);
      break;
    }
  }
//-----------------------------------------------------------------------------
// calculate DIO weights once per event
//-----------------------------------------------------------------------------
  if (fEvtPar.fPartE > 0) {
    fEvtPar.fDioLOWt = TStntuple::DioWeightAl   (fEvtPar.fPartE);
    fEvtPar.fDioLLWt = TStntuple::DioWeightAl_LL(fEvtPar.fPartE);
  }
//-----------------------------------------------------------------------------
// may want to revisit the definition of fSimp in the future
//-----------------------------------------------------------------------------
  fSimPar.fParticle = fEvtPar.fSimp;
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
  int nvdhits = fVDetBlock->NStepPoints();
  for (int i=0; i<nvdhits; i++) {
    TStepPointMC* vdhit = fVDetBlock->StepPointMC(i);
    if (vdhit->PDGCode() == fEvtPar.fSimp->fPdgCode) {
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
//-----------------------------------------------------------------------------
// init calorimeter
//-----------------------------------------------------------------------------
  fNClusters  = fClusterBlock->NClusters();
  fNCalHits   = fCalDataBlock->NHits();
//-----------------------------------------------------------------------------
// init calorimeter clusters - remember, the first one not necessarily is the 
// most energetic
//-----------------------------------------------------------------------------
  fNClusters = fClusterBlock->NClusters();
  if (fNClusters == 0) fCluster = 0;
  else                 fCluster = fClusterBlock->Cluster(0);
//-----------------------------------------------------------------------------
// loop over tracks and calculate needed parameters
//-----------------------------------------------------------------------------
  fNStrawHits          = GetHeaderBlock()->fNStrawHits;

  fEvtPar.fNTracks[0]  = fTrackBlock->NTracks();
  if (fEvtPar.fNTracks[0] == 0) fTrack = 0;
  else                          fTrack = fTrackBlock->Track(0);

  fEvtPar.fNGoodTracks[0] = 0;
  fEvtPar.fNGoodTracks[1] = 0;
  fNMatchedTracks         = 0;
//-----------------------------------------------------------------------------
// determine the number of CalPatRec tracks - this assumes that the list of 
// tracks has been created by MergePatRec
//-----------------------------------------------------------------------------
  int ntrk = fTrackBlock->NTracks();

  InitTrackPar(fTrackBlock,fClusterBlock,fTrackPar,&fSimPar);
//-----------------------------------------------------------------------------
// additional initializations - helices and time clusters
//-----------------------------------------------------------------------------
  fN300 = 0;
  for (int itrk=0; itrk<ntrk; itrk++) {
    TStnTrack*   trk = fTrackBlock->Track(itrk);
    TrackPar_t*  tp  = fTrackPar+itrk;

    tp->fHelix       = nullptr;
    tp->fTimeCluster = nullptr;

    int ih = trk->fHelixIndex;
    if (ih < fEvtPar.fNHelices) {
      tp->fHelix = fHelixBlock->Helix(ih);
    }
//-----------------------------------------------------------------------------
// fN300 : the number of late (T0>300) good tracks above 60 MeV/c
//-----------------------------------------------------------------------------
    if ((tp->fIDWord[fBestID] == 0) and (tp->fP > 60) and (trk->T0() > 300)) { 
      fN300++;
    }
//-----------------------------------------------------------------------------
// currently, at least in the KK case, an TStnHelix doesn't have a helix index stored
// find a time cluster closest to teh track in time
//-----------------------------------------------------------------------------
    TStnTimeCluster* closest_tc(nullptr);
    float min_dt = 1.e6;
    for (int itc=0; itc<fEvtPar.fNTimeClusters; itc++) {
      TStnTimeCluster* tc = fTimeClusterBlock->TimeCluster(itc);
      float dt = trk->T0()-tc->T0();
      if (fabs(dt) < fabs(min_dt)) {
        closest_tc = tc;
        min_dt     = dt;
      }
    }

    tp->fTimeCluster = closest_tc;
    tp->fDtTrackTc   = min_dt;
  }

  fEventPassedSelections = 0;

  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TPipenuAnaModule::Debug() {

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
      GetHeaderBlock()->Print(Form("bit_000: All p = %10.3f",trk->Momentum()->P()));
    }
//-----------------------------------------------------------------------------
// bit 1: IDWord =0 0 tracks
//-----------------------------------------------------------------------------
    if (GetDebugBit(1) == 1) {
      if (tp->fIDWord[fBestID] == 0) {
	GetHeaderBlock()->Print(Form("bit_001: IDWord=0 p = %10.3f",trk->Momentum()->P()));
      }
    }
//-----------------------------------------------------------------------------
// bit 3: good tracks with P > 70 T0 > 300
//-----------------------------------------------------------------------------
    if (GetDebugBit(3) == 1) {
      if (fEventPassedSelections == 1) {
        if ((tp->fP > 70) and (trk->fT0 > 300)) {
	  GetHeaderBlock()->Print(Form("large P: %f",tp->fP));
	}
      }
    }
//-----------------------------------------------------------------------------
// bit 4: good events with 65 < P < 70 , T0 > 300
//-----------------------------------------------------------------------------
    if (GetDebugBit(4) == 1) {
      if (tp->fIDWord[fBestID] == 0) {
        if ((tp->fDnTrackTc <= fDnMax) and (fabs(tp->fDtTcTc) > 100)) {
          if ((tp->fP > 65) and (tp->fP < 700) and (trk->fT0 > 300)) {
            GetHeaderBlock()->Print(Form(":bit_04: candidate: i=%i trk->P0=%10.3f tp->fP=%10.3f",itrk,trk->P0(),tp->fP));
          }
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
int TPipenuAnaModule::EndJob() {
  printf("----- end job: ---- %s\n",GetName());
  return 0;
}

//_____________________________________________________________________________
void TPipenuAnaModule::Test001() {
}

}
