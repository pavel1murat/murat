//////////////////////////////////////////////////////////////////////////////
// 2020-12-26 P.Murat
// DAR fits only, compare downstream-only reco for electron and muon events 
//
// use of tmp:
//
// use of debug bits: bits 0-2 are reserved
//  0  : all events
//  1  : passed events
//  2  : rejected events
//
// call: "temu_ana(28,4)
///////////////////////////////////////////////////////////////////////////////
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TEnv.h"
#include "TSystem.h"

#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/obj/TStnNode.hh"
#include "Stntuple/alg/TStntuple.hh"
#include "Stntuple/val/stntuple_val_functions.hh"

#include <string>
#include "math.h"

#include "murat/ana/TEmuModule.hh"

using std::string;
using std::vector;

ClassImp(murat::TEmuModule)

namespace murat {
//-----------------------------------------------------------------------------
TEmuModule::TEmuModule(const char* name, const char* title): TAnaModule(name,title) {

  fPdgCode           = 11;		// electron
  fGeneratorCode     =  2;              // 2:ConversionGun 28:StoppedParticleReactionGun

  fTrackBlockName[0] = "TrackBlockDarDe";
  fTrackBlockName[1] = "TrackBlockDarDmu";
//-----------------------------------------------------------------------------
// TrackID[0]     : Michael's cuts
// fTrackID[1-19] : cut on MVA-based TrkQual with the step of 0.05
//-----------------------------------------------------------------------------
  fBestID[0] = 16;                      // default best for PAR tracks ( > 0.05*15 = 0.8)
  fBestID[1] = 0;                       // default best for DAR tracks 

  fNMVA      = 1; 			// number of MVA's used (if used at all)

  fNID       = 1;                      // fNID has to be < TAnaModule::kMaxNTrackID

  for (int i=0; i<fNID; i++) {
    fTrackID[i] = new TStnTrackID();

    fTrackID[i]->SetMaxMomErr (100);
    fTrackID[i]->SetMaxT0Err  (100);
    fTrackID[i]->SetMinNActive( -1);
    fTrackID[i]->SetMinFitCons( -1);
    fTrackID[i]->SetMinTanDip (0.5);
    fTrackID[i]->SetMaxTanDip (1.0);
    fTrackID[i]->SetMinTrkQual(0.2);
    fTrackID[i]->SetMinD0     (-100.);
    fTrackID[i]->SetMaxD0     ( 100.);

    int mask = TStnTrackID::kT0Bit | TStnTrackID::kTanDipBit | TStnTrackID::kD0Bit | TStnTrackID::kTrkQualBit  ;

    fTrackID[0]->SetUseMask(mask);
  }
}

//-----------------------------------------------------------------------------
TEmuModule::~TEmuModule() {
  for (int i=0; i<fNID; i++) delete fTrackID[i];
}

//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TEmuModule::BeginJob() {

  TAnaModule::BeginJob();
//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
  RegisterDataBlock(fTrackBlockName[0].Data(), "TStnTrackBlock"   , &fTrackBlock[0] );
  RegisterDataBlock(fTrackBlockName[1].Data(), "TStnTrackBlock"   , &fTrackBlock[1] );
  RegisterDataBlock("ClusterBlock"           , "TStnClusterBlock" , &fClusterBlock  );
  RegisterDataBlock("SimpBlock"              , "TSimpBlock"       , &fSimpBlock     );
  RegisterDataBlock("GenpBlock"              , "TGenpBlock"       , &fGenpBlock     );
  RegisterDataBlock("SpmcBlockVDet"          , "TStepPointMCBlock", &fSpmcBlockVDet );
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();
//-----------------------------------------------------------------------------
// force using track quality MVA (DAR-only - code 1070)
//-----------------------------------------------------------------------------
  const char* training_dataset = "fele2s51b1";
  int         training_code    = 1070;

  printf(" [TEmuModule::BeginJob] TrainingDataset:%s MvaType:%i\n",training_dataset,training_code);

  TAnaModule::fUseMVA = 1;

  for (int i=0; i<2; i++) { 
    if (fTrkQualMVA[i]) delete fTrkQualMVA[i];
    fTrkQualMVA[i] = new mva_data(training_dataset,training_code);
  }
//-----------------------------------------------------------------------------
// init track ID pointers in TrackPar - do it just once
//  .. assume TrkQual to be precalculated, do only beed to cut on
//-----------------------------------------------------------------------------
  for (int i=0; i<2; i++) {
    for (int ip=0; ip<kNTrackPar; ip++) {
      TrackPar_t* tp = &fTrackPar[i][ip];
      for (int id=0; id<fNID; id++) {
	fTrackID[id]->SetLocTrkQual(0);              // use TStnTrack::fTmp[0] to store MVA-calculated TrkQual
	tp->fTrackID[id] = fTrackID[id];
      }
      tp->fLogLH      = TAnaModule::fLogLH;          // do it just once
    }

    if (fUseMVA) {
      if (fTrkQualMVA[i]) {
	fBestID[i] = fTrkQualMVA[i]->BestID();		  // Dave's default: DaveTrkQual > 0.8 ; or else
      }
    }
  }

  return 0;
}


//_____________________________________________________________________________
int TEmuModule::BeginRun() {
  int rn = GetHeaderBlock()->RunNumber();
  TStntuple::Init(rn);
  return 0;
}


//_____________________________________________________________________________
void TEmuModule::BookHistograms() {

  //  char name [200];
  //  char title[200];

  TFolder* fol;
  TFolder* hist_folder;
  char     folder_name[200];

  DeleteHistograms();

  TH1::SetDefaultSumw2(kTRUE);	

  hist_folder = (TFolder*) GetFolder()->FindObject("Hist");
//-----------------------------------------------------------------------------
// book event histograms
//-----------------------------------------------------------------------------
  int book_event_histset[kNEventHistSets];
  for (int i=0; i<kNEventHistSets; i++) book_event_histset[i] = 0;

  book_event_histset[ 0] = 1;		// all events
  book_event_histset[ 1] = 1;		// events with EclMax > fMinETrig and TClMax > 550
  book_event_histset[ 2] = 1;           // *** fill in

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
  TString*  track_selection   [kNTrackHistSets];

  for (int i=0; i<kNTrackHistSets; i++) { track_selection[i] = NULL; }

  track_selection[  0] = new TString("good [0] tracks");
  track_selection[  1] = new TString("good [1] tracks in events with no good [0] tracks TrkQual>0.3");
  track_selection[  2] = new TString("good [1] tracks in events with no good [0] tracks TrkQual>0.2");

  track_selection[100] = new TString("[0] all tracks");
  track_selection[101] = new TString("[0] BestTrackID");
  track_selection[102] = new TString("[0] BestTrackID, cluster");
  track_selection[103] = new TString("[0] BestTrackID, cluster, first disk");
  track_selection[104] = new TString("[0] BestTrackID, cluster, second disk");
  track_selection[105] = new TString("[0] BestTrackID tracks events Ecl > 60");
  track_selection[106] = new TString("[0] tracks with dtz0 < -10");
  track_selection[109] = new TString("[0] tracks with XDpF > +10 MeV");

  track_selection[110] = new TString("[0] TrackID[0]");

  track_selection[200] = new TString("[1] all tracks");
  track_selection[201] = new TString("[1] BestTrackID");

  track_selection[202] = new TString("[1] BestTrackID, cluster");
  track_selection[203] = new TString("[1] BestTrackID, cluster, first disk");
  track_selection[204] = new TString("[1] BestTrackID, cluster, second disk");

  track_selection[205] = new TString("[1] BestTrackID tracks events Ecl > 60");
  track_selection[206] = new TString("[1] tracks with dtz0 < -10");
  track_selection[209] = new TString("[1] tracks with XDpF > +10 MeV");

  track_selection[210] = new TString("[1] TrackID[0]");

  const char* folder_title;
  for (int i=0; i<kNTrackHistSets; i++) {
    if (track_selection[i] != 0) {
      sprintf(folder_name,"trk_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      folder_title = folder_name;
      folder_title = track_selection[i]->Data();
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_title);
      fHist.fTrack[i] = new TrackHist_t;
      BookTrackHistograms(fHist.fTrack[i],Form("Hist/%s",folder_name));
    }
  }

}

//_____________________________________________________________________________

void TEmuModule::FillHistograms() {

  TStnTrack*   trk;
  TrackPar_t*  tp;
  int          ihist, n_setc_tracks[2];

  int bid_par = fBestID[kELE];
  int bid_dar = fBestID[kMUO];
//-----------------------------------------------------------------------------
// event histograms
// EVT_0: all events
//-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0],&fEvtPar);
//-----------------------------------------------------------------------------
// EVT_1: events with 50 MeV+ cluster and T0 > 550
//-----------------------------------------------------------------------------
  if ((fEClMax > 50.) && (fTClMax > 550)) {
    FillEventHistograms(fHist.fEvent[1],&fEvtPar);
  }
//-----------------------------------------------------------------------------
// EVT_2: events with 50 MeV+ cluster and both tracks T0 > 550
//-----------------------------------------------------------------------------
  if ((fEClMax > fMinETrig) && (fTClMax > 550)) {
    FillEventHistograms(fHist.fEvent[2],&fEvtPar);
  }
//-----------------------------------------------------------------------------
// what does kMUO add ?
// TRK_0 : kELE tracks, BEST_ID
//-----------------------------------------------------------------------------
  if ((fTrackBlock[kELE]->NTracks() > 0) && (fTrackPar[kELE][0].fIDWord[bid_par] == 0)) {
    TStnTrack* trk = fTrackBlock[kELE]->Track(0);
    TrackPar_t* tp = &fTrackPar[kELE][0];
    FillTrackHistograms(fHist.fTrack[0],trk,tp,&fSimPar);
  }
  else {
    if (fTrackBlock[kMUO]->NTracks() > 0) {
      if (fTrackPar[kMUO][0].fIDWord[bid_dar] == 0) {
//-----------------------------------------------------------------------------
// TRK_1: have KDAR track TrkQual > 0.3 and no KPAR TrkQual > 0.4 track
//-----------------------------------------------------------------------------
	TStnTrack* trk = fTrackBlock[kMUO]->Track(0);
	TrackPar_t* tp = &fTrackPar[kMUO][0];
	FillTrackHistograms(fHist.fTrack[1],trk,tp,&fSimPar);
      }
      if (fTrackPar[1][0].fIDWord[2] == 0) {
//-----------------------------------------------------------------------------
// TRK_2: have KDAR track TrkQual > 0.2 and no KPAR TrkQual > 0.4 track
//-----------------------------------------------------------------------------
	TStnTrack* trk = fTrackBlock[1]->Track(0);
	TrackPar_t* tp = &fTrackPar[1][0];
	FillTrackHistograms(fHist.fTrack[2],trk,tp,&fSimPar);
      }
    }
  }
//-----------------------------------------------------------------------------
// KPAR and KDAR histograms, inclusive, ihist defines the offset
// i=0:KPAR, i=1:KDAR
//-----------------------------------------------------------------------------
  for (int i=0; i<2; i++) {
    n_setc_tracks[i] = 0;
    ihist            = 100*(i+1);
    int best_id      = fBestID[i];

    for (int itrk=0; itrk<fNTracks[i]; itrk++) {
      trk = fTrackBlock[i]->Track(itrk);
      tp  = fTrackPar  [i]+itrk;
//-----------------------------------------------------------------------------
// TRK_100, TRK_200: all reconstructed tracks
//-----------------------------------------------------------------------------
      FillTrackHistograms(fHist.fTrack[ihist+0],trk,tp,&fSimPar);
//-----------------------------------------------------------------------------
// TRK_101: BestID 
// TRK_102: BestID + cluster
// TRK_103: BestID + cluster + 1st disk
// TRK_104: BestID + cluster + 2nd disk
//-----------------------------------------------------------------------------
      if (tp->fIDWord[best_id] == 0) {
	FillTrackHistograms(fHist.fTrack[ihist+1],trk,tp,&fSimPar);
	n_setc_tracks[i] += 1;

	if (tp->fDiskID >= 0) {
	  FillTrackHistograms(fHist.fTrack[ihist+2],trk,tp,&fSimPar);
	  if      (tp->fDiskID == 0) FillTrackHistograms(fHist.fTrack[ihist+3],trk,tp,&fSimPar);
	  else if (tp->fDiskID == 1) FillTrackHistograms(fHist.fTrack[ihist+4],trk,tp,&fSimPar);
	}
      }
//-----------------------------------------------------------------------------
// IHIST+4: add   Ecl > 60 requirement
//-----------------------------------------------------------------------------
      if (fEClMax > 60.) {
	if (tp->fIDWord[0] == 0) {
	  FillTrackHistograms(fHist.fTrack[ihist+5],trk,tp,&fSimPar);
	}
      }
//-----------------------------------------------------------------------------
// IHIST+6: oddly, tracks with seemingly wrong reconstruced T0
//-----------------------------------------------------------------------------
      if (tp->fDtZ0 < -10) {
	FillTrackHistograms(fHist.fTrack[ihist+6],trk,tp,&fSimPar);
      }
//-----------------------------------------------------------------------------
// IHIST+9: tracks with XDpF > 10
//-----------------------------------------------------------------------------
      if (tp->fXDpF   > 10) FillTrackHistograms(fHist.fTrack[ihist+9],trk,tp,&fSimPar);
//-----------------------------------------------------------------------------
// different cuts on track quality variable
//-----------------------------------------------------------------------------
      for (int idd=0; idd<fNID; idd++) {
	if (tp->fIDWord[idd] == 0) FillTrackHistograms(fHist.fTrack[ihist+10+idd],trk,tp,&fSimPar);
      }
    }
  }
}

//-----------------------------------------------------------------------------
// 2014-04-30: it looks that reading the straw hits takes a lot of time - 
//              turn off by default by commenting it out
//-----------------------------------------------------------------------------
int TEmuModule::Event(int ientry) {

  fTrackBlock[0]->GetEntry(ientry);
  fTrackBlock[1]->GetEntry(ientry);
  fClusterBlock->GetEntry(ientry);
  fSimpBlock->GetEntry(ientry);
  fGenpBlock->GetEntry(ientry);
  fSpmcBlockVDet->GetEntry(ientry);
//-----------------------------------------------------------------------------
// assume electron in the first particle, otherwise the logic will need to
// be changed
//-----------------------------------------------------------------------------
  fEvtPar.fNGenp  = fGenpBlock->NParticles();

  fNSimp       = fSimpBlock->NParticles();
  fNClusters   = fClusterBlock->NClusters();
  // fNHelices    = fHelixBlock->NHelices();
  // fNTrackSeeds = fTrackSeedBlock->NTrackSeeds();

  fCluster = NULL;
  fEClMax  = -1;
  fTClMax  = -1;
  if (fNClusters > 0) {
    fCluster = fClusterBlock->Cluster(0);
    fEClMax  = fCluster->Energy();
    fTClMax  = fCluster->Time  ();
  }
//-----------------------------------------------------------------------------
// MC generator info
//-----------------------------------------------------------------------------
  TGenParticle* genp;
  int           pdg_code, generator_code;

  fEvtPar.fParticle = NULL;
  for (int i=fEvtPar.fNGenp-1; i>=0; i--) {
    genp           = fGenpBlock->Particle(i);
    pdg_code       = genp->GetPdgCode();
    generator_code = genp->GetStatusCode();
    if ((abs(pdg_code) == fPdgCode) && (generator_code == fGeneratorCode)) {
      fEvtPar.fParticle = genp;
      break;
    }
  }

  //  fProcess = -1;
  fWeight  =  1.;
  //  fPhotonE = -1.;

  // if (fEvtPar.fNGenp > 0) {
  //   TGenParticle* p0 = fGenpBlock->Particle(0);
  // }
//-----------------------------------------------------------------------------
// may want to revisit the definition of fSimPar in future
//-----------------------------------------------------------------------------
  fSimPar.fParticle = fSimpBlock->Particle(0);
  fSimPar.fTFront   = NULL;
  fSimPar.fTMid     = NULL;
  fSimPar.fTBack    = NULL;
  fSimPar.fGenp     = NULL;
//-----------------------------------------------------------------------------
// virtual detectors - for fSimPar need parameters at the tracker front
//-----------------------------------------------------------------------------
  int nsteps = fSpmcBlockVDet->NStepPoints();
  
  for (int i=0; i<nsteps; i++) {
    TStepPointMC* step = fSpmcBlockVDet->StepPointMC(i);
    if (step->PDGCode() == fSimPar.fParticle->PDGCode()) {
      if ((step->VolumeID() == 13) || (step->VolumeID() == 14)) {
	fSimPar.fTFront = step;
      }
      else if ((step->VolumeID() == 11) || (step->VolumeID() == 12)) {
	fSimPar.fTMid = step;
      }
    }
  }
//-----------------------------------------------------------------------------
// initialize additional track parameters
//-----------------------------------------------------------------------------
  for (int i=0; i<2; i++) {
    fNTracks    [i] = fTrackBlock[i]->NTracks();
    fNGoodTracks[i] = 0;

    for (int itrk=0; itrk<fNTracks[i]; itrk++) {
      TrackPar_t* tp  = fTrackPar[i]+itrk;
      TStnTrack*  trk = fTrackBlock[i]->Track(itrk);

      tp->fAlg        = i;
      if (fTrkQualMVA[i] != NULL) {
					// use calculated on the fly MVA estimator
	trk->SetITmp(0,0);
      }
    }

    InitTrackPar(fTrackBlock[i],fClusterBlock,fTrackPar[i],&fSimPar);
  }

  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
// looking mostly at the KDAR tracks
//-----------------------------------------------------------------------------
// diagnostics is related to KDAR
//-----------------------------------------------------------------------------
void TEmuModule::Debug() {

  // TStnTrack* trk;
  // TrackPar_t* tp(NULL);
  //  char        text[500];
//-----------------------------------------------------------------------------
// bit 0: All Events
//-----------------------------------------------------------------------------
  if (GetDebugBit(0) == 1) {
    GetHeaderBlock()->Print(Form("TEmuModule :bit000:"));
    printf("KPAR:\n");
    fTrackBlock[kELE]->Print();

    for (int i=0; i<fNTracks[kELE]; i++) { 
      PrintTrack(fTrackBlock[kELE]->Track(i),&fTrackPar[kELE][i],"data");
    }

    printf("KDAR:\n");
    fTrackBlock[kMUO]->Print();

    for (int i=0; i<fNTracks[kMUO]; i++) { 
      PrintTrack(fTrackBlock[kMUO]->Track(i),&fTrackPar[kMUO][i],"data");
    }
  }
//-----------------------------------------------------------------------------
// bit 3: KMUO tracks with large DPF > 5 MeV
//-----------------------------------------------------------------------------
  // TStnTrackBlock* cprb = fTrackBlock[kMUO];
  // int ntrk = cprb->NTracks();
  // for (int itrk=0; itrk<ntrk; itrk++) {
  //   trk = cprb->Track(itrk);
  //   tp  = &fTrackPar[kMUO][itrk];
  //   //    if ((GetDebugBit(3) == 1) && (tp->fDpF > 5.) && (trk->NActive() > 25) && (trk->Chi2Dof() < 4)) {
  //   if ((GetDebugBit(3) == 1) && (tp->fIDWord[fBestID[kMUO]] == 0) && (tp->fDpF > 1.)) {
  //     GetHeaderBlock()->Print(Form("TEmuModule bit003: tp->DpF = %10.3f trk->fP = %10.3f trk->fPFront = %10.3f nactv: %4i chi2d: %10.3f",
  // 				   tp->fDpF, tp->fP,tp->fPFront,trk->NActive(),trk->Chi2Dof()));
  //   }
  //  }
}

//-----------------------------------------------------------------------------
int TEmuModule::EndJob() {
  TAnaModule::EndJob();
  return 0;
}

//_____________________________________________________________________________
void TEmuModule::Test001() {
}


}
