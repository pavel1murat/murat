///////////////////////////////////////////////////////////////////////////////
// use of tmp:
//
// Tmp(0) : nax seg
// Tmp(1) : nst seg
// 
// use of debug bits: bits 0-2 are reserved
// 0  : all events
// 1  : passed events
// 2  : rejected events
// 
// 3  : UNUSED
// 4  : events with NProtonStrawHits >= 20
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
#include "Stntuple/obj/TStrawHitData.hh"
#include "Stntuple/val/stntuple_val_functions.hh"
//------------------------------------------------------------------------------
// Mu2e offline includes
//-----------------------------------------------------------------------------
// #include "CalorimeterGeom/inc/HexMap.hh"

#include "ana/TTriggerAnaModule.hh"

ClassImp(murat::TTriggerAnaModule)

namespace murat {
//-----------------------------------------------------------------------------
TTriggerAnaModule::TTriggerAnaModule(const char* name, const char* title): TAnaModule(name,title)
{
  fPdgCode        = 11;
  fGeneratorCode  = -1;
  fMinTrigMom     = 0.;
  fTrackBlockName = "TrackBlockPar";
}

//-----------------------------------------------------------------------------
TTriggerAnaModule::~TTriggerAnaModule() {
}


//-----------------------------------------------------------------------------
void TTriggerAnaModule::BookTimeClusterHistograms(HistBase_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];
  TimeClusterHist_t* hist = (TimeClusterHist_t*) Hist;

  HBook1F(hist->fEnergy    ,"e"     ,Form("%s: cluster energy",Folder), 200,    0,  200,Folder);
  HBook1F(hist->fNHits     ,"nhits" ,Form("%s: nhits"         ,Folder), 200,    0,  200,Folder);
  HBook1F(hist->fR         ,"r"     ,Form("%s: radius"        ,Folder), 200,    0, 1000,Folder);
  HBook1F(hist->fTime      ,"time"  ,Form("%s: time"          ,Folder), 200,    0, 2000,Folder);
  
}

//-----------------------------------------------------------------------------
void TTriggerAnaModule::BookHelixHistograms(HistBase_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];
  //  TrackSeedHist_t* hist = (TrackSeedHist_t*) Hist;
}

//-----------------------------------------------------------------------------
void TTriggerAnaModule::BookTrackSeedHistograms(HistBase_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];
  TrackSeedHist_t* hist = (TrackSeedHist_t*) Hist;

  HBook1F(hist->fP      ,"p"    ,Form("%s: mom"    ,Folder), 1100,    0,  110,Folder);
  HBook1F(hist->fNHits  ,"nhits",Form("%s: N(hits)",Folder),  200,    0,  200,Folder);
  HBook1F(hist->fChi2Dof,"chi2d",Form("%s: Chi2/DOF" ,Folder),  200,    0,   40,Folder);
  HBook1F(hist->fD0     ,"d0"   ,Form("%s: D0"       ,Folder),  200, -200,  200,Folder);
}

// //-----------------------------------------------------------------------------
// void TTriggerAnaModule::BookTrackHistograms(HistBase_t* Hist, const char* Folder) {
//   //  char name [200];
//   //  char title[200];
//   TrackHist_t* hist = (TrackHist_t*) Hist;

//   HBook1F(hist->fP      ,"p"    ,Form("%s: momentum" ,Folder), 1100,    0,  110,Folder);
//   HBook1F(hist->fNActive,"nactv",Form("%s: N(Active)",Folder),  200,    0,  200,Folder);
//   HBook1F(hist->fChi2Dof,"chi2d",Form("%s: Chi2/DOF" ,Folder),  200,    0,   40,Folder);
//   HBook1F(hist->fT0     ,"t0"   ,Form("%s: T0"       ,Folder), 2000,    0, 2000,Folder);
//   HBook1F(hist->fD0     ,"d0"   ,Form("%s: D0"       ,Folder),  200, -200,  200,Folder);
//   HBook1F(hist->fZ0     ,"z0"   ,Form("%s: z0"       ,Folder),  200,-2000, 2000,Folder);
//   HBook1F(hist->fTanDip ,"tdip" ,Form("%s: TanDip"   ,Folder),  250, -2.5,  2.5,Folder);
//   HBook1F(hist->fAlgMask,"amask",Form("%s: AlgMask"  ,Folder),   10,    0,   10,Folder);
// }

//-----------------------------------------------------------------------------
void TTriggerAnaModule::BookTriggerHistograms(HistBase_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];
  TriggerHist_t* hist = (TriggerHist_t*) Hist;

  HBook1F(hist->fBits     ,"bits"    ,Form("%s: fired trigegr bits" ,Folder), 50,    0,  50,Folder);
}

// //-----------------------------------------------------------------------------
// void TTriggerAnaModule::BookEventHistograms(HistBase_t* Hist, const char* Folder) {
//   //  char name [200];
//   //  char title[200];
//   EventHist_t* hist = (EventHist_t*) Hist;

//   HBook1F(hist->fEventNumber,"evtnum",Form("%s: Event Number"  ,Folder), 1000, 0,  1.e4,Folder);
//   HBook1F(hist->fRunNumber  ,"runnum",Form("%s: Run   Number"  ,Folder), 1000, 0,  1.e6,Folder);
//   HBook1F(hist->fPassed     ,"passed",Form("%s: event passed"  ,Folder),   10, 0,    10,Folder);
//   HBook1F(hist->fNTimeClusters ,"ntcl",Form("%s: N(time clusters)",Folder),   20, 0,    20,Folder);
//   HBook1F(hist->fNHelices ,"nhel",Form("%s: N(helices)",Folder),   10, 0,    10,Folder);
//   HBook1F(hist->fNTrackSeeds[0],"nts_0",Form("%s: N(track seeds)[0]",Folder),   10, 0,    10,Folder);
//   HBook1F(hist->fNTrackSeeds[1],"nts_1",Form("%s: N(track seeds)[1]",Folder),   10, 0,    10,Folder);
//   HBook1F(hist->fNGoodSeeds ,"ngts"  ,Form("%s: N(good seeds)" ,Folder),   10, 0,    10,Folder);
//   HBook1F(hist->fNTracks    ,"ntrk"  ,Form("%s: N(tracks)"     ,Folder),   10, 0,    10,Folder);
//   HBook1F(hist->fMcMom     ,"mc_mom"   ,Form("%s: MC Particle Momentum"            ,Folder),1000,  0,200,Folder);
//   HBook1F(hist->fMcCosTh   ,"mc_costh" ,Form("%s: MC Particle Cos(Theta) Lab"      ,Folder),100,-1,1,Folder);
// }

//_____________________________________________________________________________
void TTriggerAnaModule::BookHistograms() {

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
  book_event_histset[ 1] = 1;		// events with N(seeds) > 0
  book_event_histset[ 2] = 1;		// events with N(good seeds) > 0

  book_event_histset[11] = 1;		// events with N(good seeds) > 0, DIO-weighted

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
// book track seed histograms
//-----------------------------------------------------------------------------
  int book_tpeak_histset[kNTimeClusterHistSets];
  for (int i=0; i<kNTimeClusterHistSets; i++) book_tpeak_histset[i] = 0;

  book_tpeak_histset[ 0] = 1;		// all time clusters

  for (int i=0; i<kNTimeClusterHistSets; i++) {
    if (book_tpeak_histset[i] != 0) {
      sprintf(folder_name,"tpeak_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fTimeCluster[i] = new TimeClusterHist_t;
      BookTimeClusterHistograms(fHist.fTimeCluster[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book track seed histograms
//-----------------------------------------------------------------------------
  int book_tseed_histset[kNTrackSeedHistSets];
  for (int i=0; i<kNTrackSeedHistSets; i++) book_tseed_histset[i] = 0;

  book_tseed_histset[ 0] = 1;		// all track seeds

  book_tseed_histset[ 1] = 1;		// track seeds   |d0| < 200
  book_tseed_histset[ 2] = 1;		// track seeds   |d0| < 200 , P > 80

  book_tseed_histset[11] = 1;		// all track seeds, DIO-weighted

  for (int i=0; i<kNTrackSeedHistSets; i++) {
    if (book_tseed_histset[i] != 0) {
      sprintf(folder_name,"tseed_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fTrackSeed[i] = new TrackSeedHist_t;
      BookTrackSeedHistograms(fHist.fTrackSeed[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book track histograms
//-----------------------------------------------------------------------------
  int book_track_histset[kNTrackHistSets];
  for (int i=0; i<kNTrackHistSets; i++) book_track_histset[i] = 0;

  book_track_histset[  0] = 1;		// all    events
  book_track_histset[  1] = 1;		// tracks |D0| < 200

  book_track_histset[100] = 1;		// events with fNGoodSeeds > 0
  book_track_histset[101] = 1;		// NGoodSeeds > 0, tracks |D0| < 200

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
// book trigger histograms
//-----------------------------------------------------------------------------
  int book_trigger_histset[kNTriggerHistSets];
  for (int i=0; i<kNTriggerHistSets; i++) book_trigger_histset[i] = 0;

  book_trigger_histset[  0] = 1;		// all    events

  for (int i=0; i<kNTriggerHistSets; i++) {
    if (book_trigger_histset[i] != 0) {
      sprintf(folder_name,"trig_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fTrigger[i] = new TriggerHist_t;
      BookTriggerHistograms(fHist.fTrigger[i],Form("Hist/%s",folder_name));
    }
  }

}

//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
// void TTriggerAnaModule::FillEventHistograms(HistBase_t* Hist, double Weight) {
// //   double            cos_th, xv, yv, rv, zv, p;
// //   TLorentzVector    mom;

//   EventHist_t* hist = (EventHist_t*) Hist;

//   int event_number = GetHeaderBlock()->EventNumber();
//   int run_number   = GetHeaderBlock()->RunNumber();

//   hist->fEventNumber->Fill(event_number,Weight);
//   hist->fRunNumber->Fill(run_number,Weight);
//   hist->fPassed->Fill(fPassed,Weight);
//   hist->fNTimeClusters->Fill(fNTimeClusters,Weight);
//   hist->fNHelices->Fill(fNHelices,Weight);
//   hist->fNTrackSeeds[0]->Fill(fNTrackSeeds[0],Weight);
//   hist->fNTrackSeeds[1]->Fill(fNTrackSeeds[1],Weight);
//   hist->fNGoodSeeds->Fill(fNGoodSeeds,Weight);
//   hist->fNTracks->Fill(fNTracks,Weight);
//   hist->fMcMom->Fill(fMcMom,Weight);
//   hist->fMcCosTh->Fill(fMcCosTh,Weight);
// }

//-----------------------------------------------------------------------------
void TTriggerAnaModule::FillHelixHistograms(HistBase_t* Hist, TStnHelix* TPeak) {
  //  HelixHist_t* hist = (HelixHist_t*) Hist;
}

//-----------------------------------------------------------------------------
void TTriggerAnaModule::FillTimeClusterHistograms(HistBase_t* Hist, TStnTimeCluster* TPeak) {
  TimeClusterHist_t* hist = (TimeClusterHist_t*) Hist;

  hist->fNHits->Fill(TPeak->fNHits);
  hist->fEnergy->Fill(TPeak->fClusterEnergy);
  hist->fTime->Fill(TPeak->fClusterTime);

  double x = TPeak->fClusterX;
  double y = TPeak->fClusterY;
  double r = sqrt(x*x+y*y);
  hist->fR->Fill(r);
}

//-----------------------------------------------------------------------------
void TTriggerAnaModule::FillTrackSeedHistograms(HistBase_t* Hist, TStnTrackSeed* Seed, double Weight) {
  TrackSeedHist_t* hist = (TrackSeedHist_t*) Hist;
  hist->fP->Fill (Seed->fP,Weight);		// corrected momentum in the first point
  hist->fNHits->Fill(Seed->fNHits,Weight);
  hist->fChi2Dof->Fill(Seed->fChi2/(Seed->fNHits-5.),Weight);
  hist->fD0->Fill(Seed->fD0,Weight);
}

//-----------------------------------------------------------------------------
// for DIO : ultimately, one would need to renormalize the distribution
//-----------------------------------------------------------------------------
// void TTriggerAnaModule::FillTrackHistograms(HistBase_t* Hist, TStnTrack* Track) {

//   TLorentzVector  mom;
//   int             itrk;
//   TrackPar_t*     tp;

//   TrackHist_t* hist = (TrackHist_t*) Hist;
// 					// pointer to local track parameters
//   itrk = Track->Number();
//   tp   = fTrackPar+itrk;

//   float na = Track->NActive();

//   hist->fP->Fill (tp->fP);		// corrected momentum in the first point
//   hist->fNActive->Fill(na);
//   hist->fChi2Dof->Fill(Track->fChi2/(na-5.));
//   hist->fT0->Fill(Track->fT0);
//   hist->fD0->Fill(Track->fD0);
//   hist->fZ0->Fill(Track->fZ0);
//   hist->fTanDip->Fill(Track->fTanDip);
//   hist->fAlgMask->Fill(Track->AlgMask());
// }

//-----------------------------------------------------------------------------
void TTriggerAnaModule::FillTriggerHistograms(HistBase_t* Hist) {
  TriggerHist_t* hist = (TriggerHist_t*) Hist;

  int nb = fTriggerBlock->Paths()->GetNBits();

  for (int i=0; i<nb; i++) {
    int passed = fTriggerBlock->PathPassed(i);
    if (passed) hist->fBits->Fill(i);
  }
}

//_____________________________________________________________________________
void TTriggerAnaModule::FillHistograms() {

//-----------------------------------------------------------------------------
// time peak histograms
//
// TPEAK_0: all track seeds
//-----------------------------------------------------------------------------
  for (int i=0; i<fNTimeClusters; i++) {
    TStnTimeCluster* tpeak = fTimeClusterBlock->TimeCluster(i);
    FillTimeClusterHistograms(fHist.fTimeCluster[0],tpeak);
  }
//-----------------------------------------------------------------------------
// track seed histograms
//
// TSEED_0 : all track seeds
// TSEED_1 : track seeds  |D0| < 200
// TSEED_11: track seeds  |D0| < 200, DIO weight
//-----------------------------------------------------------------------------
  for (int i=0; i<fNTrackSeeds[0]; i++) {
    TStnTrackSeed* tseed = fTrackSeedBlock->TrackSeed(i);
    FillTrackSeedHistograms(fHist.fTrackSeed[0],tseed);
    if (fabs(tseed->fD0) < 200.) {
       FillTrackSeedHistograms(fHist.fTrackSeed[1],tseed);
       if (tseed->fP > fMinTrigMom) {
	 FillTrackSeedHistograms(fHist.fTrackSeed[2],tseed);
       }
       FillTrackSeedHistograms(fHist.fTrackSeed[11],tseed,fWeight);
    }
  }
//-----------------------------------------------------------------------------
// track histograms
//
// TRK_0: all tracks
// TRK_0: tracks    |D0| < 200 mm
//-----------------------------------------------------------------------------
  for (int i=0; i<fNTracks; i++) {
    TStnTrack* trk = fTrackBlock->Track(i);
    TrackPar_t* tp = fTrackPar+i;
    FillTrackHistograms(fHist.fTrack[0],trk,tp,&fSimPar);
    if (fabs(trk->fD0) < 200.) {
      FillTrackHistograms(fHist.fTrack[1],trk,tp,&fSimPar);
    }
  }
//-----------------------------------------------------------------------------
// TRK_1xx : fNGoodSeeds > 0 (passed events) 
//
// TRK_0: all tracks
// TRK_0: tracks    |D0| < 200 mm
//-----------------------------------------------------------------------------
  if (fNGoodSeeds > 0) {
    for (int i=0; i<fNTracks; i++) {
      TStnTrack* trk = fTrackBlock->Track(i);
      TrackPar_t* tp = fTrackPar+i;
      FillTrackHistograms(fHist.fTrack[100],trk,tp,&fSimPar);
      if (fabs(trk->fD0) < 200.) {
	FillTrackHistograms(fHist.fTrack[101],trk,tp,&fSimPar);
      }
    }
  }
//-----------------------------------------------------------------------------
// event histograms
//
// EVT_0: all events
//-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0],&fEvtPar);
  if (fNTrackSeeds[0] > 0) FillEventHistograms(fHist.fEvent[1],&fEvtPar);
  if (fNGoodSeeds     > 0) FillEventHistograms(fHist.fEvent[2],&fEvtPar);

  if (fNTrackSeeds[0] > 0) FillEventHistograms(fHist.fEvent[11],&fEvtPar);
//-----------------------------------------------------------------------------
// trigger histograms
//
// TRIG_0: all events
//-----------------------------------------------------------------------------
  FillTriggerHistograms(fHist.fTrigger[0]);
}

//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TTriggerAnaModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
  RegisterDataBlock("GenpBlock"           , "TGenpBlock"          , &fGenpBlock       );
  RegisterDataBlock("TimeClusterBlock"    , "TStnTimeClusterBlock", &fTimeClusterBlock);
  RegisterDataBlock("HelixBlock"          , "TStnHelixBlock"      , &fHelixBlock      );
  RegisterDataBlock("SimpBlock"           , "TSimpBlock"          , &fSimpBlock       );
  RegisterDataBlock(fTrackBlockName.Data(), "TStnTrackBlock"      , &fTrackBlock      );
  RegisterDataBlock("TrackSeedBlock"      , "TStnTrackSeedBlock"  , &fTrackSeedBlock  );
  RegisterDataBlock("ClusterBlock"        , "TStnClusterBlock"    , &fClusterBlock    );
  RegisterDataBlock("TriggerBlock"        , "TStnTriggerBlock"    , &fTriggerBlock    );
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();

  return 0;
}

//_____________________________________________________________________________
int TTriggerAnaModule::BeginRun() {
  int rn = GetHeaderBlock()->RunNumber();
  TStntuple::Init(rn);
  return 0;
}

//_____________________________________________________________________________
int TTriggerAnaModule::Event(int ientry) {

  //  double                p;
  //  TLorentzVector        mom;

  fGenpBlock->GetEntry(ientry);
  fSimpBlock->GetEntry(ientry);
  fTimeClusterBlock->GetEntry(ientry);
  fHelixBlock->GetEntry(ientry);
  fTrackBlock->GetEntry(ientry);
  fTrackSeedBlock->GetEntry(ientry);
  fClusterBlock->GetEntry(ientry);
  fTriggerBlock->GetEntry(ientry);
//-----------------------------------------------------------------------------
// MC generator info
//-----------------------------------------------------------------------------
  TGenParticle* genp;
  int           pdg_code, generator_code;

  TStntuple* stnt = TStntuple::Instance();

//-----------------------------------------------------------------------------
// event parameters - standardized
//-----------------------------------------------------------------------------
  fEvtPar.fNCrvClusters     = -1;
  fEvtPar.fNCrvPulses       = -1;
  fEvtPar.fNCrvCoincidences = -1;
  fEvtPar.fNGenp            = fGenpBlock->NParticles();
  fEvtPar.fParticle         = nullptr;

  for (int i=fEvtPar.fNGenp-1; i>=0; i--) {
    genp           = fGenpBlock->Particle(i);
    pdg_code       = genp->GetPdgCode();
    generator_code = genp->GetStatusCode();
    if ((abs(pdg_code) == fPdgCode) && (generator_code == fGeneratorCode)) {
      fEvtPar.fParticle = genp;
      break;
    }
  }

  fGeneratorCode = -1;
  fWeight        =  1;

  TLorentzVector    mom (1.,0.,0.,0);

  if (fEvtPar.fParticle) {
    if ((fEvtPar.fParticle->GetStatusCode() == 7) && (fEvtPar.fParticle->GetPdgCode() == 11)) {
//-----------------------------------------------------------------------------
// flat electrons
//-----------------------------------------------------------------------------
      fGeneratorCode = 7;

      fEvtPar.fParticle->Momentum(mom);
      double e     = mom.E();
      fWeight      = stnt->DioWeightAl(e);
    }
  }

  fMcMom       = mom.P();
  fMcCosTh     = mom.Pz()/fMcMom;
//-----------------------------------------------------------------------------
// may want to revisit the definition of fSimPar in future
//-----------------------------------------------------------------------------
  fSimPar.fParticle = fSimpBlock->Particle(0);
  fSimPar.fTFront   = NULL;
  fSimPar.fTMid     = NULL;
  fSimPar.fTBack    = NULL;
  fSimPar.fGenp     = NULL;
//-----------------------------------------------------------------------------
// consider an event passed, if there is a track seed with loosely defined quality
//-----------------------------------------------------------------------------
  fPassed        = 0;
  fNTimeClusters = fTimeClusterBlock->NTimeClusters(); // all
  fNHelices      = fHelixBlock->NHelices(); // all
  fNTracks       = fTrackBlock->NTracks();
  fNGoodSeeds    = 0;

  fNTrackSeeds[0] = fTrackSeedBlock->NTrackSeeds(); // all
  fNTrackSeeds[1] = 0;

  for (int i=0; i<fNTrackSeeds[0]; i++) {
    TStnTrackSeed* tseed = fTrackSeedBlock->TrackSeed(i);
    if (fabs(tseed->fD0) < 200.) {
      fNTrackSeeds[1] += 1;
      if (tseed->fP > fMinTrigMom) {
	fNGoodSeeds += 1;
      }
    }
  }

  if (fNGoodSeeds > 0) fPassed = 1;

  fNGoodTracks = 0;
  for (int i=0; i<fNTracks; i++) {
    TStnTrack* trk = fTrackBlock->Track(i);
    TrackPar_t* tp = fTrackPar+i;

    tp->fTrackID[0] = TAnaModule::fTrackID_BOX;
    tp->fTrackID[1] = TAnaModule::fTrackID_MVA;
    tp->fLogLH      = TAnaModule::fLogLH;
					// straw man definition of a good track, for trigger efficiency
    if (fabs(trk->fD0) < 200.) {
      if (trk->fP > 100) {
	if ((trk->Chi2Dof() < 5) && (trk->NActive() >= 20)) {
	  fNGoodTracks += 1;
	}
      }
    }
  }

  InitTrackPar(fTrackBlock,fClusterBlock,fTrackPar,&fSimPar);

  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TTriggerAnaModule::Debug() {

//-----------------------------------------------------------------------------
// bit 4: events with NProtonStrawHits >= 20
//-----------------------------------------------------------------------------
  // if (GetDebugBit(4) == 1) {
  //   if (fNProtonStrawHits >= 20) {
  //     GetHeaderBlock()->Print(Form("NStrawHits = %5i, NProtonStrawHits = %5i",
  // 				   fNStrawHits,fNProtonStrawHits));
  //   }
  // }
}

//_____________________________________________________________________________
int TTriggerAnaModule::EndJob() {
  printf("----- end job: ---- %s\n",GetName());
  return 0;
}

//_____________________________________________________________________________
void TTriggerAnaModule::Test001() {

  // mu2e::HexMap* hmap      = new mu2e::HexMap();

  // mu2e::HexLK hex_index(0,0);

  // for (int i=0; i<40; i++) {
  //   hex_index = hmap->lk(i);
  //   printf(" i,l,k = %5i %5i %5i\n",i,hex_index._l,hex_index._k);
  // }
}

}
