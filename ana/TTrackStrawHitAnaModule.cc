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

#include "ana/TTrackStrawHitAnaModule.hh"

ClassImp(TTrackStrawHitAnaModule)
//-----------------------------------------------------------------------------
TTrackStrawHitAnaModule::TTrackStrawHitAnaModule(const char* name, const char* title):
  TStnModule(name,title)
{
  fHitBlockName = "TrackHitBlock";
}

//-----------------------------------------------------------------------------
TTrackStrawHitAnaModule::~TTrackStrawHitAnaModule() {
}


//-----------------------------------------------------------------------------
void TTrackStrawHitAnaModule::BookTrackStrawHitHistograms(HistBase_t* Hist, const char* Folder) {
  //     char name [200];
  //     char title[200];
  TrackStrawHitHist_t* hist = (TrackStrawHitHist_t*) Hist;

  HBook1F(hist->fMcDoca       ,"mcdoca"     ,Form("%s: MC Doca"         ,Folder), 300,-3, 3   ,Folder);
  HBook1F(hist->fEnergy       ,"energy"     ,Form("%s: Hit Charge"      ,Folder), 200, 0, 0.02,Folder);
  HBook1F(hist->fDriftRadius  ,"rdrift"     ,Form("%s: Drift Radius"    ,Folder), 300, 0, 3  ,Folder);
  HBook1F(hist->fMcMomentum   ,"mc_mom"     ,Form("%s: MC Momentum"     ,Folder), 200, 0, 200,Folder);
  HBook1F(hist->fTime         ,"time"       ,Form("%s: Time"            ,Folder), 200, 0, 2000,Folder);

  HBook2F(hist->fDriftRadiusVsMcDoca,"rdrift_vs_mcdoca",Form("%s: R(drift) vs abs(mc doca)",Folder), 
	  300, 0, 3   , 300, 0,3, Folder);

  HBook1F(hist->fDr  ,"dr"     ,Form("%s: Drift - abs(mcdoca)" ,Folder), 250, -1, 4 ,Folder);
}


//-----------------------------------------------------------------------------
void TTrackStrawHitAnaModule::BookEventHistograms(HistBase_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];
  EventHist_t* hist = (EventHist_t*) Hist;

  HBook1F(hist->fEventNumber,"evtnum",Form("%s: Event Number",Folder), 1000, 0,   1.e4,Folder);
  HBook1F(hist->fRunNumber  ,"runnum",Form("%s: Run   Number",Folder), 1000, 0,   1.e6,Folder);
  HBook1F(hist->fNTracks    ,"ntracks" ,Form("%s: N(tracks)",Folder),  10,   0,   10  ,Folder);
  //  HBook1F(hist->fNTracks    ,"ntracks" ,Form("%s: N(tracks)",Folder),  200,  0,  200  ,Folder);
}

//-----------------------------------------------------------------------------
void TTrackStrawHitAnaModule::BookTrackHistograms(HistBase_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];
  TrackHist_t* hist = (TrackHist_t*) Hist;

  HBook1F(hist->fNHits  ,"nhits",Form("%s: NHits"  ,Folder), 100, 0,   100, Folder);
  HBook1F(hist->fNActive,"nactv",Form("%s: NActive",Folder), 100, 0,   100, Folder);
  HBook1F(hist->fN500   ,"n500" ,Form("%s: N500"   ,Folder),  20, 0,    20, Folder);
}

//_____________________________________________________________________________
void TTrackStrawHitAnaModule::BookHistograms() {

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
  book_event_histset[ 1] = 0;		// just an example - this should be the default

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

  book_track_histset[ 0] = 1;		// all tracks

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
// book straw hit histograms
//-----------------------------------------------------------------------------
  int book_trk_strawhit_histset[kNTrkStrawHitHistSets];
  for (int i=0; i<kNTrkStrawHitHistSets; i++) book_trk_strawhit_histset[i] = 0;

  book_trk_strawhit_histset[0] = 1;		// all hits
  book_trk_strawhit_histset[1] = 1;		// active hits
  book_trk_strawhit_histset[2] = 1;		// hits with dr > 0.5mm

  for (int i=0; i<kNTrkStrawHitHistSets; i++) {
    if (book_trk_strawhit_histset[i] != 0) {
      sprintf(folder_name,"tsh_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fTrackStrawHit[i] = new TrackStrawHitHist_t;
      BookTrackStrawHitHistograms(fHist.fTrackStrawHit[i],Form("Hist/%s",folder_name));
    }
  }
}

//-----------------------------------------------------------------------------
void TTrackStrawHitAnaModule::FillTrackStrawHitHistograms(HistBase_t* Hist, TTrackStrawHitData* Hit) {

  TrackStrawHitHist_t* hist = (TrackStrawHitHist_t*) Hist;
  
  float mcdoca = Hit->McDoca();
  float rdrift = Hit->DriftRadius();
  float dr     = rdrift-fabs(mcdoca);

  hist->fMcDoca->Fill(mcdoca);
  hist->fEnergy->Fill(Hit->Energy());
  hist->fDriftRadius->Fill(rdrift);
  hist->fMcMomentum->Fill(Hit->McMomentum());
  hist->fTime->Fill(Hit->Time());
  hist->fDriftRadiusVsMcDoca->Fill(fabs(mcdoca),rdrift);
  hist->fDr->Fill(dr);
}

//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
void TTrackStrawHitAnaModule::FillEventHistograms(HistBase_t* Hist) {
//   double            cos_th, xv, yv, rv, zv, p;
//   TLorentzVector    mom;

  EventHist_t* hist = (EventHist_t*) Hist;

  int event_number = GetHeaderBlock()->EventNumber();
  int run_number   = GetHeaderBlock()->RunNumber();

  hist->fEventNumber->Fill(event_number);
  hist->fRunNumber->Fill(run_number);
  hist->fNTracks->Fill(fNTracks);
}

//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
void TTrackStrawHitAnaModule::FillTrackHistograms(HistBase_t* Hist, TrackParBase_t* Tp) {
//   double            cos_th, xv, yv, rv, zv, p;
//   TLorentzVector    mom;

  TrackHist_t* hist = (TrackHist_t*) Hist;
  TrackPar_t*  tp   = (TrackPar_t* ) Tp;

  hist->fNHits->Fill(tp->fNHits);
  hist->fNActive->Fill(tp->fNActive);
  hist->fN500->Fill(tp->fN500);
}

//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TTrackStrawHitAnaModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
  RegisterDataBlock(fHitBlockName.Data() ,"TTrackStrawDataBlock" ,&fTrackStrawHitDataBlock);
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();

  return 0;
}


//_____________________________________________________________________________
void TTrackStrawHitAnaModule::FillHistograms() {

//-----------------------------------------------------------------------------
// event histograms
//
// EVT_0: all events
//-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0]);
//-----------------------------------------------------------------------------
// track histograms
//-----------------------------------------------------------------------------
  int ntrk = fTrackStrawHitDataBlock->NTracks();

  for (int it=0; it<ntrk; it++) {
    FillTrackHistograms(fHist.fTrack[0],&fTp[it]);
  }

//-----------------------------------------------------------------------------
// straw hit histograms
//-----------------------------------------------------------------------------
  TTrackStrawHitData* hit;
  float               dr;

  for (int it=0; it<ntrk; it++) {
    int nhits = fTrackStrawHitDataBlock->NTrackHits(it);
    int first = fTrackStrawHitDataBlock->First(it);
    for (int i=0; i<nhits; i++) {
      hit = fTrackStrawHitDataBlock->Hit(first+i);
      dr  = hit->DriftRadius()-fabs(hit->McDoca());

      FillTrackStrawHitHistograms(fHist.fTrackStrawHit[0],hit);

      if (hit->fActive) FillTrackStrawHitHistograms(fHist.fTrackStrawHit[1],hit);
      if (dr > 0.5    ) FillTrackStrawHitHistograms(fHist.fTrackStrawHit[2],hit);
    }
  }
}

//_____________________________________________________________________________
int TTrackStrawHitAnaModule::BeginRun() {
  int rn = GetHeaderBlock()->RunNumber();
  TStntuple::Init(rn);
  return 0;
}


//_____________________________________________________________________________
int TTrackStrawHitAnaModule::Event(int ientry) {

  //  double                p;
  //  TLorentzVector        mom;

  fTrackStrawHitDataBlock->GetEntry(ientry);
//-----------------------------------------------------------------------------
// assume electron in the first particle, otherwise the logic will need to 
// be changed
// if there are several hits, use the first one
//-----------------------------------------------------------------------------
  fNTracks    = fTrackStrawHitDataBlock->NTracks();

  for (int it=0; it<fNTracks; it++) {
    int nhits    = fTrackStrawHitDataBlock->NTrackHits(it);
    int first    = fTrackStrawHitDataBlock->First(it);

    fTp[it].fNHits   = nhits;
    fTp[it].fNActive = 0;
    fTp[it].fN500    = 0;

    TTrackStrawHitData* hit;
    float               dr;

    for (int i=0; i<nhits; i++) {
      hit = fTrackStrawHitDataBlock->Hit(first+i);
      if (hit->Active()) fTp[it].fNActive += 1;
      dr = hit->DriftRadius()-fabs(hit->McDoca());
      if (dr > 0.5) fTp[it].fN500 += 1;
    }
  }

  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TTrackStrawHitAnaModule::Debug() {

//-----------------------------------------------------------------------------
// bit 4: events with NProtonStrawHits >= 20
//-----------------------------------------------------------------------------
  if (GetDebugBit(4) == 1) {
    GetHeaderBlock()->Print(Form("NTracks = %5i", fNTracks));
  }
}

//_____________________________________________________________________________
int TTrackStrawHitAnaModule::EndJob() {
  printf("----- end job: ---- %s\n",GetName());
  return 0;
}

//_____________________________________________________________________________
void TTrackStrawHitAnaModule::Test001() {

  // mu2e::HexMap* hmap      = new mu2e::HexMap();

  // mu2e::HexLK hex_index(0,0);

  // for (int i=0; i<40; i++) {
  //   hex_index = hmap->lk(i);
  //   printf(" i,l,k = %5i %5i %5i\n",i,hex_index._l,hex_index._k);
  // }
}

