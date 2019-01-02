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

#include "ana/TStrawHitAnaModule.hh"

ClassImp(TStrawHitAnaModule)
//-----------------------------------------------------------------------------
TStrawHitAnaModule::TStrawHitAnaModule(const char* name, const char* title):
  TStnModule(name,title)
{
  fNStations = 20;
}

//-----------------------------------------------------------------------------
TStrawHitAnaModule::~TStrawHitAnaModule() {
}


//-----------------------------------------------------------------------------
void TStrawHitAnaModule::BookStrawHitHistograms(HistBase_t* Hist, const char* Folder) {
  //     char name [200];
  //     char title[200];
  StrawHitHist_t* hist = (StrawHitHist_t*) Hist;

  HBook1F(hist->fGeneratorCode,"gen_code"   ,Form("%s: Generator code"  ,Folder),  200, 0   ,  200,Folder);
  HBook1F(hist->fPdgCode      ,"pdg_code"   ,Form("%s: PDG code"        ,Folder), 1000,-5000, 5000,Folder);
  HBook1F(hist->fMotherPdgCode,"mother_code",Form("%s: Mother PDG code" ,Folder), 200,-1000, 1000,Folder);
  HBook1F(hist->fMcMomentum   ,"mc_mom"     ,Form("%s: MC Momentum"     ,Folder), 200, 0, 200,Folder);
  HBook1F(hist->fEnergy       ,"energy"     ,Form("%s: Hit Charge"      ,Folder), 200, 0, 0.02,Folder);
  HBook1F(hist->fTime         ,"time"       ,Form("%s: Time"            ,Folder), 200, 0, 2000,Folder);
  HBook1F(hist->fDt           ,"dt"         ,Form("%s: DeltaT = T1-T2"  ,Folder), 200, -10, 10,Folder);
  HBook1F(hist->fStation      ,"station"    ,Form("%s: Station"         ,Folder),  40,   0, 40,Folder);
  HBook1F(hist->fPanel        ,"panel"      ,Form("%s: Panel"           ,Folder),  10,   0, 10,Folder);
  HBook1F(hist->fLayer        ,"layer"      ,Form("%s: Layer"           ,Folder),   5,   0,  5,Folder);
  HBook1F(hist->fStraw        ,"straw"      ,Form("%s: Straw"           ,Folder), 100,   0,100,Folder);
  HBook1F(hist->fPreamp       ,"preamp"     ,Form("%s: Preamp"          ,Folder), 100,   0,100,Folder);
  HBook1F(hist->fFace         ,"face"       ,Form("%s: Face"            ,Folder),  10,   0, 10,Folder);
}


//-----------------------------------------------------------------------------
void TStrawHitAnaModule::BookEventHistograms(HistBase_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];
  EventHist_t* hist = (EventHist_t*) Hist;

  HBook1F(hist->fEventNumber      ,"evtnum"     ,Form("%s: Event Number"           ,Folder), 1000, 0,  1.e4,Folder);
  HBook1F(hist->fRunNumber        ,"runnum"     ,Form("%s: Run   Number"           ,Folder), 1000, 0,  1.e6,Folder);
  HBook1F(hist->fNStrawHits [0]   ,"nsh_0"      ,Form("%s: N(Straw Hits)[0]"       ,Folder),  200, 0,   200,Folder);
  HBook1F(hist->fNStrawHits [1]   ,"nsh_1"      ,Form("%s: N(Straw Hits)[1]"       ,Folder), 1000, 0, 10000,Folder);
  HBook1F(hist->fNStrawHits [2]   ,"nsh_2"      ,Form("%s: N(Straw Hits)[2] T>500" ,Folder), 1000, 0, 10000,Folder);
  HBook1F(hist->fNProtonHits[0]   ,"n_prot_sh_0",Form("%s: N(Proton Straw Hits)[0]",Folder),  200, 0,   200,Folder);
  HBook1F(hist->fNProtonHits[1]   ,"n_prot_sh_1",Form("%s: N(Proton Straw Hits)[1]",Folder), 1000, 0, 10000,Folder);
  HBook1F(hist->fNStationsWithHits,"nst_w_hits" ,Form("%s: N(stations w/hits)"     ,Folder),   50, 0,    50,Folder);
  HBook1F(hist->fDeltaSt          ,"dst"        ,Form("%s: stations: Last-First"   ,Folder),   20, 0,    20,Folder);
}

//_____________________________________________________________________________
void TStrawHitAnaModule::BookHistograms() {

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
// book GENP histograms
//-----------------------------------------------------------------------------
  int book_strawhit_histset[kNStrawHitHistSets];
  for (int i=0; i<kNStrawHitHistSets; i++) book_strawhit_histset[i] = 0;

  book_strawhit_histset[0] = 1;		// all hits
  book_strawhit_histset[1] = 1;		// proton hits
  book_strawhit_histset[2] = 1;		// hits with T > 500 ns

  for (int i=0; i<kNStrawHitHistSets; i++) {
    if (book_strawhit_histset[i] != 0) {
      sprintf(folder_name,"strh_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fStrawHit[i] = new StrawHitHist_t;
      BookStrawHitHistograms(fHist.fStrawHit[i],Form("Hist/%s",folder_name));
    }
  }
}

//-----------------------------------------------------------------------------
void TStrawHitAnaModule::FillStrawHitHistograms(HistBase_t* Hist, TStrawHitData* Hit) {

  StrawHitHist_t* hist = (StrawHitHist_t*) Hist;
  
  hist->fGeneratorCode->Fill(Hit->GeneratorCode());
  hist->fPdgCode->Fill(Hit->PdgCode());
  hist->fMotherPdgCode->Fill(Hit->MotherPdgCode());
  hist->fEnergy->Fill(Hit->Energy());
  hist->fTime->Fill(Hit->Time());
  hist->fDt->Fill(Hit->Dt());
  hist->fStation->Fill(Hit->Station());
  hist->fPanel->Fill(Hit->Panel());
  hist->fStraw->Fill(Hit->Straw());
  hist->fFace->Fill(Hit->Face());
  hist->fLayer->Fill(Hit->Layer());
  hist->fPreamp->Fill(Hit->Preamp());
}

//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
void TStrawHitAnaModule::FillEventHistograms(HistBase_t* Hist) {
//   double            cos_th, xv, yv, rv, zv, p;
//   TLorentzVector    mom;

  EventHist_t* hist = (EventHist_t*) Hist;

  int event_number = GetHeaderBlock()->EventNumber();
  int run_number   = GetHeaderBlock()->RunNumber();

  hist->fEventNumber->Fill(event_number);
  hist->fRunNumber->Fill(run_number);
  hist->fNStrawHits[0]->Fill(fNStrawHits);
  hist->fNStrawHits[1]->Fill(fNStrawHits);
  hist->fNStrawHits[2]->Fill(fNsh500);
  hist->fNProtonHits[0]->Fill(fNProtonStrawHits);
  hist->fNProtonHits[1]->Fill(fNProtonStrawHits);
  hist->fNStationsWithHits->Fill(fNStationsWithHits);
  int dst = fLastStation-fFirstStation;
  hist->fDeltaSt->Fill(dst);
}

//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TStrawHitAnaModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
  RegisterDataBlock("StrawDataBlock" ,"TStrawDataBlock" ,&fStrawHitDataBlock);
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();

  return 0;
}


//_____________________________________________________________________________
void TStrawHitAnaModule::FillHistograms() {

//-----------------------------------------------------------------------------
// event histograms
//
// EVT_0: all events
//-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0]);
//-----------------------------------------------------------------------------
// straw hit histograms
//-----------------------------------------------------------------------------
  TStrawHitData* hit;

  for (int i=0; i<fNStrawHits; i++) {
    hit = fStrawHitDataBlock->Hit(i);
    FillStrawHitHistograms(fHist.fStrawHit[0],hit);
    if (hit->PdgCode() == 2212) FillStrawHitHistograms(fHist.fStrawHit[1],hit);
    if (hit->Time()    >  500.) FillStrawHitHistograms(fHist.fStrawHit[2],hit);
  }
}



//_____________________________________________________________________________
int TStrawHitAnaModule::BeginRun() {
  int rn = GetHeaderBlock()->RunNumber();
  TStntuple::Init(rn);
  return 0;
}


//_____________________________________________________________________________
int TStrawHitAnaModule::Event(int ientry) {

  //  double                p;
  //  TLorentzVector        mom;

  fStrawHitDataBlock->GetEntry(ientry);
//-----------------------------------------------------------------------------
// assume electron in the first particle, otherwise the logic will need to 
// be changed
// if there are several hits, use the first one
//-----------------------------------------------------------------------------
  fNStrawHits       = fStrawHitDataBlock->NHits();
  fNProtonStrawHits = 0;
  fNsh500           = 0;

  fNStationsWithHits = 0;
  fFirstStation      = -1;
  fLastStation       = -1;

  for (int i=0; i<fNStations; i++) fNHitsPerStation[i] = 0;

  TStrawHitData* hit;
  for (int i=0; i<fNStrawHits; i++) {
    hit = fStrawHitDataBlock->Hit(i);
    if (hit->PdgCode() == 2212) fNProtonStrawHits += 1;
    if (hit->Time()    >   500) fNsh500 += 1;

    int ist = hit->Station();
    fNHitsPerStation[ist]++;
  }

  for (int i=0; i<fNStations; i++) {
    if (fNHitsPerStation[i] > 0) {
      fNStationsWithHits++;
      if (fFirstStation < 0) fFirstStation = i;
      if (fLastStation  < i) fLastStation  = i;
    }
  }

  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TStrawHitAnaModule::Debug() {

//-----------------------------------------------------------------------------
// bit 4: events with NProtonStrawHits >= 20
//-----------------------------------------------------------------------------
  if (GetDebugBit(4) == 1) {
    if (fNProtonStrawHits >= 20) {
      GetHeaderBlock()->Print(Form("NStrawHits = %5i, NProtonStrawHits = %5i",
				   fNStrawHits,fNProtonStrawHits));
    }
  }
}

//_____________________________________________________________________________
int TStrawHitAnaModule::EndJob() {
  printf("----- end job: ---- %s\n",GetName());
  return 0;
}

//_____________________________________________________________________________
void TStrawHitAnaModule::Test001() {

  // mu2e::HexMap* hmap      = new mu2e::HexMap();

  // mu2e::HexLK hex_index(0,0);

  // for (int i=0; i<40; i++) {
  //   hex_index = hmap->lk(i);
  //   printf(" i,l,k = %5i %5i %5i\n",i,hex_index._l,hex_index._k);
  // }
}

