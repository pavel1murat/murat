///////////////////////////////////////////////////////////////////////////////
// use of tmp:
//
// Tmp(0) : nax seg
// Tmp(1) : nst seg
// 
// use of debug bits: bits 0-2 are reserved
// bit 0 : all events
// bit 1 : passed events
// bit 2 : rejected events
// 
// bit 3 : UNUSED
// bit 4 : events with NProtonStrawHits >= 20
// bit 5 : hits 1120 < T < 1140: search for repetitions
// bit 6 : hits 1315 < SppTime < 1360: search for repetitions
// bit 7 : print hits with SppTime > 1 ms
///////////////////////////////////////////////////////////////////////////////
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TEnv.h"
#include "TSystem.h"

#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/alg/TStntuple.hh"
#include "Stntuple/obj/TStrawHit.hh"
#include "Stntuple/val/stntuple_val_functions.hh"

#include "ana/TStrawHitAnaModule.hh"

ClassImp(TStrawHitAnaModule)
//-----------------------------------------------------------------------------
TStrawHitAnaModule::TStrawHitAnaModule(const char* name, const char* title): TStnModule(name,title)
{
  fNStations = 20;
}

//-----------------------------------------------------------------------------
TStrawHitAnaModule::~TStrawHitAnaModule() {
}


//-----------------------------------------------------------------------------
void TStrawHitAnaModule::BookEventHistograms(HistBase_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];
  EventHist_t* hist = (EventHist_t*) Hist;

  HBook1F(hist->fEventNumber      ,"evtnum"     ,Form("%s: Event Number"           ,Folder), 1000, 0,  1.e4,Folder);
  HBook1F(hist->fRunNumber        ,"runnum"     ,Form("%s: Run   Number"           ,Folder), 1000, 1e5, 1.e5+1000,Folder);
  HBook1F(hist->fInstLum          ,"inst_lum"   ,Form("%s: Inst Lum"               ,Folder), 1000, 0,  20e7,Folder);
  HBook1F(hist->fNStrawHits [0]   ,"nsh_0"      ,Form("%s: N(Straw Hits)[0]"       ,Folder),  200, 0,   200,Folder);
  HBook1F(hist->fNStrawHits [1]   ,"nsh_1"      ,Form("%s: N(Straw Hits)[1]"       ,Folder), 1000, 0, 10000,Folder);
  HBook1F(hist->fNStrawHits [2]   ,"nsh_2"      ,Form("%s: N(Straw Hits)[2] T>200" ,Folder), 1000, 0, 10000,Folder);
  HBook1F(hist->fNStrawHits [3]   ,"nsh_3"      ,Form("%s: N(Straw Hits)[2] T>500" ,Folder), 1000, 0, 10000,Folder);
  HBook1F(hist->fNProtonHits[0]   ,"n_prot_sh_0",Form("%s: N(Proton Straw Hits)[0]",Folder),  200, 0,   200,Folder);
  HBook1F(hist->fNProtonHits[1]   ,"n_prot_sh_1",Form("%s: N(Proton Straw Hits)[1]",Folder), 1000, 0, 10000,Folder);
  HBook1F(hist->fNStationsWithHits,"nst_w_hits" ,Form("%s: N(stations w/hits)"     ,Folder),   50, 0,    50,Folder);
  HBook1F(hist->fDeltaSt          ,"dst"        ,Form("%s: stations: Last-First"   ,Folder),   20, 0,    20,Folder);
}

//-----------------------------------------------------------------------------
void TStrawHitAnaModule::BookStrawHitHistograms(HistBase_t* Hist, const char* Folder) {
  //     char name [200];
  //     char title[200];
  StrawHitHist_t* hist = (StrawHitHist_t*) Hist;

  HBook1F(hist->fGeneratorCode,"gen_id"     ,Form("%s: Generator ID"         ,Folder),  200,    0,   200,Folder);
  HBook1F(hist->fSimID        ,"sim_id"     ,Form("%s: SimP ID"              ,Folder), 1000,    0,  1000,Folder);
  HBook1F(hist->fPdgCode      ,"pdg_id"     ,Form("%s: PDG ID"               ,Folder), 1000,-5000,  5000,Folder);
  HBook1F(hist->fMotherPdgCode,"mother_id"  ,Form("%s: Mother PDG ID"        ,Folder),  200,-1000,  1000,Folder);
  HBook1F(hist->fMcMomentum   ,"mc_mom"     ,Form("%s: MC Momentum"          ,Folder),  200,    0,   200,Folder);
  HBook1F(hist->fEnergy       ,"energy"     ,Form("%s: Hit Charge"           ,Folder),  200,    0,  0.02,Folder);
  HBook1F(hist->fTime         ,"time"       ,Form("%s: Time"                 ,Folder), 1000,    0, 10000,Folder);
  HBook1F(hist->fSppTime[0]   ,"spptime_0"  ,Form("%s: Simp prod Time[0], ns",Folder),  200,    0,  2000,Folder);
  HBook1F(hist->fSppTime[1]   ,"spptime_1"  ,Form("%s: Simp prod Time[1], us",Folder), 5000,    0, 10000,Folder);
  HBook1F(hist->fSppTime[2]   ,"spptime_2"  ,Form("%s: Simp prod Time[2], ms",Folder),  200,    0,    20,Folder);
  HBook1F(hist->fDt           ,"dt"         ,Form("%s: DeltaT = T1-T2"       ,Folder),  200, -100,   100,Folder);
  HBook1F(hist->fStation      ,"station"    ,Form("%s: Station"              ,Folder),   40,    0,    40,Folder);
  HBook1F(hist->fPanel        ,"panel"      ,Form("%s: Panel"                ,Folder),   10,    0,    10,Folder);
  HBook1F(hist->fLayer        ,"layer"      ,Form("%s: Layer"                ,Folder),    5,    0,     5,Folder);
  HBook1F(hist->fStraw        ,"straw"      ,Form("%s: Straw"                ,Folder),  100,    0,   100,Folder);
  HBook1F(hist->fPreamp       ,"preamp"     ,Form("%s: Preamp"               ,Folder),  100,    0,   100,Folder);
  HBook1F(hist->fFace         ,"face"       ,Form("%s: Face"                 ,Folder),   10,    0,    10,Folder);
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
  TString*  hit_selection  [kNStrawHitHistSets];
  for (int i=0; i<kNStrawHitHistSets; i++) hit_selection[i] = nullptr;

  hit_selection[  0] = new TString("all      hits");
  hit_selection[  1] = new TString("all      hits in events with N>=20 hits");
  hit_selection[ 10] = new TString("e+/-     hits");
  hit_selection[ 11] = new TString("e+/-     hits, P <  20 MeV/c");
  hit_selection[ 12] = new TString("e+/-     hits, P >= 20 MeV/c");
  hit_selection[ 20] = new TString("mu+/-    hits");
  hit_selection[ 30] = new TString("pi+/-    hits");
  hit_selection[ 40] = new TString("proton   hits");
  hit_selection[200] = new TString("all hits T>200ns");
  hit_selection[500] = new TString("all hits T>500ns");

  for (int i=0; i<kNStrawHitHistSets; i++) {
    if (hit_selection[i] != 0) {
      sprintf(folder_name,"strh_%i",i);
      fol          = (TFolder*) hist_folder->FindObject(folder_name);
      const char* title = hit_selection[i]->Data();
      if (! fol) fol = hist_folder->AddFolder(folder_name,title);
      fHist.fStrawHit[i] = new StrawHitHist_t;
      BookStrawHitHistograms(fHist.fStrawHit[i],Form("Hist/%s",folder_name));
    }
  }
}

//-----------------------------------------------------------------------------
void TStrawHitAnaModule::FillStrawHitHistograms(HistBase_t* Hist, TStrawHit* Hit, StrawHitPar_t* Shp) {

  StrawHitHist_t* hist = (StrawHitHist_t*) Hist;
  
  hist->fGeneratorCode->Fill(Hit->GenID());
  hist->fSimID->Fill(Hit->SimID());
  hist->fPdgCode->Fill(Hit->PdgID());
  hist->fMotherPdgCode->Fill(Hit->MotherPdgID());
  hist->fEnergy->Fill(Hit->EDep());
  hist->fMcMomentum->Fill(Hit->McMom());
  hist->fTime->Fill(Hit->Time(0));
  hist->fSppTime[0]->Fill(Shp->fSppTime);               // ns
  hist->fSppTime[1]->Fill(Shp->fSppTime/1.e3);          // us
  hist->fSppTime[2]->Fill(Shp->fSppTime/1.e6);          // ms
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

  int   event_number = GetHeaderBlock()->EventNumber();
  int   run_number   = GetHeaderBlock()->RunNumber();
  float inst_lum     = GetHeaderBlock()->InstLum();

  hist->fEventNumber->Fill(event_number);
  hist->fRunNumber->Fill(run_number);
  hist->fInstLum->Fill(inst_lum);
  hist->fNStrawHits[0]->Fill(fNStrawHits);
  hist->fNStrawHits[1]->Fill(fNStrawHits);
  hist->fNStrawHits[2]->Fill(fNsh200);
  hist->fNStrawHits[3]->Fill(fNsh500);
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
  RegisterDataBlock("SimpBlock"     ,"TSimpDataBlock" ,&fSimpBlock);
  RegisterDataBlock("StrawHitBlock" ,"TStrawHitBlock" ,&fStrawHitBlock);
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
  for (int i=0; i<fNStrawHits; i++) {
    if (i > kMaxNStrawHits) {
      GetHeaderBlock()->Print(Form("ERROR: hit number %i greater than kMaxNStrawHits = %i, SKIP.",i,kMaxNStrawHits));
      break;
    }

    TStrawHit* hit = fStrawHitBlock->Hit(i);
    StrawHitPar_t* shp = fStrawHitPar+i;
    int pdg_id = hit->PdgID();

    FillStrawHitHistograms(fHist.fStrawHit[0],hit,shp);

    if (fNStrawHits >= 20) FillStrawHitHistograms(fHist.fStrawHit[1],hit,shp);

    if      (abs(pdg_id) ==    11) { 
      FillStrawHitHistograms(fHist.fStrawHit[ 10],hit,shp);
      if (hit->McMom() < 20) FillStrawHitHistograms(fHist.fStrawHit[ 11],hit,shp);
      else                   FillStrawHitHistograms(fHist.fStrawHit[ 12],hit,shp);
    }
    else if (abs(pdg_id) ==    13) FillStrawHitHistograms(fHist.fStrawHit[ 20],hit,shp);
    else if (abs(pdg_id) ==   211) FillStrawHitHistograms(fHist.fStrawHit[ 30],hit,shp);
    else if (abs(pdg_id) ==  2212) FillStrawHitHistograms(fHist.fStrawHit[ 40],hit,shp);

    if      (hit->Time(0) >  200.) FillStrawHitHistograms(fHist.fStrawHit[200],hit,shp);
    if      (hit->Time(0) >  500.) FillStrawHitHistograms(fHist.fStrawHit[500],hit,shp);

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

  fStrawHitBlock->GetEntry(ientry);
  if (fSimpBlock) fSimpBlock->GetEntry(ientry);
//-----------------------------------------------------------------------------
// assume electron in the first particle, otherwise the logic will need to 
// be changed
// if there are several hits, use the first one
//-----------------------------------------------------------------------------
  fNStrawHits        = fStrawHitBlock->NHits();
  fNProtonStrawHits  = 0;
  fNsh200            = 0;
  fNsh500            = 0;

  fNStationsWithHits = 0;
  fFirstStation      = -1;
  fLastStation       = -1;

  for (int i=0; i<fNStations; i++) fNHitsPerStation[i] = 0;

  TStrawHit* hit;
  for (int i=0; i<fNStrawHits; i++) {
    if (i > kMaxNStrawHits) {
      GetHeaderBlock()->Print(Form("ERROR: hit number %i greater than kMaxNStrawHits = %i, skip.",i,kMaxNStrawHits));
      break;
    }

    StrawHitPar_t* shp = fStrawHitPar+i;

    hit = fStrawHitBlock->Hit(i);
    if (hit->PdgID() == 2212) fNProtonStrawHits += 1;
    if (hit->Time(0) >   200) fNsh200 += 1;
    if (hit->Time(0) >   500) fNsh500 += 1;

    int ist = hit->Station();
    fNHitsPerStation[ist]++;

    int sim_id = hit->SimID();
    shp->fSppTime = -1;
//-----------------------------------------------------------------------------
// MC information may be missing
//-----------------------------------------------------------------------------
    if (fSimpBlock and (fSimpBlock->NParticles() > 0)) { 
      TSimParticle* simp(nullptr);
    
      simp = fSimpBlock->FindParticle(sim_id);
      if (simp) {
	shp->fSppTime = simp->StartPos()->T();
      }
      else {
	GetHeaderBlock()->Print(Form(" ERROR in TStrawHitAnaModule::Event: simp=NULL for hit=%i and sim_id=%i",i,sim_id));
      }
    }
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
//-----------------------------------------------------------------------------  
// bit 5 : print hits with 1120 < T < 1140
//-----------------------------------------------------------------------------
  if (GetDebugBit(5) == 1) {
    int banner_printed(0);
    for (int i=0; i<fNStrawHits; i++) {
      StrawHitPar_t* shp = fStrawHitPar+i;
      TStrawHit* hit = fStrawHitBlock->Hit(i);
      if ((hit->Time(0) > 1120) && (hit->Time(0) < 1140)) {
	if (banner_printed == 0) { 
	  GetHeaderBlock()->Print("bit 005: Event of interest");
	  PrintStrawHit(hit,shp,"banner") ; 
	  banner_printed = 1; 
	}
	PrintStrawHit(hit,shp,"data") ;
      }
    }
  }
//-----------------------------------------------------------------------------  
// bit 6 : print hits with 1315 < SppTime < 1360
//-----------------------------------------------------------------------------
  if (GetDebugBit(6) == 1) {
    int banner_printed(0);
    for (int i=0; i<fNStrawHits; i++) {
      TStrawHit* hit = fStrawHitBlock->Hit(i);
      StrawHitPar_t* shp = fStrawHitPar+i;
      if ((shp->fSppTime > 1315) && (shp->fSppTime <= 1360)) {
	if (banner_printed == 0) { 
	  GetHeaderBlock()->Print("bit 006: Event of interest");
	  PrintStrawHit(hit,shp,"banner") ; 
	  banner_printed = 1; 
	}
	PrintStrawHit(hit,shp,"data") ;
      }
    }
  }
//-----------------------------------------------------------------------------  
// bit 7 : print hits with SppTime > 1 ms
//-----------------------------------------------------------------------------
  if (GetDebugBit(7) == 1) {
    int nhits = 0;
    for (int i=0; i<fNStrawHits; i++) {
      // TStrawHitData* hit = fStrawHitDataBlock->Hit(i);
      StrawHitPar_t* shp = fStrawHitPar+i;
      if (shp->fSppTime > 1e6) {
	nhits += 1;
      }
    }
    if (nhits > 0) GetHeaderBlock()->Print(Form("bit 007: hits with spptime > 1e6"));
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

//_____________________________________________________________________________
void TStrawHitAnaModule::PrintStrawHit(TStrawHit* Hit, StrawHitPar_t* Shp, const char* Option) {

  TString opt = Option;
  opt.ToLower();

  if ((opt == "") || (opt.Index("banner") >= 0)) {
    printf("----------------------------------------------------------------------------------------------------------------------------------\n");
    printf("  Index    Time       Dt       Energy    PdgCode    PdgCode(M)  GenCode  SimID    McMom   parentID pGenID    pMom    Shp->fSppTime\n");
    printf("----------------------------------------------------------------------------------------------------------------------------------\n");
  }

  if (opt == "banner") return;

  float pmom     = -1;
  int   pgenID   = -1;
  int   parentID = -1;
  
  if (fSimpBlock and (fSimpBlock->NParticles() > 0)) {
    int simID = Hit->SimID();
    TSimParticle* sim = fSimpBlock->FindParticle(simID);

    parentID = sim->ParentID();

    TSimParticle* parent = sim;
    
    while (parent != nullptr) {
      parent   = fSimpBlock->FindParticle(parentID);
      if (parent) {
	pmom     = parent->StartMom()->P();
	pgenID   = parent->GeneratorID();
	parentID = parent->ParentID();
      }
    }
  }

  printf("%6i %9.3f %9.3f %5.1f %5.1f %10.3f %10i %10i %8i %8i %9.3f %10i %5i %9.3f %12.3e\n",
	 Hit->StrawID(), 
	 Hit->Time(0),
	 Hit->Time(1),
	 Hit->TOT (0),
	 Hit->TOT (1),
	 Hit->EDep(),
	 Hit->GenID(),
	 Hit->SimID(),
	 Hit->PdgID(),
	 Hit->MotherPdgID(),
	 Hit->McMom(),
	 parentID, 
	 pgenID  , 
	 pmom    ,
	 Shp->fSppTime
	 );
  
}

