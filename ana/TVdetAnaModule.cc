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
// 4  : events with NHitsTF > 1
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
#include "Stntuple/obj/TDisk.hh"
#include "Stntuple/val/stntuple_val_functions.hh"
//------------------------------------------------------------------------------
// Mu2e offline includes
//-----------------------------------------------------------------------------
// #include "CalorimeterGeom/inc/HexMap.hh"

#include "ana/TVdetAnaModule.hh"

ClassImp(TVdetAnaModule)
//-----------------------------------------------------------------------------
TVdetAnaModule::TVdetAnaModule(const char* name, const char* title):
  TStnModule(name,title)
{
  fPdgCode       = 11;
  fGeneratorCode = 28;
}

//-----------------------------------------------------------------------------
TVdetAnaModule::~TVdetAnaModule() {
}


//-----------------------------------------------------------------------------
void TVdetAnaModule::BookVdetHitHistograms(VdetHitHist_t* Hist, const char* Folder) {
  //     char name [200];
  //     char title[200];
  //-----------------------------------------------------------------------------
  //  
  //-----------------------------------------------------------------------------

  HBook1F(Hist->fIndex  ,"index",Form("%s: VD index",Folder),  1000, 0, 1000,Folder);
  HBook1F(Hist->fPdgCode,"pdg_code",Form("%s: PDG code",Folder),  200,-1000, 1000,Folder);
  HBook1F(Hist->fGenCode,"gen_code",Form("%s: generator code",Folder),  100, -10, 90,Folder);
  HBook1F(Hist->fEnergy ,"energy",Form("%s: Hit Energy" ,Folder),200, 0, 200,Folder);
  HBook1F(Hist->fTime   ,"time"  ,Form("%s: Hit Time  " ,Folder),200, 0,2000,Folder);
}


//-----------------------------------------------------------------------------
void TVdetAnaModule::BookEventHistograms(EventHist_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];

  HBook1F(Hist->fNVdetHits ,"nhvd"   ,Form("%s: N VDET Hits"         ,Folder), 1000, 0,  1000,Folder);
  HBook1F(Hist->fNHitsTF   ,"nhtf"   ,Form("%s: N Hits Tracker Front",Folder), 100, 0,  100,Folder);
  HBook1F(Hist->fNHitsTB   ,"nhtb"   ,Form("%s: N Hits Tracker Back ",Folder), 100, 0,  100,Folder);

  HBook1F(Hist->fMomTF   ,"ptf"   ,Form("%s: Mom  Tracker Front ",Folder), 1000, 0,  200,Folder);
  HBook1F(Hist->fMomTB   ,"ptb"   ,Form("%s: Mom  Tracker Back " ,Folder), 1000, 0,  200,Folder);
  HBook1F(Hist->fMomLoss ,"dpfb"   ,Form("%s: Mom  Loss Front-Back " ,Folder), 1000, 0,  10,Folder);
}

//_____________________________________________________________________________
void TVdetAnaModule::BookHistograms() {

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
  book_event_histset[ 1] = 0;		// 

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
// book virtual detector hit histograms
//-----------------------------------------------------------------------------
  int book_vdethit_histset[kNVdetHitHistSets];
  for (int i=0; i<kNVdetHitHistSets; i++) book_vdethit_histset[i] = 0;

  book_vdethit_histset[0] = 1;		// all clusters

  for (int i=0; i<kNVdetHitHistSets; i++) {
    if (book_vdethit_histset[i] != 0) {
      sprintf(folder_name,"vdt_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fVdetHit[i] = new VdetHitHist_t;
      BookVdetHitHistograms(fHist.fVdetHit[i],Form("Hist/%s",folder_name));
    }
  }
}

//-----------------------------------------------------------------------------
void TVdetAnaModule::FillVdetHitHistograms(VdetHitHist_t* Hist, TVdetHitData* Hit) {

  Hist->fIndex->Fill(Hit->Index());
  Hist->fPdgCode->Fill(Hit->PdgCode());
  Hist->fGenCode->Fill(Hit->GeneratorCode());
  Hist->fEnergy->Fill(Hit->Energy());
  Hist->fTime  ->Fill(Hit->Time  ());
}

//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
void TVdetAnaModule::FillEventHistograms(EventHist_t* Hist) {
//   double            cos_th, xv, yv, rv, zv, p;
//   TLorentzVector    mom;

  float  loss (-1);

  Hist->fNVdetHits->Fill(fNVdetHits);
  Hist->fNHitsTF->Fill(fNHitsTF);
  Hist->fNHitsTB->Fill(fNHitsTB);
  Hist->fMomTF->Fill(fMomTF);
  Hist->fMomTB->Fill(fMomTB);
  if (fMomTF > 0) loss = fMomTF-fMomTB;
  Hist->fMomLoss->Fill(loss);
}

//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TVdetAnaModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
  RegisterDataBlock("VdetBlock" ,"TVdetDataBlock" ,&fVdetDataBlock);
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();

  return 0;
}


//_____________________________________________________________________________
void TVdetAnaModule::FillHistograms() {

//-----------------------------------------------------------------------------
// event histograms
//
// EVT_0: all events
//-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0]);
//-----------------------------------------------------------------------------
// straw hit histograms
//-----------------------------------------------------------------------------
//  int            nh;
  TVdetHitData* hit;

  for (int i=0; i<fNVdetHits; i++) {
    hit = fVdetDataBlock->Hit(i);
    FillVdetHitHistograms(fHist.fVdetHit[0],hit);
  }
//-----------------------------------------------------------------------------
// fill GENP histograms
// GEN_0: all particles
//-----------------------------------------------------------------------------
//   TGenParticle* genp;
//   for (int i=0; i<fNGenp; i++) {
//     genp = fGenpBlock->Particle(i);
//     FillGenpHistograms(fHist.fGenp[0],genp);
//  }
}



//_____________________________________________________________________________
int TVdetAnaModule::BeginRun() {
  int rn = GetHeaderBlock()->RunNumber();
  TStntuple::Init(rn);
  return 0;
}


//_____________________________________________________________________________
int TVdetAnaModule::Event(int ientry) {

  //  double                p;
  //  TLorentzVector        mom;

  fVdetDataBlock->GetEntry(ientry);
  //  fGenpBlock->GetEntry(ientry);
//-----------------------------------------------------------------------------
// assume electron in the first particle, otherwise the logic will need to 
// be changed
// if there are several hits, use the first one
//-----------------------------------------------------------------------------
//  fNGenp      = fGenpBlock->NParticles();
  TVdetHitData* hit;

  fNHitsTF   = 0;
  fNHitsTB   = 0;
  fMomTF     = -1.;
  fMomTB     = -1.;

  fNVdetHits = fVdetDataBlock->NHits();

  for (int i=0; i<fNVdetHits; i++) {
    hit = fVdetDataBlock->Hit(i);
    if (hit->Index() == 13) {
					// tracker FRONT
      fNHitsTF += 1;
      if (fMomTF < 0) fMomTF = hit->McMomentum();
    }
    else if (hit->Index() == 15) {
					// tracker END
      fNHitsTB += 1;
      if (fMomTB < 0) fMomTB = hit->McMomentum();
    }
  }

  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TVdetAnaModule::Debug() {

//-----------------------------------------------------------------------------
// bit 4: events with NHitsTF > 1
//-----------------------------------------------------------------------------
  if (GetDebugBit(4) == 1) {
    if (fNHitsTF > 1) {
      GetHeaderBlock()->Print(Form("NHits(TF) = %5i",fNHitsTF));
    }
  }
}

//_____________________________________________________________________________
int TVdetAnaModule::EndJob() {
  printf("----- end job: ---- %s\n",GetName());
  return 0;
}

//_____________________________________________________________________________
void TVdetAnaModule::Test001() {

  // mu2e::HexMap* hmap      = new mu2e::HexMap();

  // mu2e::HexLK hex_index(0,0);

  // for (int i=0; i<40; i++) {
  //   hex_index = hmap->lk(i);
  //   printf(" i,l,k = %5i %5i %5i\n",i,hex_index._l,hex_index._k);
  // }
}

