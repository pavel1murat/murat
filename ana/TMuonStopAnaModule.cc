///////////////////////////////////////////////////////////////////////////////
// 2017-05-24: analyse STNTUPLE produced by running on the output of Stage3 simulation
//             (output stream selected by tgtStopFilter)
// stopped muon is the last TSimParticle, its momentum in the end shoudl be equal to zero
// 
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
// 3  : events with P(muon gp) > 8000
// 4  : events with electron or muon gparent=NULL
// 5  : events with electron or muon grandparent PDG code > 2500
///////////////////////////////////////////////////////////////////////////////
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TEnv.h"
#include "TSystem.h"

#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/alg/TStntuple.hh"
#include "Stntuple/val/stntuple_val_functions.hh"
//------------------------------------------------------------------------------
// Mu2e offline includes
//-----------------------------------------------------------------------------
// #include "CalorimeterGeom/inc/HexMap.hh"

#include "ana/TMuonStopAnaModule.hh"

ClassImp(TMuonStopAnaModule)
//-----------------------------------------------------------------------------
TMuonStopAnaModule::TMuonStopAnaModule(const char* name, const char* title):
  TStnModule(name,title)
{
}

//-----------------------------------------------------------------------------
TMuonStopAnaModule::~TMuonStopAnaModule() {
}


//-----------------------------------------------------------------------------
void TMuonStopAnaModule::BookSimpHistograms(HistBase_t* Hist, const char* Folder) {
  SimpHist_t* hist = (SimpHist_t*) Hist;

  HBook1F(hist->fVolumeID   ,"vol_id"   ,Form("%s: Volume ID"   ,Folder), 200,  2400, 2600,Folder);
  HBook1F(hist->fGeneratorID,"gen_id"   ,Form("%s: Generator ID",Folder), 200,   -10,  190,Folder);
  HBook1F(hist->fTime       ,"time"     ,Form("%s: Stop Time"   ,Folder), 200,     0, 2000,Folder);
  HBook1F(hist->fParentMom  ,"pmom"     ,Form("%s: Parent Mom"  ,Folder), 200,     0, 2000,Folder);
  HBook1F(hist->fParentPDG  ,"ppdg"     ,Form("%s: Parent PDG"  ,Folder), 200, -1000, 1000,Folder);

  HBook2F(hist->fYVsX       ,"y_vs_x"   ,Form("%s: Y vs X (all)"    ,Folder), 200,  -200, 200, 200, -200, 200, Folder);
  HBook2F(hist->fYVsX_2480  ,"y_vs_x_2480",Form("%s: Y vs X [2480]" ,Folder), 200,  -200, 200, 200, -200, 200, Folder);
  HBook2F(hist->fYVsX_2513  ,"y_vs_x_2513",Form("%s: Y vs X [2513]" ,Folder), 200,  -200, 200, 200, -200, 200, Folder);
}


//-----------------------------------------------------------------------------
void TMuonStopAnaModule::BookVdetHistograms(HistBase_t* Hist, const char* Folder) {

  VdetHist_t* hist = (VdetHist_t*) Hist;

  HBook1F(hist->fIndex   ,"index"   ,Form("%s: VD index"      ,Folder),1000, 0, 1000,Folder);
  HBook1F(hist->fPDGCode ,"pdg_code",Form("%s: PDG code"      ,Folder),2000,-1000, 1000,Folder);
  HBook1F(hist->fGenCode ,"gen_code",Form("%s: generator code",Folder), 100, -10, 90,Folder);
  HBook1F(hist->fMomentum,"mom"     ,Form("%s: Momentum"      ,Folder), 200, 0, 200,Folder);
  HBook1F(hist->fTime    ,"time"    ,Form("%s: Hit Time  "    ,Folder), 200, 0,2000,Folder);
  HBook2F(hist->fYVsX    ,"y_vs_x"  ,Form("%s: Y vs X (all)"  ,Folder), 200,  -1000, 1000, 200, -1000,1000, Folder);
}


//-----------------------------------------------------------------------------
void TMuonStopAnaModule::BookEventHistograms(HistBase_t* Hist, const char* Folder) {
  EventHist_t* hist = (EventHist_t*) Hist;

  HBook1F(hist->fEventNumber ,"evtnum" ,Form("%s: Event Number"  ,Folder), 1000, 0, 1.e4,Folder);
  HBook1F(hist->fRunNumber   ,"runnum" ,Form("%s: Run   Number"  ,Folder), 1000, 0, 1.e6,Folder);
  HBook1F(hist->fNVdetHits   ,"nvdh"   ,Form("%s: N VDET hits"   ,Folder),  100, 0,  100,Folder);
  HBook1F(hist->fNVdetHits_9 ,"nvdh_9" ,Form("%s: N VDET= 9 hits",Folder),  100, 0,  100,Folder);
  HBook1F(hist->fNVdetHits_13,"nvdh_13",Form("%s: N VDET=13 hits",Folder),  100, 0,  100,Folder);
  HBook1F(hist->fETot_13     ,"etot_13",Form("%s: PTot VDET=13"  ,Folder),  200, 0,  100,Folder);
}

//_____________________________________________________________________________
void TMuonStopAnaModule::BookHistograms() {

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
// book SimParticle histograms
//-----------------------------------------------------------------------------
  int book_simp_histset[kNSimpHistSets];
  for (int i=0; i<kNSimpHistSets; i++) book_simp_histset[i] = 0;

  book_simp_histset[0] = 1;		// all steps

  for (int i=0; i<kNSimpHistSets; i++) {
    if (book_simp_histset[i] != 0) {
      sprintf(folder_name,"simp_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fSimp[i] = new SimpHist_t;
      BookSimpHistograms(fHist.fSimp[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book Vdet histograms - so far need only VDET=9 (in front of the stopping target)
//-----------------------------------------------------------------------------
  int book_vdet_histset[kNVdetHistSets];
  for (int i=0; i<kNVdetHistSets; i++) book_vdet_histset[i] = 0;

  book_vdet_histset[  9] = 1;		// all particles, VDET= 9 , ST_In
  book_vdet_histset[ 13] = 1;		// all particles, VDET=13 , TT_FrontHollow
  book_vdet_histset[ 14] = 1;		// all particles, VDET=13 , TT_FrontHollow, r > 40

  book_vdet_histset[109] = 1;		// e-  , VDET=9
  book_vdet_histset[209] = 1;		// e+  , VDET=9
  book_vdet_histset[309] = 1;		// mu- , VDET=9
  book_vdet_histset[409] = 1;		// mu+ , VDET=9

  for (int i=0; i<kNVdetHistSets; i++) {
    if (book_vdet_histset[i] != 0) {
      sprintf(folder_name,"vdet_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fVdet[i] = new VdetHist_t;
      BookVdetHistograms(fHist.fVdet[i],Form("Hist/%s",folder_name));
    }
  }
}

//-----------------------------------------------------------------------------
void TMuonStopAnaModule::FillSimpHistograms(HistBase_t* Hist, TSimParticle* Simp, SimpData_t* Sd) {

  SimpHist_t* hist = (SimpHist_t*) Hist;
  
  hist->fVolumeID->Fill(Simp->fEndVolumeIndex);
  hist->fGeneratorID->Fill(Simp->fGeneratorID);
  
  float x = Simp->EndPos()->X()+3904.;
  float y = Simp->EndPos()->Y();
  float t = Simp->EndPos()->T();

  hist->fTime->Fill(t);
  hist->fParentMom->Fill(fParent->StartMom()->P());
  hist->fParentPDG->Fill(fParent->PDGCode());

  hist->fYVsX->Fill(x,y);

  if (Simp->fEndVolumeIndex == 2480) hist->fYVsX_2480->Fill(x,y);
  if (Simp->fEndVolumeIndex == 2513) hist->fYVsX_2513->Fill(x,y);
}

//-----------------------------------------------------------------------------
void TMuonStopAnaModule::FillVdetHistograms(HistBase_t* Hist, TStepPointMC* Step) {

  VdetHist_t* hist = (VdetHist_t*) Hist;
  
  hist->fIndex   ->Fill(Step->VolumeID());
  hist->fPDGCode ->Fill(Step->PDGCode());
  hist->fGenCode ->Fill(Step->GenIndex());
  hist->fMomentum->Fill(Step->Mom()->Mag());
  hist->fTime    ->Fill(Step->Time());
  hist->fYVsX    ->Fill(Step->Pos()->X()+3904.,Step->Pos()->Y());
}

//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
void TMuonStopAnaModule::FillEventHistograms(HistBase_t* Hist) {

  EventHist_t* hist = (EventHist_t*) Hist;

  int event_number = GetHeaderBlock()->EventNumber();
  int run_number   = GetHeaderBlock()->RunNumber();

  hist->fEventNumber->Fill(event_number);
  hist->fRunNumber->Fill(run_number);
  hist->fNVdetHits->Fill(fNVdetHits);
  hist->fNVdetHits_9 ->Fill(fNVdetHits_9);
  hist->fNVdetHits_13->Fill(fNVdetHits_13);
  hist->fETot_13     ->Fill(fETot_13);

}

//-----------------------------------------------------------------------------
void TMuonStopAnaModule::FillHistograms() {
//-----------------------------------------------------------------------------
// event histograms
//
// EVT_0: all events
//-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0]);
//-----------------------------------------------------------------------------
// SimParticle histograms: there should be only one stopped muon per event
//-----------------------------------------------------------------------------
  SimpData_t*   sd =  fSimpData;

  FillSimpHistograms(fHist.fSimp[0],fMuon,sd);
//-----------------------------------------------------------------------------
// VDET histograms
//-----------------------------------------------------------------------------
  for (int i=0; i<fNVdetHits; i++) {
    TStepPointMC* step = fVdetBlock->StepPointMC(i);

    if (step->VolumeID() ==  9) FillVdetHistograms(fHist.fVdet[ 9],step);
    if (step->VolumeID() == 13) {
      FillVdetHistograms(fHist.fVdet[13],step);
      float x = step->Pos()->X()+3904.;
      float y = step->Pos()->Y();
      float r = sqrt(x*x+y*y);
      if ((r >= 400) && (r < 800)) { 
	FillVdetHistograms(fHist.fVdet[14],step);
      }
      
    }

    if (step->PDGCode() == 11) {
      if (step->VolumeID() == 9) FillVdetHistograms(fHist.fVdet[109],step);
    }
    if (step->PDGCode() == -11) {
      if (step->VolumeID() == 9) FillVdetHistograms(fHist.fVdet[209],step);
    }
    if (step->PDGCode() == 13) {
      if (step->VolumeID() == 9) FillVdetHistograms(fHist.fVdet[309],step);
    }
    if (step->PDGCode() == -13) {
      if (step->VolumeID() == 9) FillVdetHistograms(fHist.fVdet[409],step);
    }
  }
}



//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TMuonStopAnaModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
//  RegisterDataBlock("SpmcBlock","TStepPointMCBlock",&fStepPointMCBlock);
  RegisterDataBlock("SimpBlock","TSimpBlock"   ,&fSimpBlock       );
  RegisterDataBlock("VdetBlock","TSpmcBlock"   ,&fVdetBlock       );
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();

  return 0;
}


//_____________________________________________________________________________
int TMuonStopAnaModule::BeginRun() {
  int rn = GetHeaderBlock()->RunNumber();
  TStntuple::Init(rn);
  return 0;
}


//-----------------------------------------------------------------------------
// if erething went well, the stopped muon should be the last particle
//-----------------------------------------------------------------------------
int TMuonStopAnaModule::Event(int ientry) {

  //  double                p;
  //  TLorentzVector        mom;
  //  fStepPointMCBlock->GetEntry(ientry);

  fSimpBlock->GetEntry(ientry);
  fVdetBlock->GetEntry(ientry);

  fNVdetHits    = fVdetBlock->NStepPoints();
  fNVdetHits_9  = 0;
  fNVdetHits_13 = 0;
  fETot_13      = 0;

  for (int i=0; i<fNVdetHits; i++) {
    TStepPointMC* step = fVdetBlock->StepPointMC(i);
    if (step->VolumeID() ==  9) fNVdetHits_9 += 1;
    if (step->VolumeID() == 13) {
      float x = step->Pos()->X()+3904.;
      float y = step->Pos()->Y();
      float r = sqrt(x*x+y*y);
      if ((r >= 400) && (r < 800)) {
	fNVdetHits_13 += 1;
	fETot_13      += step->Mom()->Mag();
      }
    }
  }

  int np = fSimpBlock->NParticles();
  fProton = fSimpBlock->Particle(0);
  fMuon  = fSimpBlock->Particle(np-1);

  TSimParticle* parent = fMuon;
//-----------------------------------------------------------------------------
// loop over "steps" in SimpBlock
//-----------------------------------------------------------------------------
  int loc = np-1;
  while (loc >= 0) {
    int parent_id = parent->ParentID();
    parent = fSimpBlock->FindParticle(parent_id);
//-----------------------------------------------------------------------------
// sometimes history tree includes a scattered proton, which produces a pion,
// decaying into a muon. In this case call pion a 'grandparent'
//-----------------------------------------------------------------------------
    if (parent && (parent->PDGCode() != 13)) break;
  }

  fParent = parent;

  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TMuonStopAnaModule::Debug() {

//-----------------------------------------------------------------------------
// bit 4: events with NHitsTF > 1
//-----------------------------------------------------------------------------
  // if (GetDebugBit(4) == 1) {
  //   GetHeaderBlock()->Print(Form("NHits(TF) = %5i",fNGenp));
  // }
}

//_____________________________________________________________________________
int TMuonStopAnaModule::EndJob() {
  printf("----- end job: ---- %s\n",GetName());
  return 0;
}

//_____________________________________________________________________________
void TMuonStopAnaModule::Test001() {
}

