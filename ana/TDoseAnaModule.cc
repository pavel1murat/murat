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
#include "Stntuple/val/stntuple_val_functions.hh"
//------------------------------------------------------------------------------
// Mu2e offline includes
//-----------------------------------------------------------------------------
#include "ana/TDoseAnaModule.hh"

ClassImp(TDoseAnaModule)
//-----------------------------------------------------------------------------
TDoseAnaModule::TDoseAnaModule(const char* name, const char* title):
  TStnModule(name,title)
{
  fPdgCode       = 11;
  fGeneratorCode = 28;
  fSpmcBlockName = "SpmcBlock" ; // default for Stntuple
}

//-----------------------------------------------------------------------------
TDoseAnaModule::~TDoseAnaModule() {
}


//-----------------------------------------------------------------------------
void TDoseAnaModule::BookSpmcHistograms(HistBase_t* Hist, const char* Folder) {
  //     char name [200];
  //     char title[200];
  SpmcHist_t* hist = (SpmcHist_t*) Hist;

  HBook1F(hist->fVolumeID,"vol_id"   ,Form("%s: VolumeID"       ,Folder), 200,    0, 10000,Folder);
  HBook1F(hist->fGenIndex,"gen_index",Form("%s: GenIndex"       ,Folder), 200,    0, 10000,Folder);
  HBook1F(hist->fSimID   ,"sim_id"   ,Form("%s: SimID"          ,Folder), 200,    0,  1000,Folder);

  HBook1F(hist->fPDGCode[0] ,"pdg_0" ,Form("%s: PDG code"       ,Folder),  700,  -350,   350,Folder);
  HBook1F(hist->fPDGCode[1] ,"pdg_1" ,Form("%s: PDG code"       ,Folder),  700, -3500,  3500,Folder);

  HBook1F(hist->fCreationCode   ,"cr_code" ,Form("%s: Creation code",Folder), 200,   0,   200,Folder);
  HBook1F(hist->fParentSimID    ,"psim_id" ,Form("%s: Parent SimID",Folder), 200,   0,  1000,Folder);
  HBook1F(hist->fParentPDGCode  ,"ppdg"    ,Form("%s: Parent PDG code",Folder), 700, -350,   350,Folder);
  HBook1F(hist->fEndProcessCode ,"end_code",Form("%s: End process code",Folder), 200,   0,   200,Folder);
  HBook1F(hist->fEDepTot        ,"edep_tot",Form("%s: EDEP tot"        ,Folder), 200,   0,   10 ,Folder);
  HBook1F(hist->fEDepNio        ,"edep_nio",Form("%s: EDEP NIO"        ,Folder), 200,   0,   10 ,Folder);
  HBook1F(hist->fTime           ,"time"    ,Form("%s: Time"            ,Folder), 200,   0,  2000,Folder);
  HBook1F(hist->fStepLength     ,"step"    ,Form("%s: Ltep Length"     ,Folder), 200,   0,   100,Folder);
  HBook1F(hist->fMomentum       ,"mom"     ,Form("%s: Momentum"        ,Folder), 500,   0,   250,Folder);
  HBook1F(hist->fCosTh          ,"cos_th"  ,Form("%s: Cos_TH"          ,Folder), 200,  -1,   1  ,Folder);

  HBook2F(hist->fYVsX           ,"y_vs_x"     ,Form("%s: Y vs X"       ,Folder), 500,  -1000, 1000, 500, -500,  500, Folder);
  HBook2F(hist->fYVsZ           ,"y_vs_z"     ,Form("%s: Y vs Z"       ,Folder), 500,  -500,   500, 500,  -500, 500, Folder);
}


//-----------------------------------------------------------------------------
void TDoseAnaModule::BookEventHistograms(HistBase_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];

  //  char name [200];
  //  char title[200];
  EventHist_t* hist = (EventHist_t*) Hist;

  HBook1F(hist->fEventNumber,"evtnum",Form("%s: Event Number",Folder), 1000, 0,  1.e4,Folder);
  HBook1F(hist->fRunNumber  ,"runnum",Form("%s: Run   Number",Folder), 1000, 0,  1.e6,Folder);
  HBook1F(hist->fNSteps     ,"nsteps",Form("%s: NSteps"      ,Folder), 1000, 0,  1.e3,Folder);
}

//_____________________________________________________________________________
void TDoseAnaModule::BookHistograms() {

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
// book StepPointMC histograms
//-----------------------------------------------------------------------------
  int book_spmc_histset[kNSpmcHistSets];
  for (int i=0; i<kNSpmcHistSets; i++) book_spmc_histset[i] = 0;

  book_spmc_histset[0] = 1;		// all clusters

  for (int i=0; i<kNSpmcHistSets; i++) {
    if (book_spmc_histset[i] != 0) {
      sprintf(folder_name,"spmc_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fSpmc[i] = new SpmcHist_t;
      BookSpmcHistograms(fHist.fSpmc[i],Form("Hist/%s",folder_name));
    }
  }
}

// //-----------------------------------------------------------------------------
// // 
// //-----------------------------------------------------------------------------
// void TDoseAnaModule::FillEventHistograms(HistBase_t* Hist) {
// //   double            cos_th, xv, yv, rv, zv, p;
// //   TLorentzVector    mom;
// }

//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TDoseAnaModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
  RegisterDataBlock(fSpmcBlockName.Data(),"TStepPointMCBlock" ,&fSpmcBlock);
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();

  return 0;
}


//-----------------------------------------------------------------------------
void TDoseAnaModule::FillSpmcHistograms(HistBase_t* Hist, TStepPointMC* Step, SpmcData_t* Sd) {

  SpmcHist_t* hist = (SpmcHist_t*) Hist;

  int vol_id = Step->VolumeID();
  
  hist->fVolumeID->Fill(vol_id);
  hist->fGenIndex->Fill(Step->GenIndex());
  hist->fSimID   ->Fill(Step->SimID());
  hist->fPDGCode[0] ->Fill(Step->PDGCode());
  hist->fPDGCode[1] ->Fill(Step->PDGCode());
  hist->fCreationCode->Fill(Step->CreationCode());
  hist->fParentSimID ->Fill(Step->ParentSimID());
  hist->fParentPDGCode->Fill(Step->ParentPDGCode());
  hist->fEndProcessCode->Fill(Step->EndProcessCode());

  hist->fEDepTot->Fill(Step->EDepTot());
  hist->fEDepNio->Fill(Step->EDepNio());
  hist->fTime   ->Fill(Step->Time());
  hist->fStepLength->Fill(Step->StepLength());

  double p = Step->Mom()->Mag();
  hist->fMomentum->Fill(p);

  float x = Step->Pos()->X();
  float y = Step->Pos()->Y();
  float z = Step->Pos()->Z();

  hist->fYVsX->Fill(x,y);		// useful for stage 2
  hist->fYVsZ->Fill(z,y);		// useful for stage 1

}

//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
void TDoseAnaModule::FillEventHistograms(HistBase_t* Hist) {
//   double            cos_th, xv, yv, rv, zv, p;
//   TLorentzVector    mom;

  EventHist_t* hist = (EventHist_t*) Hist;

  int event_number = GetHeaderBlock()->EventNumber();
  int run_number   = GetHeaderBlock()->RunNumber();

  hist->fEventNumber->Fill(event_number);
  hist->fRunNumber->Fill(run_number);

  hist->fNSteps->Fill(fNSteps);

}

//_____________________________________________________________________________
void TDoseAnaModule::FillHistograms() {

//-----------------------------------------------------------------------------
// event histograms
//
// EVT_0: all events
//-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0]);
//-----------------------------------------------------------------------------
// StepPointMC histograms
//-----------------------------------------------------------------------------
  int nsteps = fSpmcBlock->NStepPoints();

  for (int i=0; i<nsteps; i++) {
    TStepPointMC* s  = fSpmcBlock->StepPointMC(i);
    SpmcData_t*   sd = fSpmcData+i;

    FillSpmcHistograms(fHist.fSpmc[0],s,sd);
  }
}



//_____________________________________________________________________________
int TDoseAnaModule::BeginRun() {
  int rn = GetHeaderBlock()->RunNumber();
  TStntuple::Init(rn);
  return 0;
}


//_____________________________________________________________________________
int TDoseAnaModule::Event(int ientry) {

  //  double                p;
  //  TLorentzVector        mom;

  fSpmcBlock->GetEntry(ientry);
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
//  TStepPointMC* spmc;

  int nsteps = fSpmcBlock->NStepPoints();

  for (int i=0; i<nsteps; i++) {
    //    TStepPointMC* s       = fSpmcBlock->StepPointMC(i);
    SpmcData_t* spmc_data = fSpmcData+i;
    //    int id                = s->GetUniqueID();
//-----------------------------------------------------------------------------
// find particle in SimpBlock
//-----------------------------------------------------------------------------
    spmc_data->fParticle  = NULL; // fSimpBlock->FindParticle(id);
    spmc_data->fParent    = NULL;
    spmc_data->fGParent   = NULL;
  }

  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TDoseAnaModule::Debug() {

// //-----------------------------------------------------------------------------
// // bit 4: events with NHitsTF > 1
// //-----------------------------------------------------------------------------
//   if (GetDebugBit(4) == 1) {
//     if (fNHitsTF > 1) {
//       GetHeaderBlock()->Print(Form("NHits(TF) = %5i",fNHitsTF));
//     }
//   }
}

//_____________________________________________________________________________
int TDoseAnaModule::EndJob() {
  printf("----- end job: ---- %s\n",GetName());
  return 0;
}

//_____________________________________________________________________________
void TDoseAnaModule::Test001() {
}

