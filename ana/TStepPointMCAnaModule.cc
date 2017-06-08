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
// #include "CalorimeterGeom/inc/HexMap.hh"

#include "ana/TStepPointMCAnaModule.hh"

ClassImp(TStepPointMCAnaModule)
//-----------------------------------------------------------------------------
TStepPointMCAnaModule::TStepPointMCAnaModule(const char* name, const char* title):
  TStnModule(name,title)
{
  // fPdgCode       = 11;
  // fGeneratorCode = 28;
  fSpmcBlockName = "SpmcBlock";
}

//-----------------------------------------------------------------------------
TStepPointMCAnaModule::~TStepPointMCAnaModule() {
}


//-----------------------------------------------------------------------------
void TStepPointMCAnaModule::BookStepPointMCHistograms(HistBase_t* Hist, const char* Folder) {
  //     char name [200];
  //     char title[200];
  StepPointMCHist_t* hist = (StepPointMCHist_t*) Hist;

  HBook1F(hist->fVolumeID,"vol_id"   ,Form("%s: VolumeID"       ,Folder), 200,    0, 10000,Folder);
  HBook1F(hist->fGenIndex,"gen_index",Form("%s: GenIndex"       ,Folder), 200,    0, 10000,Folder);
  HBook1F(hist->fSimID   ,"sim_id"   ,Form("%s: SimID"          ,Folder), 200,    0,  1000,Folder);

  HBook1F(hist->fPDGCode[0] ,"pdg_0" ,Form("%s: PDG code"       ,Folder),  500,  -250,   250,Folder);
  HBook1F(hist->fPDGCode[1] ,"pdg_1" ,Form("%s: PDG code"       ,Folder), 2000, -10000, 10000,Folder);

  HBook1F(hist->fCreationCode   ,"cr_code",Form("%s: Creation code",Folder), 200,   0,   200,Folder);
  HBook1F(hist->fParentSimID    ,"psim_id",Form("%s: Parent SimID",Folder), 200,   0,  1000,Folder);
  HBook1F(hist->fParentPDGCode  ,"ppdg"   ,Form("%s: Parent PDG code",Folder), 200, -100,   100,Folder);
  HBook1F(hist->fEndProcessCode ,"end_code",Form("%s: End process code",Folder), 200,   0,   200,Folder);
  HBook1F(hist->fEDepTot        ,"edep_tot",Form("%s: EDEP tot"        ,Folder), 200,   0,   10 ,Folder);
  HBook1F(hist->fEDepNio        ,"edep_nio",Form("%s: EDEP NIO"        ,Folder), 200,   0,   10 ,Folder);
  HBook1F(hist->fTime           ,"time"    ,Form("%s: Time"            ,Folder), 200,   0,  2000,Folder);
  HBook1F(hist->fStepLength     ,"step"    ,Form("%s: Ltep Length"     ,Folder), 200,   0,   100,Folder);
  HBook1F(hist->fMomentum       ,"mom"     ,Form("%s: Momentum"        ,Folder), 500,   0,   250,Folder);

  HBook2F(hist->fYVsX           ,"y_vs_x"     ,Form("%s: Y vs X"       ,Folder), 100, -250,  250, 100, -250, 250, Folder);
  HBook2F(hist->fYVsZ           ,"y_vs_z"     ,Form("%s: Y vs Z"       ,Folder), 500, -250,  250, 500, -250, 250, Folder);

}


//-----------------------------------------------------------------------------
void TStepPointMCAnaModule::BookEventHistograms(HistBase_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];
  EventHist_t* hist = (EventHist_t*) Hist;

  HBook1F(hist->fEventNumber,"evtnum",Form("%s: Event Number",Folder), 1000, 0,  1.e4,Folder);
  HBook1F(hist->fRunNumber  ,"runnum",Form("%s: Run   Number",Folder), 1000, 0,  1.e6,Folder);
}

//_____________________________________________________________________________
void TStepPointMCAnaModule::BookHistograms() {

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
// book StepPointMC histograms
//-----------------------------------------------------------------------------
  int book_spmc_histset[kNStepPointMCHistSets];
  for (int i=0; i<kNStepPointMCHistSets; i++) book_spmc_histset[i] = 0;

  book_spmc_histset[0] = 1;		// all steps
  book_spmc_histset[1] = 1;		// electrons
  book_spmc_histset[2] = 1;		// positrons
  book_spmc_histset[3] = 1;		// mu-
  book_spmc_histset[4] = 1;		// mu+
  book_spmc_histset[5] = 1;		// photons
  book_spmc_histset[6] = 1;		// negative pions
  book_spmc_histset[7] = 1;		// positive pions
  book_spmc_histset[8] = 1;		// protons+antiprotons
  book_spmc_histset[9] = 1;		// everything else

  book_spmc_histset[101] = 1;		// electrons with P > 60 MeV/c
  book_spmc_histset[102] = 1;		// electrons with T > 400 ns
  book_spmc_histset[103] = 1;		// electrons with p < 2 MeV/c
  book_spmc_histset[104] = 1;		// electrons with 2 < p < 3 MeV/c
  book_spmc_histset[105] = 1;		// electrons with p > 3 MeV/c

  book_spmc_histset[203] = 1;		// positrons with p < 2 MeV/c
  book_spmc_histset[204] = 1;		// positrons with 2 < p < 3 MeV/c
  book_spmc_histset[205] = 1;		// positrons with p > 3 MeV/c

  for (int i=0; i<kNStepPointMCHistSets; i++) {
    if (book_spmc_histset[i] != 0) {
      sprintf(folder_name,"spmc_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fStepPointMC[i] = new StepPointMCHist_t;
      BookStepPointMCHistograms(fHist.fStepPointMC[i],Form("Hist/%s",folder_name));
    }
  }
}

//-----------------------------------------------------------------------------
void TStepPointMCAnaModule::FillStepPointMCHistograms(HistBase_t* Hist, TStepPointMC* Step) {

  StepPointMCHist_t* hist = (StepPointMCHist_t*) Hist;
  
  hist->fVolumeID->Fill(Step->VolumeID());
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
void TStepPointMCAnaModule::FillEventHistograms(HistBase_t* Hist) {
//   double            cos_th, xv, yv, rv, zv, p;
//   TLorentzVector    mom;

  EventHist_t* hist = (EventHist_t*) Hist;

  int event_number = GetHeaderBlock()->EventNumber();
  int run_number   = GetHeaderBlock()->RunNumber();

  hist->fEventNumber->Fill(event_number);
  hist->fRunNumber->Fill(run_number);

}

//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TStepPointMCAnaModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks 'SpmcBlock' or 'StepPointMCBlock' (old)
//-----------------------------------------------------------------------------
  RegisterDataBlock(fSpmcBlockName.Data(),"TStepPointMCBlock" ,&fStepPointMCBlock);
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();

  return 0;
}


//_____________________________________________________________________________
void TStepPointMCAnaModule::FillHistograms() {

//-----------------------------------------------------------------------------
// event histograms
//
// EVT_0: all events
//-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0]);
//-----------------------------------------------------------------------------
// straw hit histograms
//-----------------------------------------------------------------------------
  TStepPointMC* spmc;

  int nsteps = fStepPointMCBlock->NStepPoints();
  for (int i=0; i<nsteps; i++) {
    spmc    = fStepPointMCBlock->StepPointMC(i);
    float p = spmc->Mom()->Mag();
    float t = spmc->Time();

    FillStepPointMCHistograms(fHist.fStepPointMC[0],spmc);

    if      (spmc->PDGCode() ==   11) {
      FillStepPointMCHistograms(fHist.fStepPointMC[1],spmc);
      if (p >  60) FillStepPointMCHistograms(fHist.fStepPointMC[101],spmc);
      if (t > 400) FillStepPointMCHistograms(fHist.fStepPointMC[102],spmc);
      if (p <   2)              FillStepPointMCHistograms(fHist.fStepPointMC[103],spmc);
      if ((p >=  2) && (p < 3)) FillStepPointMCHistograms(fHist.fStepPointMC[104],spmc);
      if (p >   3)              FillStepPointMCHistograms(fHist.fStepPointMC[105],spmc);
    }
    else if (spmc->PDGCode() ==  -11) {
      FillStepPointMCHistograms(fHist.fStepPointMC[2],spmc);
      if (p <   2)              FillStepPointMCHistograms(fHist.fStepPointMC[203],spmc);
      if ((p >=  2) && (p < 3)) FillStepPointMCHistograms(fHist.fStepPointMC[204],spmc);
      if (p >   3)              FillStepPointMCHistograms(fHist.fStepPointMC[205],spmc);
    }
    else if (spmc->PDGCode() ==   13) FillStepPointMCHistograms(fHist.fStepPointMC[3],spmc);
    else if (spmc->PDGCode() ==  -13) FillStepPointMCHistograms(fHist.fStepPointMC[4],spmc);
    else if (spmc->PDGCode() ==   22) FillStepPointMCHistograms(fHist.fStepPointMC[5],spmc);
    else if (spmc->PDGCode() == -211) FillStepPointMCHistograms(fHist.fStepPointMC[6],spmc);
    else if (spmc->PDGCode() ==  211) FillStepPointMCHistograms(fHist.fStepPointMC[7],spmc);
    else if (abs(spmc->PDGCode() == 2212)) FillStepPointMCHistograms(fHist.fStepPointMC[8],spmc);
    else                              FillStepPointMCHistograms(fHist.fStepPointMC[9],spmc);

  }
}



//_____________________________________________________________________________
int TStepPointMCAnaModule::BeginRun() {
  int rn = GetHeaderBlock()->RunNumber();
  TStntuple::Init(rn);
  return 0;
}


//_____________________________________________________________________________
int TStepPointMCAnaModule::Event(int ientry) {

  //  double                p;
  //  TLorentzVector        mom;

  fStepPointMCBlock->GetEntry(ientry);
//-----------------------------------------------------------------------------
// assume electron in the first particle, otherwise the logic will need to 
// be changed
// if there are several hits, use the first one
//-----------------------------------------------------------------------------
  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TStepPointMCAnaModule::Debug() {

//-----------------------------------------------------------------------------
// bit 4: events with NHitsTF > 1
//-----------------------------------------------------------------------------
  // if (GetDebugBit(4) == 1) {
  //   GetHeaderBlock()->Print(Form("NHits(TF) = %5i",fNGenp));
  // }
}

//_____________________________________________________________________________
int TStepPointMCAnaModule::EndJob() {
  printf("----- end job: ---- %s\n",GetName());
  return 0;
}

//_____________________________________________________________________________
void TStepPointMCAnaModule::Test001() {
}

