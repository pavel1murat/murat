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
// 4  : hits in Pbar absorber with X>0.1 mm 
// 5  : hits in Coll31        with X>825 mm 
// 6  : UNUSED
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
#include "ana/TColl3DoseAnaModule.hh"

ClassImp(TColl3DoseAnaModule)
//-----------------------------------------------------------------------------
TColl3DoseAnaModule::TColl3DoseAnaModule(const char* name, const char* title):
  TStnModule(name,title)
{
  fPdgCode       = 11;
  fGeneratorCode = 28;
  fNPOT          = 1.2e20;           // "per year of running"
  fNPerPOT       = 1.;
}

//-----------------------------------------------------------------------------
TColl3DoseAnaModule::~TColl3DoseAnaModule() {
}


//-----------------------------------------------------------------------------
void TColl3DoseAnaModule::BookColl31SpmcHistograms(HistBase_t* Hist, const char* Folder) {
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

  HBook2F(hist->fYVsX           ,"y_vs_x"     ,Form("%s: Y vs X"       ,Folder), 1000,     0, 1000,  60, -300, 300, Folder);
  HBook2F(hist->fYVsXWtE        ,"y_vs_x_wte" ,Form("%s: Y vs X wt(E)" ,Folder), 1000,     0, 1000,  60, -300, 300, Folder);

  HBook2F(hist->fYVsZ           ,"y_vs_z"     ,Form("%s: Y vs Z"       ,Folder),  60,  -300,  300,  60, -300, 300, Folder);
  HBook2F(hist->fYVsZWtE        ,"y_vs_z_wte" ,Form("%s: Y vs Z wt(E)" ,Folder),  60,  -300,  300,  60, -300, 300, Folder);

  for (int i=0; i<kNSlices; i++) {
    HBook2F(hist->fYVsZWtESlice[i],Form("y_vs_z_wte_%i",i) ,Form("%s: Y vs Z wt(E) slice %i" ,Folder,i),  
	    60,  -300,  300,  60, -300, 300, Folder);

    HBook2F(hist->fYVsZDose[i],Form("y_vs_z_dose_%i",i) ,Form("%s: Y vs Z Doze slice %i" ,Folder,i),  
	    60,  -300,  300,  60, -300, 300, Folder);
  }
}

//-----------------------------------------------------------------------------
void TColl3DoseAnaModule::BookColl32SpmcHistograms(HistBase_t* Hist, const char* Folder) {
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

  HBook2F(hist->fYVsX           ,"y_vs_x"     ,Form("%s: Y vs X"       ,Folder), 1000,  -1000,   0,  60, -300, 300, Folder);
  HBook2F(hist->fYVsXWtE        ,"y_vs_x_wte" ,Form("%s: Y vs X wt(E)" ,Folder), 1000,  -1000,   0,  60, -300, 300, Folder);

  HBook2F(hist->fYVsZ           ,"y_vs_z"     ,Form("%s: Y vs Z"       ,Folder),  60,   -300, 300,  60, -300, 300, Folder);
  HBook2F(hist->fYVsZWtE        ,"y_vs_z_wte" ,Form("%s: Y vs Z wt(E)" ,Folder),  60,   -300, 300,  60, -300, 300, Folder);

  for (int i=0; i<kNSlices; i++) {
    HBook2F(hist->fYVsZWtESlice[i],Form("y_vs_z_wte_%i",i) ,Form("%s: Y vs Z wt(E) slice %i" ,Folder,i),  
	    60,  -300,  300,  60, -300, 300, Folder);

    HBook2F(hist->fYVsZDose[i],Form("y_vs_z_dose_%i",i) ,Form("%s: Y vs Z Doze slice %i" ,Folder,i),  
	    60,  -300,  300,  60, -300, 300, Folder);
  }
}

//-----------------------------------------------------------------------------
// histograms with the same name but for different volumes could have different number of bins
//-----------------------------------------------------------------------------
void TColl3DoseAnaModule::BookPbrAbsSpmcHistograms(HistBase_t* Hist, const char* Folder) {
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

  HBook2F(hist->fYVsX           ,"y_vs_x"     ,Form("%s: Y vs X"       ,Folder),  100,   -0.5,  0.5, 300,  -300, 300, Folder);
  HBook2F(hist->fYVsXWtE        ,"y_vs_x_wte" ,Form("%s: Y vs X wt(E)" ,Folder),  100,   -0.5,  0.5, 300,  -300, 300, Folder);

  HBook2F(hist->fYVsZ           ,"y_vs_z"     ,Form("%s: Y vs Z"       ,Folder),  60,  -300,   300,  60,  -300, 300, Folder);
  HBook2F(hist->fYVsZWtE        ,"y_vs_z_wte" ,Form("%s: Y vs Z wt(E)" ,Folder),  60,  -300,   300,  60,  -300, 300, Folder);

  for (int i=0; i<kNSlices; i++) {
    HBook2F(hist->fYVsZWtESlice[i],Form("y_vs_z_wte_%i",i) ,Form("%s: Y vs Z wt(E) slice %i" ,Folder,i),  
	    60,  -300,  300,  60, -300, 300, Folder);

    HBook2F(hist->fYVsZDose[i],Form("y_vs_z_dose_%i",i) ,Form("%s: Y vs Z Doze slice %i" ,Folder,i),  
	    60,  -300,  300,  60, -300, 300, Folder);
  }
}


//-----------------------------------------------------------------------------
void TColl3DoseAnaModule::BookEventHistograms(HistBase_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];

  //  char name [200];
  //  char title[200];
  EventHist_t* hist = (EventHist_t*) Hist;

  HBook1F(hist->fEventNumber,"evtnum",Form("%s: Event Number",Folder), 1000, 0,  1.e4,Folder);
  HBook1F(hist->fRunNumber  ,"runnum",Form("%s: Run   Number",Folder), 1000, 0,  1.e6,Folder);
  HBook1F(hist->fNColl31Steps,"ncoll31steps",Form("%s: N Coll31 Steps"      ,Folder), 1000, 0,  1.e3,Folder);
  HBook1F(hist->fNColl32Steps,"ncoll32steps",Form("%s: N Coll32 Steps"      ,Folder), 1000, 0,  1.e3,Folder);
  HBook1F(hist->fNPbrAbsSteps,"npbrabssteps",Form("%s: N PbrAbs Steps"      ,Folder), 1000, 0,  1.e3,Folder);
}

//_____________________________________________________________________________
void TColl3DoseAnaModule::BookHistograms() {

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
// book Coll31 StepPointMC histograms
//-----------------------------------------------------------------------------
  int book_coll31_spmc_histset[kNSpmcHistSets];
  for (int i=0; i<kNSpmcHistSets; i++) book_coll31_spmc_histset[i] = 0;

  book_coll31_spmc_histset[0] = 1;		// all clusters

  for (int i=0; i<kNSpmcHistSets; i++) {
    if (book_coll31_spmc_histset[i] != 0) {
      sprintf(folder_name,"coll31_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fColl31Spmc[i] = new SpmcHist_t;
      BookColl31SpmcHistograms(fHist.fColl31Spmc[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book Coll32 StepPointMC histograms
//-----------------------------------------------------------------------------
  int book_coll32_spmc_histset[kNSpmcHistSets];
  for (int i=0; i<kNSpmcHistSets; i++) book_coll32_spmc_histset[i] = 0;

  book_coll32_spmc_histset[0] = 1;		// all clusters

  for (int i=0; i<kNSpmcHistSets; i++) {
    if (book_coll32_spmc_histset[i] != 0) {
      sprintf(folder_name,"coll32_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fColl32Spmc[i] = new SpmcHist_t;
      BookColl32SpmcHistograms(fHist.fColl32Spmc[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book PbarAbs StepPointMC histograms
//-----------------------------------------------------------------------------
  int book_pbrabs_spmc_histset[kNSpmcHistSets];
  for (int i=0; i<kNSpmcHistSets; i++) book_pbrabs_spmc_histset[i] = 0;

  book_pbrabs_spmc_histset[0] = 1;		// all clusters

  for (int i=0; i<kNSpmcHistSets; i++) {
    if (book_pbrabs_spmc_histset[i] != 0) {
      sprintf(folder_name,"pbrabs_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fPbrAbsSpmc[i] = new SpmcHist_t;
      BookPbrAbsSpmcHistograms(fHist.fPbrAbsSpmc[i],Form("Hist/%s",folder_name));
    }
  }
}

//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TColl3DoseAnaModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
  RegisterDataBlock("Coll31SpmcBlock" ,"TStepPointMCBlock" ,&fColl31SpmcBlock);
  RegisterDataBlock("Coll32SpmcBlock" ,"TStepPointMCBlock" ,&fColl32SpmcBlock);
  RegisterDataBlock("PbarAbsSpmcBlock","TStepPointMCBlock" ,&fPbrAbsSpmcBlock);
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();

  return 0;
}


//-----------------------------------------------------------------------------
void TColl3DoseAnaModule::FillSpmcHistograms(HistBase_t* Hist, TStepPointMC* Step, SpmcData_t* Sd) {

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

  double step = Step->StepLength();

  double p = Step->Mom()->Mag();
  hist->fMomentum->Fill(p);

  TVector3 unit = Step->Mom()->Unit();
					// assume energy is deposited in the end of step
  float x = Step->Pos()->X()+unit.X()*step;
  float y = Step->Pos()->Y()+unit.Y()*step;
  float z = Step->Pos()->Z()+unit.Z()*step;

  hist->fYVsX->Fill(x,y);		      // useful for stage 2
  hist->fYVsXWtE->Fill(x,y,Step->EDepTot());  // useful for stage 1

  hist->fYVsZ->Fill(z,y);		      // useful for stage 1
  hist->fYVsZWtE->Fill(z,y,Step->EDepTot());  // useful for stage 1

  if ((Sd->fBinX >= 0) && (Sd->fBinX <= kNSlices)) {
    hist->fYVsZWtESlice[Sd->fBinX]->Fill(z,y,Step->EDepTot());
  }

}

//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
void TColl3DoseAnaModule::FillEventHistograms(HistBase_t* Hist) {
//   double            cos_th, xv, yv, rv, zv, p;
//   TLorentzVector    mom;

  EventHist_t* hist = (EventHist_t*) Hist;

  int event_number = GetHeaderBlock()->EventNumber();
  int run_number   = GetHeaderBlock()->RunNumber();

  hist->fEventNumber->Fill(event_number);
  hist->fRunNumber->Fill(run_number);

  hist->fNColl31Steps->Fill(fNColl31Steps);
  hist->fNColl32Steps->Fill(fNColl32Steps);
  hist->fNPbrAbsSteps->Fill(fNPbrAbsSteps);

}

//_____________________________________________________________________________
void TColl3DoseAnaModule::FillHistograms() {

//-----------------------------------------------------------------------------
// event histograms
// EVT_0: all events
//-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0]);
//-----------------------------------------------------------------------------
// StepPointMC histograms
//-----------------------------------------------------------------------------
  for (int i=0; i<fNColl31Steps; i++) {
    if (i >= kMaxNSteps) break;

    TStepPointMC* s  = fColl31SpmcBlock->StepPointMC(i);

    SpmcData_t*   sd = fColl31SpmcData+i;

    sd->fBinX = (825.-s->Pos()->X());  // should be greater than zero

    if (sd->fBinX < 0) {
      printf(" ERROR : X = %12.5f\n",s->Pos()->X());
    }

    FillSpmcHistograms(fHist.fColl31Spmc[0],s,sd);
  }

  for (int i=0; i<fNColl32Steps; i++) {
    if (i >= kMaxNSteps) break;
    TStepPointMC* s  = fColl32SpmcBlock->StepPointMC(i);
    SpmcData_t*   sd = fColl32SpmcData+i;

    sd->fBinX = (-25.-s->Pos()->X());  // should be greater than zero

    FillSpmcHistograms(fHist.fColl32Spmc[0],s,sd);
  }

  for (int i=0; i<fNPbrAbsSteps; i++) {
    if (i >= kMaxNSteps) break;
    TStepPointMC* s  = fPbrAbsSpmcBlock->StepPointMC(i);
    SpmcData_t*   sd = fPbrAbsSpmcData+i;

    sd->fBinX = (0.11-s->Pos()->X())/0.01;  // all in mm, should be greater than zero - bin=10um

    FillSpmcHistograms(fHist.fPbrAbsSpmc[0],s,sd);
  }
}



//_____________________________________________________________________________
int TColl3DoseAnaModule::BeginRun() {
  int rn = GetHeaderBlock()->RunNumber();
  TStntuple::Init(rn);
  return 0;
}


//_____________________________________________________________________________
int TColl3DoseAnaModule::Event(int ientry) {

  //  double                p;
  //  TLorentzVector        mom;

  fColl31SpmcBlock->GetEntry(ientry);
  fColl32SpmcBlock->GetEntry(ientry);
  fPbrAbsSpmcBlock->GetEntry(ientry);
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
//  TStepPointMC* spmc;

  fNColl31Steps = fColl31SpmcBlock->NStepPoints();

  for (int i=0; i<fNColl31Steps; i++) {
    if (i >= kMaxNSteps) {
      printf(" >>> COLL31 Error: too many steps: %5i\n",fNColl31Steps);
      break;
    }
    //    TStepPointMC* s       = fSpmcBlock->StepPointMC(i);
    SpmcData_t* spmc_data = fColl31SpmcData+i;
    //    int id                = s->GetUniqueID();
//-----------------------------------------------------------------------------
// find particle in SimpBlock
//-----------------------------------------------------------------------------
    spmc_data->fParticle  = NULL; // fSimpBlock->FindParticle(id);
    spmc_data->fParent    = NULL;
    spmc_data->fGParent   = NULL;
  }

  fNColl32Steps = fColl32SpmcBlock->NStepPoints();
  for (int i=0; i<fNColl32Steps; i++) {
    if (i >= kMaxNSteps) {
      printf(" >>> COLL32 Error: too many steps: %5i\n",fNColl32Steps);
      break;
    }
    //    TStepPointMC* s       = fSpmcBlock->StepPointMC(i);
    SpmcData_t* spmc_data = fColl32SpmcData+i;
    //    int id                = s->GetUniqueID();
//-----------------------------------------------------------------------------
// find particle in SimpBlock
//-----------------------------------------------------------------------------
    spmc_data->fParticle  = NULL; // fSimpBlock->FindParticle(id);
    spmc_data->fParent    = NULL;
    spmc_data->fGParent   = NULL;
  }

  fNPbrAbsSteps = fPbrAbsSpmcBlock->NStepPoints();
  for (int i=0; i<fNPbrAbsSteps; i++) {
    if (i >= kMaxNSteps) {
      printf(" >>> PbrAbs Error: too many steps: %5i\n",fNPbrAbsSteps);
      break;
    }
    //    TStepPointMC* s       = fSpmcBlock->StepPointMC(i);
    SpmcData_t* spmc_data = fPbrAbsSpmcData+i;
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
void TColl3DoseAnaModule::Debug() {

//-----------------------------------------------------------------------------
// bit 4: 
//-----------------------------------------------------------------------------
  if (GetDebugBit(4) == 1) {

    int banner_printed(0);

    for (int i=0; i<fNPbrAbsSteps; i++) {
      if (i >= kMaxNSteps) {
	printf(" >>> PbrAbs Error: too many steps: %5i\n",fNPbrAbsSteps);
	break;
      }
      TStepPointMC* s       = fPbrAbsSpmcBlock->StepPointMC(i);

      if (s->Pos()->X() > 0.10) {
	if (banner_printed == 0) {
	  s->Print("banner");
	  banner_printed = 1;
	}
	s->Print("data");
      }
    }
  }

//-----------------------------------------------------------------------------
// bit 5: 
//-----------------------------------------------------------------------------
  if (GetDebugBit(5) == 1) {

    int banner_printed(0);

    for (int i=0; i<fNColl31Steps; i++) {
      if (i >= kMaxNSteps) {
	printf(" >>> PbrAbs Error: too many steps: %5i\n",fNColl31Steps);
	break;
      }
      TStepPointMC* s       = fColl31SpmcBlock->StepPointMC(i);

      if (s->Pos()->X() > 820) {
	if (banner_printed == 0) {
	  s->Print("banner");
	  banner_printed = 1;
	}
	s->Print("data");
      }
    }
  }

}

//-----------------------------------------------------------------------------
// this is the time to make the dose histograms
//-----------------------------------------------------------------------------
int TColl3DoseAnaModule::EndJob() {
  printf("----- end job: ---- %s\n",GetName());
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  double mev_per_joule = 1.6e-19*1.e6;
  double krad_per_gray = 10.;
  double rhoCu         = 8.815e-3 ;	// kg/cm^3
  double rhoCapton     = 1.42e-3  ;	// kg/cm^3
//-----------------------------------------------------------------------------
// all Y:Z histograms have the same number of bins
//-----------------------------------------------------------------------------
  int nx   =  fHist.fColl31Spmc[0]->fYVsZWtESlice[0]->GetNbinsX();
  int ny   =  fHist.fColl31Spmc[0]->fYVsZWtESlice[0]->GetNbinsY();

  int n_pot_dataset = 200000; //  GetAna()->NProcessedEvents();

  for (int i=0; i<kNSlices; i++) {

    float mass = rhoCu*0.1;                          // in g/cm^2
    float sf   = (fNPOT*fNPerPOT)/(n_pot_dataset+1.e-12)*mev_per_joule/mass/krad_per_gray;

//-----------------------------------------------------------------------------
// Coll31
//-----------------------------------------------------------------------------
    TH2F* h1 = fHist.fColl31Spmc[0]->fYVsZWtESlice[i];
    TH2F* h2 = fHist.fColl31Spmc[0]->fYVsZDose    [i];

    for (int ix=1; ix<=nx; ix++) {
      for (int iy=1; iy<=ny; iy++) {
	double x    = h1->GetBinContent(ix,iy);
	double e    = h1->GetBinError  (ix,iy);

	h2->SetBinContent(ix,iy,x*sf);
	h2->SetBinContent(ix,iy,e*sf);
      }
    }
//-----------------------------------------------------------------------------
// Coll32 : the same scale factor
//-----------------------------------------------------------------------------
    h1 = fHist.fColl32Spmc[0]->fYVsZWtESlice[i];
    h2 = fHist.fColl32Spmc[0]->fYVsZDose    [i];

    for (int ix=1; ix<=nx; ix++) {
      for (int iy=1; iy<=ny; iy++) {
	double x    = h1->GetBinContent(ix,iy);
	double e    = h1->GetBinError  (ix,iy);

	h2->SetBinContent(ix,iy,x*sf);
	h2->SetBinContent(ix,iy,e*sf);
      }
    }
//-----------------------------------------------------------------------------
// Pbar absorber - capton, different Dx bin - 10 um
// unlike in thick collimator, assume that energy is distributed uniformly here
//-----------------------------------------------------------------------------
    mass   = rhoCapton*0.02;       // in g/cm^2, slice thickness 200 um
    sf     = (fNPOT*fNPerPOT)/(n_pot_dataset+1.e-12)*mev_per_joule/mass/krad_per_gray;

    h1 = fHist.fPbrAbsSpmc[0]->fYVsZWtE;
    h2 = fHist.fPbrAbsSpmc[0]->fYVsZDose    [i];

    for (int ix=1; ix<=nx; ix++) {
      for (int iy=1; iy<=ny; iy++) {
	double x    = h1->GetBinContent(ix,iy);
	double e    = h1->GetBinError  (ix,iy);

	h2->SetBinContent(ix,iy,x*sf);
	h2->SetBinContent(ix,iy,e*sf);
      }
    }
  }

  return 0;
}

//_____________________________________________________________________________
void TColl3DoseAnaModule::Test001() {
}

