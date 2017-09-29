///////////////////////////////////////////////////////////////////////////////
// Use Of Tmp:
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
#include "ana/TColl1DoseAnaModule.hh"

ClassImp(TColl1DoseAnaModule)
//-----------------------------------------------------------------------------
TColl1DoseAnaModule::TColl1DoseAnaModule(const char* name, const char* title):
  TStnModule(name,title)
{
  fPdgCode       = 11;
  fGeneratorCode = 28;
  fNPOT          = 1.2e20;           // "per year of running"
  fNPerPOT       = 1.;
}

//-----------------------------------------------------------------------------
TColl1DoseAnaModule::~TColl1DoseAnaModule() {
}

//-----------------------------------------------------------------------------
// histograms with the same name but for different volumes could have different number of bins
//-----------------------------------------------------------------------------
void TColl1DoseAnaModule::BookPbrAbsSpmcHistograms(HistBase_t* Hist, const char* Folder) {
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

  HBook2F(hist->fYVsX           ,"y_vs_x"     ,Form("%s: Y vs X"       ,Folder),  100,   3650, 4150, 100,  -250, 250, Folder);
  HBook2F(hist->fYVsXWtE        ,"y_vs_x_wte" ,Form("%s: Y vs X wt(E)" ,Folder),  100,   3650, 4150, 100,  -250, 250, Folder);
  HBook2F(hist->fYVsXDose       ,"y_vs_x_dose",Form("%s: Y vs X dose " ,Folder),  100,   3650, 4150, 100,  -250, 250, Folder);

}


//-----------------------------------------------------------------------------
void TColl1DoseAnaModule::BookEventHistograms(HistBase_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];

  //  char name [200];


  //  char title[200];
  EventHist_t* hist = (EventHist_t*) Hist;

  HBook1F(hist->fEventNumber,"evtnum",Form("%s: Event Number",Folder), 1000, 0,  1.e4,Folder);
  HBook1F(hist->fRunNumber  ,"runnum",Form("%s: Run   Number",Folder), 1000, 0,  1.e6,Folder);
  HBook1F(hist->fNPbrAbsSteps,"npbabssteps",Form("%s: N PbrAbs Steps"      ,Folder), 1000, 0,  1.e3,Folder);
}

//_____________________________________________________________________________
void TColl1DoseAnaModule::BookHistograms() {

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
      fHist.fPbarAbsSpmc[i] = new SpmcHist_t;
      BookPbrAbsSpmcHistograms(fHist.fPbarAbsSpmc[i],Form("Hist/%s",folder_name));
    }
  }
}

//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TColl1DoseAnaModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
  RegisterDataBlock("PbarAbsDiskSpmcBlock","TStepPointMCBlock" ,&fPbarAbsSpmcBlock);
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();

  return 0;
}


//-----------------------------------------------------------------------------
void TColl1DoseAnaModule::FillSpmcHistograms(HistBase_t* Hist, TStepPointMC* Step, SpmcData_t* Sd) {

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
  //  float z = Step->Pos()->Z()+unit.Z()*step;

  hist->fYVsX->Fill(x,y);		      // useful for stage 2
  hist->fYVsXWtE->Fill(x,y,Step->EDepTot());  // useful for stage 1

  // hist->fYVsZ->Fill(z,y);		      // useful for stage 1
  // hist->fYVsZWtE->Fill(z,y,Step->EDepTot());  // useful for stage 1
}

//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
void TColl1DoseAnaModule::FillEventHistograms(HistBase_t* Hist) {
//   double            cos_th, xv, yv, rv, zv, p;
//   TLorentzVector    mom;

  EventHist_t* hist = (EventHist_t*) Hist;

  int event_number = GetHeaderBlock()->EventNumber();
  int run_number   = GetHeaderBlock()->RunNumber();

  hist->fEventNumber->Fill(event_number);
  hist->fRunNumber->Fill(run_number);

  hist->fNPbrAbsSteps->Fill(fNPbrAbsSteps);

}

//_____________________________________________________________________________
void TColl1DoseAnaModule::FillHistograms() {

//-----------------------------------------------------------------------------
// event histograms
// EVT_0: all events
//-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0]);
//-----------------------------------------------------------------------------
// StepPointMC histograms
//-----------------------------------------------------------------------------
  for (int i=0; i<fNPbrAbsSteps; i++) {
    if (i >= kMaxNSteps) break;
    TStepPointMC* s  = fPbarAbsSpmcBlock->StepPointMC(i);
    SpmcData_t* sd  = fPbarAbsSpmcData+i;
    FillSpmcHistograms(fHist.fPbarAbsSpmc[0],s,sd);
  }
}



//_____________________________________________________________________________
int TColl1DoseAnaModule::BeginRun() {
  int rn = GetHeaderBlock()->RunNumber();
  TStntuple::Init(rn);
  return 0;
}


//_____________________________________________________________________________
int TColl1DoseAnaModule::Event(int ientry) {

  fPbarAbsSpmcBlock->GetEntry(ientry);
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  fNPbrAbsSteps = fPbarAbsSpmcBlock->NStepPoints();
  for (int i=0; i<fNPbrAbsSteps; i++) {
    if (i >= kMaxNSteps) {
      printf(" >>> PbrAbs Error: too many steps: %5i\n",fNPbrAbsSteps);
      break;
    }
    SpmcData_t* spmc_data = fPbarAbsSpmcData+i;
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
void TColl1DoseAnaModule::Debug() {

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
      TStepPointMC* s       = fPbarAbsSpmcBlock->StepPointMC(i);

      if (s->Pos()->X() > 0.10) {
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
// time to make the dose histograms
//-----------------------------------------------------------------------------
int TColl1DoseAnaModule::EndJob() {
  printf("----- end job: ---- %s\n",GetName());
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  double mev_per_joule = 1.6e-19*1.e6;
  //  double krad_per_gray = 10.;
  // double rhoCu         = 8.815e-3 ;	// kg/cm^3
  // double rhoCapton     = 1.42e-3  ;	// kg/cm^3
  double rhoAl         = 2.7e-3    ;	// kg/cm^3
//-----------------------------------------------------------------------------
// all Y:Z histograms have the same number of bins
//-----------------------------------------------------------------------------
  int nx   =  fHist.fPbarAbsSpmc[0]->fYVsXWtE->GetNbinsX();
  int ny   =  fHist.fPbarAbsSpmc[0]->fYVsXWtE->GetNbinsY();

  float dx = fHist.fPbarAbsSpmc[0]->fYVsXWtE->GetXaxis()->GetBinWidth(1)/10.; // in cm
  float dy = fHist.fPbarAbsSpmc[0]->fYVsXWtE->GetYaxis()->GetBinWidth(1)/10.; // in cm

  int n_pot_dataset = 1000000; //  GetAna()->NProcessedEvents();
//-----------------------------------------------------------------------------
// first Pbar absorber - aluminum 
// unlike in thick collimator, assume that energy is distributed uniformly here
//-----------------------------------------------------------------------------
  double mass   = rhoAl*0.025*dx*dy;       // in kg (per cm^2), thickness 250 um
  double sf     = (fNPOT*fNPerPOT)/(n_pot_dataset+1.e-12)*mev_per_joule/mass; // /krad_per_gray;

  TH2F* h1 = fHist.fPbarAbsSpmc[0]->fYVsXWtE;
  TH2F* h2 = fHist.fPbarAbsSpmc[0]->fYVsXDose;

  for (int ix=1; ix<=nx; ix++) {
    for (int iy=1; iy<=ny; iy++) {
      double x    = h1->GetBinContent(ix,iy);
      double e    = h1->GetBinError  (ix,iy);
      
      h2->SetBinContent(ix,iy,x*sf);
      h2->SetBinError  (ix,iy,e*sf);
    }
  }

  return 0;
}

//_____________________________________________________________________________
void TColl1DoseAnaModule::Test001() {
}

