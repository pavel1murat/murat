///////////////////////////////////////////////////////////////////////////////
// 2017-05-22: to run on output of Stage1 simulation
// proton is the first TSimParticle number 0
//
// fSpmcBlock : all particles reaching some final detector  
//          ("StepPointMC:: mubeam" ot "StepPointMC::DSVacuum" or smth else for Stage3)
// fVdetBlock        : all virtual detectors... ("StepPointMC::virtualdetectors" collection")
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

#include "ana/TBeamFlashAnaModule.hh"

ClassImp(TBeamFlashAnaModule)
//-----------------------------------------------------------------------------
TBeamFlashAnaModule::TBeamFlashAnaModule(const char* name, const char* title):
  TStnModule(name,title)
{
  // fPdgCode       = 11;
  // fGeneratorCode = 28;

  fSpmcBlockName = "SpmcBlock";
}

//-----------------------------------------------------------------------------
TBeamFlashAnaModule::~TBeamFlashAnaModule() {
}


//-----------------------------------------------------------------------------
void TBeamFlashAnaModule::BookSpmcHistograms(HistBase_t* Hist, const char* Folder) {
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

  HBook2F(hist->fYVsX           ,"y_vs_x"     ,Form("%s: Y vs X"       ,Folder), 100, -250,  250, 100, -250, 250, Folder);
  HBook2F(hist->fYVsZ           ,"y_vs_z"     ,Form("%s: Y vs Z"       ,Folder), 500, -250,  250, 500, -250, 250, Folder);

  HBook1F(hist->fGpPDGCode[0]   ,"gp_pdg_0"        ,Form("%s: Grand Parent PDG code",Folder), 700, -350, 350,Folder);
  HBook1F(hist->fGpPDGCode[1]   ,"gp_pdg_1"        ,Form("%s: Grand Parent PDG code",Folder), 700, -3500,3500,Folder);
  HBook2F(hist->fGpCosThVsMom   ,"gp_mom_vs_costh" ,Form("%s: Grand Parent Momentum",Folder), 500,  0, 10000, 100,-1,1,Folder);
  HBook1F(hist->fGpTheta        ,"gp_theta"        ,Form("%s: GP Theta"             ,Folder), 16,  -0.05,   3.15  ,Folder);

  HBook1F(hist->fGpPionTheta    ,"gpp_theta"        ,Form("%s: GP Pion Theta"     ,Folder), 16,  -0.05,   3.15  ,Folder);
  HBook1F(hist->fGpPionCosTh    ,"gpp_costh"        ,Form("%s: GP Pion Cos Theta" ,Folder),200,     -1,      1  ,Folder);
  HBook1F(hist->fGpPionMom      ,"gpp_mom"          ,Form("%s: GP Pion Momentum"  ,Folder),250,      0,   2500  ,Folder);
}


//-----------------------------------------------------------------------------
void TBeamFlashAnaModule::BookVDetHistograms(HistBase_t* Hist, const char* Folder) {
  //     char name [200];
  //     char title[200];

  VDetHist_t* hist = (VDetHist_t*) Hist;

  HBook1F(hist->fVolumeID,"vol_id"   ,Form("%s: VolumeID"       ,Folder), 200,    0, 200,Folder);

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

  HBook2F(hist->fYVsX           ,"y_vs_x"     ,Form("%s: Y vs X"       ,Folder), 100, -250,  250, 100, -250, 250, Folder);
  HBook2F(hist->fYVsZ           ,"y_vs_z"     ,Form("%s: Y vs Z"       ,Folder), 500, -250,  250, 500, -250, 250, Folder);

  HBook1F(hist->fGpPDGCode[0]   ,"gp_pdg_0"        ,Form("%s: Grand Parent PDG code",Folder), 700, -350, 350,Folder);
  HBook1F(hist->fGpPDGCode[1]   ,"gp_pdg_1"        ,Form("%s: Grand Parent PDG code",Folder), 700, -3500,3500,Folder);
  HBook2F(hist->fGpCosThVsMom   ,"gp_mom_vs_costh" ,Form("%s: Grand Parent Momentum",Folder), 500,  0, 10000, 100,-1,1,Folder);
  HBook1F(hist->fGpTheta        ,"gp_theta"        ,Form("%s: GP Theta"          ,Folder), 16,  -0.05,   3.15  ,Folder);

  HBook1F(hist->fGpPionTheta    ,"gpp_theta"        ,Form("%s: GP Pion Theta"     ,Folder), 16,  -0.05,   3.15  ,Folder);
  HBook1F(hist->fGpPionCosTh    ,"gpp_costh"        ,Form("%s: GP Pion Cos Theta" ,Folder),200,     -1,      1  ,Folder);
  HBook1F(hist->fGpPionMom      ,"gpp_mom"          ,Form("%s: GP Pion Momentum"  ,Folder),250,      0,   2500  ,Folder);
}


//-----------------------------------------------------------------------------
void TBeamFlashAnaModule::BookEventHistograms(HistBase_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];
  EventHist_t* hist = (EventHist_t*) Hist;

  HBook1F(hist->fEventNumber,"evtnum",Form("%s: Event Number",Folder), 1000, 0,  1.e4,Folder);
  HBook1F(hist->fRunNumber  ,"runnum",Form("%s: Run   Number",Folder), 1000, 0,  1.e6,Folder);
}

//_____________________________________________________________________________
void TBeamFlashAnaModule::BookHistograms() {

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
// book Spmc (on the final detector) histograms
//-----------------------------------------------------------------------------
  int book_spmc_histset[kNSpmcHistSets];
  for (int i=0; i<kNSpmcHistSets; i++) book_spmc_histset[i] = 0;

  book_spmc_histset[  0] = 1;		// all steps
  book_spmc_histset[  1] = 1;		// electrons
  book_spmc_histset[  2] = 1;		// positrons
  book_spmc_histset[  3] = 1;		// mu-
  book_spmc_histset[  4] = 1;		// mu+
  book_spmc_histset[  5] = 1;		// photons
  book_spmc_histset[  6] = 1;		// negative pions
  book_spmc_histset[  7] = 1;		// positive pions
  book_spmc_histset[  8] = 1;		// protons+antiprotons
  book_spmc_histset[  9] = 1;		// everything else

  book_spmc_histset[101] = 1;		// electrons with P > 50 MeV/c
  book_spmc_histset[102] = 1;		// electrons with T > 400 ns
  book_spmc_histset[103] = 1;		// electrons with p < 2 MeV/c
  book_spmc_histset[104] = 1;		// electrons with 2 < p < 3 MeV/c
  book_spmc_histset[105] = 1;		// electrons with p > 3 MeV/c
  book_spmc_histset[106] = 1;		// electrons with p > 20 MeV/c
  book_spmc_histset[107] = 1;		// electrons with p < 20 MeV/c

  book_spmc_histset[203] = 1;		// positrons with p < 2 MeV/c
  book_spmc_histset[204] = 1;		// positrons with 2 < p < 3 MeV/c
  book_spmc_histset[205] = 1;		// positrons with p > 3 MeV/c

  book_spmc_histset[301] = 1;		// muons with P > 50 MeV/c
  book_spmc_histset[310] = 1;		// muons with P < 50 MeV/c

  for (int i=0; i<kNSpmcHistSets; i++) {
    if (book_spmc_histset[i] != 0) {
      sprintf(folder_name,"spmc_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fSpmc[i] = new SpmcHist_t;
      BookSpmcHistograms(fHist.fSpmc[i],Form("Hist/%s",folder_name));
    }
  }

//-----------------------------------------------------------------------------
// book VDet (on all detectors) histograms
//-----------------------------------------------------------------------------
  int book_vdet_histset[kNVDetHistSets];
  for (int i=0; i<kNVDetHistSets; i++) book_vdet_histset[i] = 0;

  book_vdet_histset[ 91] = 1;		// electrons at VDET ID=9 - in front of ST
  book_vdet_histset[ 92] = 1;		// positrons
  book_vdet_histset[ 93] = 1;		// mu-
  book_vdet_histset[ 94] = 1;		// mu+

  book_vdet_histset[911] = 1;		// electrons at VDET ID=91 PBAR window
  book_vdet_histset[912] = 1;		// positrons
  book_vdet_histset[913] = 1;		// mu-
  book_vdet_histset[914] = 1;		// mu+

  book_vdet_histset[915] = 1;		// mu- with P< 2 MeV/c

  for (int i=0; i<kNVDetHistSets; i++) {
    if (book_vdet_histset[i] != 0) {
      sprintf(folder_name,"vdet_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fVDet[i] = new VDetHist_t;
      BookVDetHistograms(fHist.fVDet[i],Form("Hist/%s",folder_name));
    }
  }
}

//-----------------------------------------------------------------------------
void TBeamFlashAnaModule::FillSpmcHistograms(HistBase_t* Hist, TStepPointMC* Step, SpmcData_t* Sd) {

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

  double        gp_mom(-1), gp_cos_th(-10), gp_theta(-10.);
  int           gp_pdg_code(-1000000);
  TSimParticle* gp = Sd->fGParent;
  
  if (gp) {
    TLorentzVector* v1 = &fProton->fStartMom;
    TLorentzVector* v2 = &gp->fStartMom;

    gp_cos_th = (v1->Px()*v2->Px()+v1->Py()*v2->Py()+v1->Pz()*v2->Pz())/v1->P()/v2->P()/(1+1.e-10);
    gp_theta  = TMath::ACos(gp_cos_th);
    
    gp_pdg_code = gp->fPdgCode;
    gp_mom      = gp->fStartMom.P();
  }

  hist->fGpPDGCode[0]->Fill(gp_pdg_code);
  hist->fGpPDGCode[1]->Fill(gp_pdg_code);
  hist->fGpTheta->Fill(gp_theta);
  hist->fGpCosThVsMom->Fill(gp_mom,gp_cos_th);

  if (abs(gp_pdg_code) == 211) {
    hist->fGpPionCosTh->Fill(gp_cos_th);
    hist->fGpPionTheta->Fill(gp_theta);
    hist->fGpPionMom->Fill(gp_mom);
  }
  
}

//-----------------------------------------------------------------------------
void TBeamFlashAnaModule::FillVDetHistograms(HistBase_t* Hist, TStepPointMC* Step, VDetData_t* Vdd) {

  VDetHist_t* hist = (VDetHist_t*) Hist;

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

  hist->fYVsX->Fill(Vdd->fX,Vdd->fY);

  double gp_mom(-1), gp_cos_th(-10.), gp_theta(-10.);
  int    gp_pdg_code(-1000000);
  TSimParticle* gp = Vdd->fGParent;
  
  if (gp) {
    TLorentzVector* v1 = &fProton->fStartMom;
    TLorentzVector* v2 = &gp->fStartMom;

    gp_cos_th   = (v1->Px()*v2->Px()+v1->Py()*v2->Py()+v1->Pz()*v2->Pz())/v1->P()/v2->P()/(1+1.e-10);
    gp_theta    = TMath::ACos(gp_cos_th);
    gp_pdg_code = gp->fPdgCode;
    gp_mom      = gp->fStartMom.P();
  }

  hist->fMomentum->Fill(Vdd->fP);
  hist->fCosTh->Fill(Vdd->fCosTh);

  hist->fGpPDGCode[0]->Fill(gp_pdg_code);
  hist->fGpPDGCode[1]->Fill(gp_pdg_code);
  hist->fGpCosThVsMom->Fill(gp_mom,gp_cos_th);
  hist->fGpTheta->Fill(gp_theta);
  if (abs(gp_pdg_code) == 211) {
    hist->fGpPionCosTh->Fill(gp_cos_th);
    hist->fGpPionTheta->Fill(gp_theta);
    hist->fGpPionMom->Fill(gp_mom);
  }

}

//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
void TBeamFlashAnaModule::FillEventHistograms(HistBase_t* Hist) {
//   double            cos_th, xv, yv, rv, zv, p;
//   TLorentzVector    mom;

  EventHist_t* hist = (EventHist_t*) Hist;

  int event_number = GetHeaderBlock()->EventNumber();
  int run_number   = GetHeaderBlock()->RunNumber();

  hist->fEventNumber->Fill(event_number);
  hist->fRunNumber->Fill(run_number);

}

//_____________________________________________________________________________
void TBeamFlashAnaModule::FillHistograms() {

//-----------------------------------------------------------------------------
// event histograms
//
// EVT_0: all events
//-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0]);
//-----------------------------------------------------------------------------
// Spmc histograms - this assumes that all step points are on the same 
// virtual detector - "mubeam" or "DSVacuum" collection
//-----------------------------------------------------------------------------
  TStepPointMC* spmc;
  SpmcData_t*   sd;

  int nsteps = fSpmcBlock->NStepPoints();
  for (int i=0; i<nsteps; i++) {
    spmc = fSpmcBlock->StepPointMC(i);
    sd   = fSpmcData+i;

    float p = spmc->Mom()->Mag();
    float t = spmc->Time();

    FillSpmcHistograms(fHist.fSpmc[0],spmc,sd);

    if      (spmc->PDGCode() ==   11) {
      FillSpmcHistograms(fHist.fSpmc[1],spmc,sd);
      if (p >  50) FillSpmcHistograms(fHist.fSpmc[101],spmc,sd);
      if (t > 400) FillSpmcHistograms(fHist.fSpmc[102],spmc,sd);
      if (p <   2)              FillSpmcHistograms(fHist.fSpmc[103],spmc,sd);
      if ((p >=  2) && (p < 3)) FillSpmcHistograms(fHist.fSpmc[104],spmc,sd);
      if (p >   3)              FillSpmcHistograms(fHist.fSpmc[105],spmc,sd);
      if (p >= 20)              FillSpmcHistograms(fHist.fSpmc[106],spmc,sd);
      if (p <  20)              FillSpmcHistograms(fHist.fSpmc[107],spmc,sd);

      if (GetDebugBit(4) == 1) {
	if (sd->fGParent == NULL) {
	  GetHeaderBlock()->Print(Form("electron GP = NULL"));
	}
      }
      if (GetDebugBit(5) == 1) {
	if (sd->fGParent && (sd->fGParent->PDGCode() > 2500)) {
	  GetHeaderBlock()->Print(Form("electron GP PDG code: %i",sd->fGParent->PDGCode()));
	}
      }
    }
    else if (spmc->PDGCode() ==  -11) {
      FillSpmcHistograms(fHist.fSpmc[2],spmc,sd);
      if (p <   2)              FillSpmcHistograms(fHist.fSpmc[203],spmc,sd);
      if ((p >=  2) && (p < 3)) FillSpmcHistograms(fHist.fSpmc[204],spmc,sd);
      if (p >   3)              FillSpmcHistograms(fHist.fSpmc[205],spmc,sd);
    }
    else if (spmc->PDGCode() ==   13) {
      FillSpmcHistograms(fHist.fSpmc[3],spmc,sd);
      if (p >=  50) FillSpmcHistograms(fHist.fSpmc[301],spmc,sd);
      if (p <   50) FillSpmcHistograms(fHist.fSpmc[310],spmc,sd);

      if (GetDebugBit(3) == 1) {
	p = sd->fGParent->fStartMom.P();
	if (p > 8000) {
	  GetHeaderBlock()->Print(Form("muon GP mom : %10.3f",p));
	}
      }
      if (GetDebugBit(4) == 1) {
	if (sd->fGParent == NULL) {
	  GetHeaderBlock()->Print(Form("muon GP = NULL"));
	}
      }
      if (GetDebugBit(5) == 1) {
	if (sd->fGParent && (sd->fGParent->PDGCode() > 2500)) {
	  GetHeaderBlock()->Print(Form("muon     GP PDG code: %i",sd->fGParent->PDGCode()));
	}
      }
    }
    else if (spmc->PDGCode()     ==   -13) FillSpmcHistograms(fHist.fSpmc[4],spmc,sd);
    else if (spmc->PDGCode()     ==    22) FillSpmcHistograms(fHist.fSpmc[5],spmc,sd);
    else if (spmc->PDGCode()     ==  -211) FillSpmcHistograms(fHist.fSpmc[6],spmc,sd);
    else if (spmc->PDGCode()     ==   211) FillSpmcHistograms(fHist.fSpmc[7],spmc,sd);
    else if (abs(spmc->PDGCode() == 2212)) FillSpmcHistograms(fHist.fSpmc[8],spmc,sd);
    else                                   FillSpmcHistograms(fHist.fSpmc[9],spmc,sd);

  }
//-----------------------------------------------------------------------------
// VDet histograms - for all virtual detectors
// currently, fill only for VDET ID=9 and 91 (TS1 proton absorber window)
//-----------------------------------------------------------------------------
  TStepPointMC* spvd;
  VDetData_t*   vdd;

  int n_vd_pts = fVDetBlock->NStepPoints();
  for (int i=0; i<n_vd_pts; i++) {
    spvd = fVDetBlock->StepPointMC(i);
    vdd  = fVDetData+i;

  //   float p = spvd->Mom()->Mag();
  //   float t = spvd->Time();

    int vol_id = spvd->VolumeID();

    if      (vol_id == 9) {
       if      (spvd->PDGCode() ==   11) FillVDetHistograms(fHist.fVDet[ 91],spvd,vdd);
       if      (spvd->PDGCode() ==  -11) FillVDetHistograms(fHist.fVDet[ 92],spvd,vdd);
       if      (spvd->PDGCode() ==   13) FillVDetHistograms(fHist.fVDet[ 93],spvd,vdd);
       if      (spvd->PDGCode() ==  -13) FillVDetHistograms(fHist.fVDet[ 94],spvd,vdd);
    }
    else if (vol_id == 91) {
       if      (spvd->PDGCode() ==   11) FillVDetHistograms(fHist.fVDet[911],spvd,vdd);
       if      (spvd->PDGCode() ==  -11) FillVDetHistograms(fHist.fVDet[912],spvd,vdd);
       if      (spvd->PDGCode() ==   13) FillVDetHistograms(fHist.fVDet[913],spvd,vdd);
       if      (spvd->PDGCode() ==  -13) FillVDetHistograms(fHist.fVDet[914],spvd,vdd);

       if (vdd->fP < 2) {
	 if (spvd->PDGCode() ==   11) FillVDetHistograms(fHist.fVDet[915],spvd,vdd);
       }
    }
  }
}



//_____________________________________________________________________________
int TBeamFlashAnaModule::BeginRun() {
  int rn = GetHeaderBlock()->RunNumber();
  TStntuple::Init(rn);
  return 0;
}


//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TBeamFlashAnaModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
  RegisterDataBlock(fSpmcBlockName.Data(),"TStepPointMCBlock",&fSpmcBlock);
  RegisterDataBlock("VdetBlock"          ,"TStepPointMCBlock",&fVDetBlock);
  //  RegisterDataBlock("StepPointMCBlock" ,"TStepPointMCBlock" ,&fSpmcBlock);
  RegisterDataBlock("SimpBlock"          ,"TSimpBlock"       ,&fSimpBlock);
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();

  return 0;
}

//_____________________________________________________________________________
int TBeamFlashAnaModule::Event(int ientry) {

  //  double                p;
  //  TLorentzVector        mom;

  fSpmcBlock->GetEntry(ientry);
  fVDetBlock->GetEntry(ientry);
  fSimpBlock->GetEntry(ientry);

  fProton = fSimpBlock->Particle(0);
//-----------------------------------------------------------------------------
// loop over "steps" in SpmcBlock - assume that this block has 
// one steppoint per particle, like "DsVolume"
//-----------------------------------------------------------------------------
  int nsteps = fSpmcBlock->NStepPoints();
  //  int nsimp  = fSimpBlock->NParticles();

  for (int i=0; i<nsteps; i++) {
    TStepPointMC* s       = fSpmcBlock->StepPointMC(i);
    SpmcData_t* spmc_data = fSpmcData+i;
    int id                = s->GetUniqueID();
//-----------------------------------------------------------------------------
// find particle in SimpBlock
//-----------------------------------------------------------------------------
    spmc_data->fParticle  = fSimpBlock->FindParticle(id);
    spmc_data->fParent    = NULL;
    spmc_data->fGParent   = NULL;
					// search for particle's oldest parent with parent code = 1
    int parent_id = s->ParentSimID();
    int pdg_code  = s->PDGCode();

    TSimParticle* parent = fSimpBlock->FindParticle(parent_id);

    while(parent) {
      parent_id = parent->ParentID();
      TSimParticle* p = fSimpBlock->FindParticle(parent_id);
      if (p            == NULL) break; 
      if ((spmc_data->fParent == NULL) && (p->PDGCode() != pdg_code)) {
	spmc_data->fParent = p;
      }

      if (p->PDGCode() == 2212) break;
      parent = p;
//-----------------------------------------------------------------------------
// sometimes history tree includes a scattered proton, which produces a pion,
// decaying into a muon. In this case call pion a 'grandparent'
//-----------------------------------------------------------------------------
      if (parent_id == 1) break;
    }
    spmc_data->fGParent = parent;
  }

  int nvds = fVDetBlock->NStepPoints();

  for (int i=0; i<nvds; i++) {
    TStepPointMC* s       = fVDetBlock->StepPointMC(i);
    int id                = s->GetUniqueID();
    VDetData_t* vdet_data = fVDetData+i;
//-----------------------------------------------------------------------------
// find particle in SimpBlock
//-----------------------------------------------------------------------------
    vdet_data->fParticle  = fSimpBlock->FindParticle(id);
    vdet_data->fParent    = NULL;
    vdet_data->fGParent   = NULL;
    vdet_data->fP         = s->Mom()->Mag();
					// search for particle's oldest parent with parent code = 1
    int parent_id = s->ParentSimID();
    int pdg_code  = s->PDGCode();

    TSimParticle* parent = fSimpBlock->FindParticle(parent_id);

    while(parent) {
      parent_id = parent->ParentID();
      TSimParticle* p = fSimpBlock->FindParticle(parent_id);
      if (p            == NULL) break; 
      if ((vdet_data->fParent == NULL) && (p->PDGCode() != pdg_code)) {
	vdet_data->fParent = p;
      }

      if (p->PDGCode() == 2212) break;
      parent = p;
//-----------------------------------------------------------------------------
// sometimes history tree includes a scattered proton, which produces a pion,
// decaying into a muon. In this case call pion a 'grandparent'
//-----------------------------------------------------------------------------
      if (parent_id == 1) break;
    }
    vdet_data->fGParent = parent;

    int vol_id = s->VolumeID();
    if (vol_id == 9) {
//-----------------------------------------------------------------------------
// TS1 proton absorber window
//-----------------------------------------------------------------------------
      vdet_data->fX     = s->Pos()->X()-3904;
      vdet_data->fY     = s->Pos()->Y() ;
      vdet_data->fCosTh = s->Mom()->Pt()/vdet_data->fP;
    }
    else if (vol_id == 91) {
//-----------------------------------------------------------------------------
// TS1 proton absorber window
//-----------------------------------------------------------------------------
      vdet_data->fX     = s->Pos()->X()-3904;
      vdet_data->fY     = s->Pos()->Y() ;
      vdet_data->fCosTh = s->Mom()->Pt()/vdet_data->fP;
    }
  }

  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TBeamFlashAnaModule::Debug() {

//-----------------------------------------------------------------------------
// bit 4: events with NHitsTF > 1
//-----------------------------------------------------------------------------
  // if (GetDebugBit(4) == 1) {
  //   GetHeaderBlock()->Print(Form("NHits(TF) = %5i",fNGenp));
  // }
}

//_____________________________________________________________________________
int TBeamFlashAnaModule::EndJob() {
  printf("----- end job: ---- %s\n",GetName());
  return 0;
}

//_____________________________________________________________________________
void TBeamFlashAnaModule::Test001() {
}

