///////////////////////////////////////////////////////////////////////////////
// 2017-05-22: to run on output of Stage1 simulation
// proton is the first TSimParticle number 0
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

#include "ana/TG4ValidationModule.hh"

ClassImp(TG4ValidationModule)
//-----------------------------------------------------------------------------
TG4ValidationModule::TG4ValidationModule(const char* name, const char* title):
  TStnModule(name,title)
{
  // fPdgCode       = 11;
  // fGeneratorCode = 28;
}

//-----------------------------------------------------------------------------
TG4ValidationModule::~TG4ValidationModule() {
}


//-----------------------------------------------------------------------------
// Momentum[1]: for comparison with HARP data
//-----------------------------------------------------------------------------
void TG4ValidationModule::BookSimpHistograms(HistBase_t* Hist, const char* Folder) {
  //     char name [200];
  //     char title[200];
  SimpHist_t* hist = (SimpHist_t*) Hist;

  HBook1F(hist->fPDGCode[0] ,"pdg_0" ,Form("%s: PDG code"       ,Folder),  700,  -350,   350,Folder);
  HBook1F(hist->fPDGCode[1] ,"pdg_1" ,Form("%s: PDG code"       ,Folder),  700, -3500,  3500,Folder);

  HBook1F(hist->fMomentum[0],"mom_0" ,Form("%s: Momentum"       ,Folder),  1000,   0,  10000,Folder);

  float lower_limit[] = {
    0., 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.60, 0.70, 0.80
  };

  HBook1F(hist->fMomentum[1],"mom_1" ,Form("%s: Momentum, GeV/c",Folder),   13,   lower_limit, Folder);

  HBook1F(hist->fCosTheta,   "costh" ,Form("%s: CosTh   "       ,Folder),   200,  -1,      1, Folder);
  HBook1F(hist->fTheta   ,   "theta" ,Form("%s: Theta, rad"     ,Folder),    16,  -0.05,    3.15, Folder);
}


//-----------------------------------------------------------------------------
void TG4ValidationModule::BookEventHistograms(HistBase_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];
  EventHist_t* hist = (EventHist_t*) Hist;

  HBook1F(hist->fEventNumber,"evtnum",Form("%s: Event Number",Folder), 1000, 0,  1.e4,Folder);
  HBook1F(hist->fRunNumber  ,"runnum",Form("%s: Run   Number",Folder), 1000, 0,  1.e6,Folder);
  HBook1F(hist->fNPiMinus   ,"npim"  ,Form("%s: N(pi minus)" ,Folder), 50  , 0,  50  ,Folder);
  HBook1F(hist->fNPiPlus    ,"npip"  ,Form("%s: N(pi plus)"  ,Folder), 50  , 0,  50  ,Folder);
  HBook1F(hist->fNPi0       ,"npi0"  ,Form("%s: N(pi0)"      ,Folder), 50  , 0,  50  ,Folder);
  HBook1F(hist->fNPions     ,"npions",Form("%s: N(pions)"    ,Folder), 50  , 0,  50  ,Folder);
  HBook1F(hist->fNProtons   ,"nprot" ,Form("%s: N(protons)"  ,Folder), 50  , 0,  50  ,Folder);
  HBook1F(hist->fNNeutrons  ,"nneut" ,Form("%s: N(neutrons)" ,Folder), 50  , 0,  50  ,Folder);
}

//_____________________________________________________________________________
void TG4ValidationModule::BookHistograms() {

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
  int book_simp_histset[kNSimpHistSets];
  for (int i=0; i<kNSimpHistSets; i++) book_simp_histset[i] = 0;

  book_simp_histset[  0] = 1;		// all particles
  book_simp_histset[100] = 1;		// pi-
  book_simp_histset[200] = 1;		// pi+
  book_simp_histset[300] = 1;		// pi0
  book_simp_histset[400] = 1;		// protons 
  book_simp_histset[500] = 1;		// neutrons

  book_simp_histset[101] = 1;		// pi-   100 < p < 150
  book_simp_histset[102] = 1;		// pi-   150 < p < 200
  book_simp_histset[103] = 1;		// pi-   200 < p < 250
  book_simp_histset[104] = 1;		// pi-   250 < p < 300
  book_simp_histset[105] = 1;		// pi-   300 < p < 350
  book_simp_histset[106] = 1;		// pi-   350 < p < 400
  book_simp_histset[107] = 1;		// pi-   400 < p < 450
  book_simp_histset[108] = 1;		// pi-   450 < p < 500

  book_simp_histset[121] = 1;		// pi-   0.35 < theta < 0.55
  book_simp_histset[122] = 1;		// pi-   0.55 < theta < 0.75
  book_simp_histset[123] = 1;		// pi-   0.75 < theta < 0.95
  book_simp_histset[124] = 1;		// pi-   0.95 < theta < 1.15
  book_simp_histset[125] = 1;		// pi-   1.15 < theta < 1.35
  book_simp_histset[126] = 1;		// pi-   1.35 < theta < 1.55
  book_simp_histset[127] = 1;		// pi-   1.55 < theta < 1.75
  book_simp_histset[128] = 1;		// pi-   1.75 < theta < 1.95
  book_simp_histset[129] = 1;		// pi-   1.95 < theta < 2.15

  for (int i=0; i<kNSimpHistSets; i++) {
    if (book_simp_histset[i] != 0) {
      sprintf(folder_name,"simp_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fSimp[i] = new SimpHist_t;
      BookSimpHistograms(fHist.fSimp[i],Form("Hist/%s",folder_name));
    }
  }
}

//-----------------------------------------------------------------------------
void TG4ValidationModule::FillSimpHistograms(HistBase_t* Hist, TSimParticle* Simp, SimpData_t* SimpData) {

  SimpHist_t* hist = (SimpHist_t*) Hist;

  hist->fPDGCode[0] ->Fill(Simp->PDGCode());
  hist->fPDGCode[1] ->Fill(Simp->PDGCode());

  const TLorentzVector* v1 = fProton->StartMom();
  const TLorentzVector* v2 = Simp->StartMom();

  double p = v2->P();

  hist->fMomentum[0]->Fill(p);

  double cos_th = (v1->Px()*v2->Px()+v1->Py()*v2->Py()+v1->Pz()*v2->Pz())/v1->P()/v2->P()/(1+1.e-10);
  hist->fCosTheta->Fill(cos_th);

  double th = TMath::ACos(cos_th);
//-----------------------------------------------------------------------------
// cross section plots: normalization of the cross sections
//-----------------------------------------------------------------------------
  // double rho(11.47);			// Ta, g/cm^3
  // double thickness(0.6);		// 6 mm
  // double NAv (6.022e23);		// 

  hist->fMomentum[1]->Fill(p/1000.);	// in GeV/c
  hist->fTheta->Fill(th);
}

//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
void TG4ValidationModule::FillEventHistograms(HistBase_t* Hist) {
//   double            cos_th, xv, yv, rv, zv, p;
//   TLorentzVector    mom;

  EventHist_t* hist = (EventHist_t*) Hist;

  int event_number = GetHeaderBlock()->EventNumber();
  int run_number   = GetHeaderBlock()->RunNumber();

  hist->fEventNumber->Fill(event_number);
  hist->fRunNumber->Fill(run_number);
  hist->fNPiMinus->Fill(fNPiMinus);
  hist->fNPiPlus->Fill(fNPiPlus);
  hist->fNPi0->Fill(fNPi0);
  hist->fNPions->Fill(fNPions);
  hist->fNProtons->Fill(fNProtons);
  hist->fNNeutrons->Fill(fNNeutrons);

}

//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TG4ValidationModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
  RegisterDataBlock("SimpBlock" ,"TSimpBlock" ,&fSimpBlock);
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();

  return 0;
}


//-----------------------------------------------------------------------------
void TG4ValidationModule::FillHistograms() {
//-----------------------------------------------------------------------------
// event histograms
//
// EVT_0: all events
//-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0]);
//-----------------------------------------------------------------------------
// SimParticle histograms
//-----------------------------------------------------------------------------
  int pdg_code(0);

  const TLorentzVector* v1 = fProton->StartMom();

  int nsimp = fSimpBlock->NParticles();
  for (int i=1; i<nsimp; i++) {
    TSimParticle* simp = fSimpBlock->Particle(i);
    SimpData_t*   sd   = fSimpData + i;

    FillSimpHistograms(fHist.fSimp[0],simp,sd);
    pdg_code = simp->PDGCode();

    const TLorentzVector* v2 = simp->StartMom();
//-----------------------------------------------------------------------------
// pi-
//-----------------------------------------------------------------------------
    if      (pdg_code == -211) {	
      FillSimpHistograms(fHist.fSimp[100],simp,sd);
      
      double p = simp->StartMom()->P();
      double cos_th = (v1->Px()*v2->Px()+v1->Py()*v2->Py()+v1->Pz()*v2->Pz())/v1->P()/v2->P()/(1+1.e-10);
      double th     = TMath::ACos(cos_th); // in radians, [0,pi]

      if ((p >= 100) && (p < 150)) FillSimpHistograms(fHist.fSimp[101],simp,sd);
      if ((p >= 150) && (p < 200)) FillSimpHistograms(fHist.fSimp[102],simp,sd);
      if ((p >= 200) && (p < 250)) FillSimpHistograms(fHist.fSimp[103],simp,sd);
      if ((p >= 250) && (p < 300)) FillSimpHistograms(fHist.fSimp[104],simp,sd);
      if ((p >= 300) && (p < 350)) FillSimpHistograms(fHist.fSimp[105],simp,sd);
      if ((p >= 350) && (p < 400)) FillSimpHistograms(fHist.fSimp[106],simp,sd);
      if ((p >= 400) && (p < 450)) FillSimpHistograms(fHist.fSimp[107],simp,sd);
      if ((p >= 450) && (p < 500)) FillSimpHistograms(fHist.fSimp[108],simp,sd);

      if ((th>= 0.35) && (th < 0.55)) FillSimpHistograms(fHist.fSimp[121],simp,sd);
      if ((th>= 0.55) && (th < 0.75)) FillSimpHistograms(fHist.fSimp[122],simp,sd);
      if ((th>= 0.75) && (th < 0.85)) FillSimpHistograms(fHist.fSimp[123],simp,sd);
      if ((th>= 0.95) && (th < 1.15)) FillSimpHistograms(fHist.fSimp[124],simp,sd);
      if ((th>= 1.15) && (th < 1.35)) FillSimpHistograms(fHist.fSimp[125],simp,sd);
      if ((th>= 1.35) && (th < 1.55)) FillSimpHistograms(fHist.fSimp[126],simp,sd);
      if ((th>= 1.55) && (th < 1.75)) FillSimpHistograms(fHist.fSimp[127],simp,sd);
      if ((th>= 1.75) && (th < 1.95)) FillSimpHistograms(fHist.fSimp[128],simp,sd);
      if ((th>= 1.95) && (th < 2.15)) FillSimpHistograms(fHist.fSimp[129],simp,sd);
    }
//-----------------------------------------------------------------------------
// pi+
//-----------------------------------------------------------------------------
    else if (pdg_code ==  211) {
      FillSimpHistograms(fHist.fSimp[200],simp,sd);
    }
//-----------------------------------------------------------------------------
// pi0
//-----------------------------------------------------------------------------
    else if (pdg_code ==  111) {
      FillSimpHistograms(fHist.fSimp[300],simp,sd);
    }
//-----------------------------------------------------------------------------
// protons
//-----------------------------------------------------------------------------
    else if (pdg_code ==  2212) {
      FillSimpHistograms(fHist.fSimp[400],simp,sd);
    }
//-----------------------------------------------------------------------------
// neutrons
//-----------------------------------------------------------------------------
    else if (pdg_code ==  2112) {
      FillSimpHistograms(fHist.fSimp[500],simp,sd);
    }
  }

}



//_____________________________________________________________________________
int TG4ValidationModule::BeginRun() {
  int rn = GetHeaderBlock()->RunNumber();
  TStntuple::Init(rn);
  return 0;
}


//_____________________________________________________________________________
int TG4ValidationModule::Event(int ientry) {

  fSimpBlock->GetEntry(ientry);

  fProton = fSimpBlock->Particle(0);
//-----------------------------------------------------------------------------
// loop over particles
//-----------------------------------------------------------------------------
  fNPiMinus  = 0;
  fNPiPlus   = 0;
  fNPi0      = 0;
  fNNeutrons = 0;
  fNProtons  = 0;

  int nsimp = fSimpBlock->NParticles();
  for (int i=1; i<nsimp; i++) {
    TSimParticle* simp = fSimpBlock->Particle(i);
    int pdg_code = simp->PDGCode();
    if      (pdg_code ==  -211) fNPiMinus++;
    else if (pdg_code ==   211) fNPiPlus++;
    else if (pdg_code ==   111) fNPi0++;
    else if (pdg_code ==  2212) fNProtons++;
    else if (pdg_code ==  2112) fNNeutrons++;
  }

  fNPions = fNPiMinus+fNPiPlus+fNPi0;

  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TG4ValidationModule::Debug() {

//-----------------------------------------------------------------------------
// bit 4: events with NHitsTF > 1
//-----------------------------------------------------------------------------
  // if (GetDebugBit(4) == 1) {
  //   GetHeaderBlock()->Print(Form("NHits(TF) = %5i",fNGenp));
  // }
}

//_____________________________________________________________________________
int TG4ValidationModule::EndJob() {
  printf("----- end job: ---- %s\n",GetName());
  return 0;
}

//_____________________________________________________________________________
void TG4ValidationModule::Test001() {
}

