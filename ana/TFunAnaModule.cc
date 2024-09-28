//////////////////////////////////////////////////////////////////////////////
// use of TStnTrack::fTmp:
//
// Tmp(0) : corrected momentum at the tracker front - not yet
// 
// use of debug bits: 
//  0  : one line per track
//  1  : passed events
//  2  : rejected events
// 
//  3  : events with set C tracks and 70mm < |dx|  < 90 mm
//  4  : events with DpF > 1 MeV : obviously, misreconstructed ones
//  5  : events with N(tracks) > 1
//  6  : events trk_41 with 0.8< E/P < 1.1 - tracks missed by CalPatRec
//  7  : events (muo) with LogLHRCal >   20
//  8  : events (ele) with LogLHRCal < - 20
//  9  : events (muo) with 0.42 < E/P < 0.46
// 10  : events (muo) with Set C track with ECL > 80 MeV
// 28  : Set C DEM tracks with E/P > 1.1
// 29  : TRK_19 (Set C DEM tracks with a cluster) and LLHR(cal) < 0
// 31  : EVT_6 events with ce_costh > 0.8 
// 32  : TRK_1 events with chi2tcm > 100. 
// 33  : DU < -80mm - study edge effects
// 34  : EVT_7: events with E_CL > 60 and no CalPatRec tracks 
// 35  : TRK_1: events with P > 106 MeV/c - misreconstruction
// 36  : TRK_23 events with P < 80: odd misidentified muons - turned out to be DIO electrons
// 37  : TRK_26 LLHR_CAL > 5
// 38  : EVT_7: events with E_CL > 60 and no CalPatRec tracks and TrkPatRec
// 39  : trk_1: events with |SIN_TC| > 0.6
// 40  : trk_0: events with P > 60 MeV/c
// 41  : events with ntrk > 1
//
// 3 different ID : 
///////////////////////////////////////////////////////////////////////////////
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TEnv.h"
#include "TSystem.h"
#include "TDirectory.h"

#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/alg/TStntuple.hh"
#include "Stntuple/obj/TStnNode.hh"
#include "Stntuple/val/stntuple_val_functions.hh"
#include "Offline/DataProducts/inc/VirtualDetectorId.hh"

//------------------------------------------------------------------------------
// Mu2e offline includes
//-----------------------------------------------------------------------------
#include "murat/ana/TFunAnaModule.hh"

ClassImp(murat::TFunAnaModule)

namespace murat {
//-----------------------------------------------------------------------------
TFunAnaModule::TFunAnaModule(const char* name, const char* title): TAnaModule(name,title)
{
  fTrackBlockName         = "TrackBlock";
//-----------------------------------------------------------------------------
// initialize Track ID
// 0: SetC  1: TrkQual>0.1 2:TrkQual>0.4
// what about number of hits ? - 3: no cuts on the number of hits
//-----------------------------------------------------------------------------
}

//-----------------------------------------------------------------------------
TFunAnaModule::~TFunAnaModule() {
}


//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TFunAnaModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
  RegisterDataBlock("SimpBlock"            ,"TSimpBlock"       ,&fSimpBlock   );
  RegisterDataBlock("SpmcBlock"            ,"TStepPointMCBlock",&fVDetBlock   );
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();

  return 0;
}

//-----------------------------------------------------------------------------
  void TFunAnaModule::BookFunHistograms(FunHist_t* Hist, const char* Folder) {
  //     char name [200];
  //     char title[200];
  //-----------------------------------------------------------------------------
  //  
  //-----------------------------------------------------------------------------
    HBook1F(Hist->fTofVD13        ,"tof_vd13",Form("%s: TOF to VD13"                 ,Folder), 200,   0, 100,Folder);
    HBook2F(Hist->fTofVD13VsCosth ,"tvd13_vs_costh",Form("%s: TOF to VD13 vs Cos(TH)",Folder), 100,  -1,   1,200, 0,  100,Folder);
    HBook2F(Hist->fTofVD13VsZ0    ,"tvd13_vs_z0"   ,Form("%s: TOF to VD13 vs Z0"     ,Folder), 120,5100,6300,200, 0,  100,Folder);
    HBook2F(Hist->fTofVD13VsPhi   ,"tvd13_vs_phi"  ,Form("%s: TOF to VD13 vs Phi"    ,Folder), 126,   0, 6.3,200, 0,  100,Folder);
    HBook2F(Hist->fTofVD13VsRho   ,"tvd13_vs_rho"  ,Form("%s: TOF to VD13 vs rho^p"  ,Folder), 100,   0, 1000,200, 0,  100,Folder);
  }

//_____________________________________________________________________________
void TFunAnaModule::BookHistograms() {

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
// book simp histograms
//-----------------------------------------------------------------------------
  int book_simp_histset[kNSimpHistSets];
  for (int i=0; i<kNSimpHistSets; i++) book_simp_histset[i] = 0;

  book_simp_histset[ 0] = 1;		// all events

  for (int i=0; i<kNSimpHistSets; i++) {
    if (book_simp_histset[i] != 0) {
      sprintf(folder_name,"sim_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fSimp[i] = new SimpHist_t;
      BookSimpHistograms(fHist.fSimp[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book fun (temporary, module specific) histograms
//-----------------------------------------------------------------------------
  int book_fun_histset[kNFunHistSets];
  for (int i=0; i<kNFunHistSets; i++) book_fun_histset[i] = 0;

  book_fun_histset[0] = 1;		// all particles
  book_fun_histset[1] = 1;		// Pz > 0 particles
  book_fun_histset[2] = 1;		// Pz < 0 particles

  for (int i=0; i<37; i++) {
    book_fun_histset[100+i] = 1;            // per foil
  }

  for (int i=0; i<kNFunHistSets; i++) {
    if (book_fun_histset[i] != 0) {
      sprintf(folder_name,"fun_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fFun[i] = new FunHist_t;
      BookFunHistograms(fHist.fFun[i],Form("Hist/%s",folder_name));
    }
  }
}

//-----------------------------------------------------------------------------
void TFunAnaModule::FillFunHistograms(FunHist_t* Hist) {
  
  Hist->fTofVD13->Fill(fEvtPar.fTofVD13);

  double px = fEvtPar.fSimp->StartMom()->Px();
  double py = fEvtPar.fSimp->StartMom()->Py();

  double costh  = fEvtPar.fSimp->StartMom()->Pz()/fEvtPar.fSimp->StartMom()->P();
  Hist->fTofVD13VsCosth->Fill(costh,fEvtPar.fTofVD13);
  Hist->fTofVD13VsZ0->Fill   (fEvtPar.fSimp->StartPos()->Z(),fEvtPar.fTofVD13);

  double phi    = atan2(fEvtPar.fSimp->StartMom()->Py(),fEvtPar.fSimp->StartMom()->Px())+M_PI;
  Hist->fTofVD13VsPhi->Fill  (phi,fEvtPar.fTofVD13);

  double x      = fEvtPar.fSimp->StartPos()->X()+3904.;
  double y      = fEvtPar.fSimp->StartPos()->Y();

  double nux =  py/fEvtPar.fSimp->StartMom()->Pt();
  double nuy = -px/fEvtPar.fSimp->StartMom()->Pt();
  // assume the Bfield around the target of 1.5 T, r in mm
  double r   = fEvtPar.fSimp->StartMom()->Pt()/3*10/1.5;

  double rhox = x-nux*r;
  double rhoy = y-nuy*r;

  double rho = sqrt(rhox*rhox+rhoy*rhoy)+r;
  
  Hist->fTofVD13VsRho->Fill  (rho,fEvtPar.fTofVD13);
}

//_____________________________________________________________________________
void TFunAnaModule::FillHistograms() {

//-----------------------------------------------------------------------------
// event histograms
//-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0],&fEvtPar);
//-----------------------------------------------------------------------------
// Simp histograms
//-----------------------------------------------------------------------------
  if (fSimp) {
    SimpData_t sd; sd.fParent = nullptr;
    FillSimpHistograms(fHist.fSimp[0],fSimp,&sd);
  }
//-----------------------------------------------------------------------------
// fill fun histograms
// fun_0: all particles
//-----------------------------------------------------------------------------
  FillFunHistograms(fHist.fFun[0]);
  if (fEvtPar.fSimp->StartMom()->Pz() > 0) FillFunHistograms(fHist.fFun[1]);
  else                                     FillFunHistograms(fHist.fFun[2]);

  // per-foil histograms. there are a few events with v0 = 3071 - why ?
  for (int i=0; i<37; i++) {
    int foil = fEvtPar.fSimp->StartVolumeIndex()-3072;
    if ((foil >= 0) and (foil < 37)) {
      FillFunHistograms(fHist.fFun[100+foil]);
    }
  }
}


//-----------------------------------------------------------------------------
// 2014-04-30: it looks that reading the straw hits takes a lot of time - 
//              turn off by default by commenting it out
//-----------------------------------------------------------------------------
int TFunAnaModule::Event(int ientry) {

  //  TDiskCalorimeter::GeomData_t disk_geom;

  fSimpBlock->GetEntry(ientry);
  fVDetBlock->GetEntry(ientry);     int nvdhits = fVDetBlock->NStepPoints();
//-----------------------------------------------------------------------------
// at some point, the TVdetBlock class got obsolete, and now the virtual detector 
// hits are stored in TStepPointMCBlock 
// dont' try to read it as that would fail
//-----------------------------------------------------------------------------
  fEvtPar.fDioLOWt          = 1.;
  fEvtPar.fDioLLWt          = 1.;
  fEvtPar.fNCrvClusters     = -1;
  fEvtPar.fNCrvPulses       = -1;
  fEvtPar.fNCrvCoincidences = -1;
  fEvtPar.fNStrawHits       = GetHeaderBlock()->NStrawHits();
  fEvtPar.fNComboHits       = GetHeaderBlock()->NComboHits();
  fEvtPar.fNGenp            = 0;
  fEvtPar.fNSimp            = fSimpBlock->NParticles();
  fEvtPar.fSimp             = NULL;
  fEvtPar.fPartE            = -1.;
  fEvtPar.fPionSurvProb     =  1.;
  fEvtPar.fTofVD13          = -1.;
//-----------------------------------------------------------------------------
// luminosity weight
//-----------------------------------------------------------------------------
  // fLumWt = GetHeaderBlock()->LumWeight();
//-----------------------------------------------------------------------------
// figure out the MC truth
//-----------------------------------------------------------------------------
  for (int i=fEvtPar.fNSimp-1; i>=0; i--) {
    TSimParticle* simp = fSimpBlock->Particle(i);
    int pdg_code       = simp->PDGCode();
    int generator_id   = simp->GeneratorID();         // MC process code

    if ((pdg_code == fPDGCode) and (generator_id == fMCProcessCode)) {
      fEvtPar.fSimp  = simp;
      fEvtPar.fPartE = simp->StartMom()->Energy();
      for (int iv=0; iv<nvdhits; iv++) {
        TStepPointMC* vdh = fVDetBlock->StepPointMC(iv);
        if ((vdh->SimID() == (int) simp->GetUniqueID()) and (vdh->VolumeID() == 13)) {
          fEvtPar.fTofVD13 = vdh->Time()-simp->StartPos()->T();
        }
      }
      break;
    }
  }
//-----------------------------------------------------------------------------
// may want to revisit the definition of fSimp in the future
//-----------------------------------------------------------------------------
  fSimp             = fEvtPar.fSimp;
  fSimPar.fParticle = fSimp;
  fSimPar.fTFront   = NULL;
  fSimPar.fTMid     = NULL;
  fSimPar.fTBack    = NULL;
//-----------------------------------------------------------------------------
// process virtual detectors - for fSimp need parameters at tracker entrance
// use the first fit
//-----------------------------------------------------------------------------
  fSimPar.fTFront   = NULL;
  fSimPar.fTMid     = NULL;
  fSimPar.fTBack    = NULL;
  for (int i=0; i<nvdhits; i++) {
    TStepPointMC* vdhit = fVDetBlock->StepPointMC(i);
    if (vdhit->SimID() == (int) fEvtPar.fSimp->GetUniqueID()) {
//-----------------------------------------------------------------------------
// 'signal' particle'
//-----------------------------------------------------------------------------
      if ((vdhit->VolumeID() == 13) || (vdhit->VolumeID() == 14)) {
	if (fDirection*vdhit->Mom()->Z() > 0) {
	  if (fSimPar.fTFront == 0) fSimPar.fTFront = vdhit;
	}
      }
      else if ((vdhit->VolumeID() == 11) || (vdhit->VolumeID() == 12)) {
	if (fDirection*vdhit->Mom()->Z() > 0) {
	  if (fSimPar.fTMid == 0) fSimPar.fTMid = vdhit;
	}
      }
      else if (vdhit->VolumeID() == mu2e::VirtualDetectorId::TT_Back) {
	if (fDirection*vdhit->Mom()->Z() > 0) {
	  if (fSimPar.fTBack == NULL) fSimPar.fTBack = vdhit;
	}
      }
    }
  }

//-----------------------------------------------------------------------------
// loop over tracks and calculate needed parameters
//-----------------------------------------------------------------------------
  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TFunAnaModule::Debug() {

  if (GetDebugBit(41) == 1) {
    GetHeaderBlock()->Print(Form("bit_041: ntrk = %i",1));
  }
}

//_____________________________________________________________________________
int TFunAnaModule::EndJob() {
  printf("----- end job: ---- %s\n",GetName());
  return 0;
}

//_____________________________________________________________________________
void TFunAnaModule::Test001() {
}

}
