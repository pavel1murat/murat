//////////////////////////////////////////////////////////////////////////////
// use of TStnTrack::fTmp:
//
// Tmp(0) : corrected momentum at the tracker front - not yet
// 
// use of debug bits: bits 0-2 are reserved
//  0  : all events
//  1  : passed events
//  2  : rejected events
//
//  3  : events with NElePos > 1
// 
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
#include "Stntuple/geom/TDisk.hh"
#include "Stntuple/obj/TStnNode.hh"
#include "Stntuple/val/stntuple_val_functions.hh"
#include "Offline/DataProducts/inc/VirtualDetectorId.hh"

#include "murat/ana/TSimpAnaModule.hh"

ClassImp(murat::TSimpAnaModule)

namespace murat {
//-----------------------------------------------------------------------------
TSimpAnaModule::TSimpAnaModule(const char* name, const char* title): TAnaModule(name,title) {
//-----------------------------------------------------------------------------
// MC truth: define which MC particle to consider as signal
// 2:conversionGun, 28:StoppedParticleReactionGun - see 
//-----------------------------------------------------------------------------
  fMinNStrawHits = 20;
}

//-----------------------------------------------------------------------------
TSimpAnaModule::~TSimpAnaModule() {
}


//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TSimpAnaModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
  RegisterDataBlock("GenpBlock"           ,"TGenpBlock"       ,&fGenpBlock);
  RegisterDataBlock("SimpBlock"           ,"TSimpBlock"       ,&fSimpBlock);
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();

  return 0;
}


//-----------------------------------------------------------------------------
void TSimpAnaModule::BookHistograms() {

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
  book_event_histset[ 1] = 1;		// events with N(electrons) > 0
  book_event_histset[ 2] = 1;		// events with N(muons)     > 0
  book_event_histset[ 3] = 1;		// events with N(e+ + e-)   > 1

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
// book Genp histograms
//-----------------------------------------------------------------------------
  int book_genp_histset[kNGenpHistSets];
  for (int i=0; i<kNGenpHistSets; i++) book_genp_histset[i] = 0;

  book_genp_histset[0] = 1;		// all particles

  for (int i=0; i<kNGenpHistSets; i++) {
    if (book_genp_histset[i] != 0) {
      sprintf(folder_name,"gen_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fGenp[i] = new GenpHist_t;
      BookGenpHistograms(fHist.fGenp[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book simp histograms
//-----------------------------------------------------------------------------
  int book_simp_histset[kNSimpHistSets];
  for (int i=0; i<kNSimpHistSets; i++) book_simp_histset[i] = 0;

  book_simp_histset[  0] = 1;		// all particles
  book_simp_histset[  1] = 1;		// all particles with N(straw hits) > fMinNStrawHits
  book_simp_histset[ 10] = 1;		// electrons with N(straw hits) > fMinNStrawHits
  book_simp_histset[ 20] = 1;		// positrons with N(straw hits) > fMinNStrawHits
  book_simp_histset[ 30] = 1;		// muons     with N(straw hits) > fMinNStrawHits
  book_simp_histset[ 40] = 1;		// pions     with N(straw hits) > fMinNStrawHits
  book_simp_histset[ 41] = 1;		// pions     with with the weight of survival probability

  for (int i=0; i<kNSimpHistSets; i++) {
    if (book_simp_histset[i] != 0) {
      sprintf(folder_name,"sim_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fSimp[i] = new SimpHist_t;
      BookSimpHistograms(fHist.fSimp[i],Form("Hist/%s",folder_name));
    }
  }
}

//_____________________________________________________________________________
void TSimpAnaModule::FillHistograms() {

//-----------------------------------------------------------------------------
// event histograms
//-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0],&fEvtPar);

  if (fNElectrons >= 1) FillEventHistograms(fHist.fEvent[1],&fEvtPar);
  if ((fNMuons[1] >= 1) and (fNMuons[0] >= 2)) FillEventHistograms(fHist.fEvent[2],&fEvtPar);
  if ((fNEle40 > 1    ) and (fNPos40 > 0    )) FillEventHistograms(fHist.fEvent[3],&fEvtPar);
//-----------------------------------------------------------------------------
// Simp histograms
//-----------------------------------------------------------------------------
  TSimParticle* simp;
  SimpData_t    sd;

  sd.fParent = nullptr;
  
  for (int i=0; i<fEvtPar.fNSimp; i++) {
    simp = fSimpBlock->Particle(i);
    FillSimpHistograms(fHist.fSimp[0],simp,&sd);
    int pdg_code = simp->PDGCode();
    if (simp->NStrawHits() > fMinNStrawHits) {
      FillSimpHistograms(fHist.fSimp[1],simp,&sd);
      if (pdg_code == 11) {
	FillSimpHistograms(fHist.fSimp[10],simp,&sd);
      }
      else if (pdg_code == -11) {
	FillSimpHistograms(fHist.fSimp[20],simp,&sd);
      }
      else if (abs(pdg_code) == 13) {
	FillSimpHistograms(fHist.fSimp[30],simp,&sd);
      }
      else if (abs(pdg_code) == 211) {
	FillSimpHistograms(fHist.fSimp[40],simp,&sd);
	FillSimpHistograms(fHist.fSimp[41],simp,&sd,fEvtPar.fPionSurvProb);
      }
    }
  }
//-----------------------------------------------------------------------------
// fill GENP histograms
// GEN_0: all particles
//-----------------------------------------------------------------------------
  TGenParticle* genp;
  for (int i=0; i<fEvtPar.fNGenp; i++) {
    genp = fGenpBlock->Particle(i);
    FillGenpHistograms(fHist.fGenp[0],genp);
  }
}



//-----------------------------------------------------------------------------
// 2014-04-30: it looks that reading the straw hits takes a lot of time - 
//              turn off by default by commenting it out
//-----------------------------------------------------------------------------
int TSimpAnaModule::Event(int ientry) {

  fGenpBlock->GetEntry(ientry);
  fSimpBlock->GetEntry(ientry);
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
  fEvtPar.fNGenp            = fGenpBlock->NParticles();
  fEvtPar.fSimp             = nullptr;
  fEvtPar.fPartE            = -1.;
  fEvtPar.fNSimp            = fSimpBlock->NParticles();
  fEvtPar.fPionSurvProb     = 1.;
//-----------------------------------------------------------------------------
// luminosity weight
//-----------------------------------------------------------------------------
//  fLumWt = GetHeaderBlock()->LumWeight();
//-----------------------------------------------------------------------------
// figure MC truth
//-----------------------------------------------------------------------------
  for (int i=fEvtPar.fNSimp-1; i>=0; i--) {
    TSimParticle* simp = fSimpBlock->Particle(i);
    int pdg_code       = simp->PDGCode();
    int generator_id   = simp->GeneratorID();         // MC process code

    if ((pdg_code == fPDGCode) and (generator_id == fMCProcessCode)) {
      fEvtPar.fSimp  = simp;
      fEvtPar.fPartE = simp->StartMom()->Energy();

      if (abs(pdg_code) == 211) {
        fEvtPar.fPionSurvProb = simp->SurvivalProb();
      }
      break;
    }
  }
//-----------------------------------------------------------------------------
// calculate DIO weights once per event
//-----------------------------------------------------------------------------
  if (fEvtPar.fPartE > 0) {
    fEvtPar.fDioLOWt = TStntuple::DioWeightAl   (fEvtPar.fPartE);
    fEvtPar.fDioLLWt = TStntuple::DioWeightAl_LL(fEvtPar.fPartE);
  }
//-----------------------------------------------------------------------------
// may want to revisit the definition of fSimp in the future
//-----------------------------------------------------------------------------
  fSimPar.fParticle = fEvtPar.fSimp;
  fSimPar.fTFront   = NULL;
  fSimPar.fTMid     = NULL;
  fSimPar.fTBack    = NULL;
  // fSimPar.fGenp     = fEvtPar.fParticle;
//-----------------------------------------------------------------------------
// process virtual detectors - for fSimp need parameters at tracker entrance
// use the first fit
//-----------------------------------------------------------------------------
  fSimPar.fTFront   = NULL;
  fSimPar.fTMid     = NULL;
  fSimPar.fTBack    = NULL;

  fNMuons[0]  = 0;
  fNMuons[1]  = 0;
  fNElectrons = 0;
  fNEle40     = 0;
  fNPos40     = 0;

  for (int i=0; i<fEvtPar.fNSimp; i++) {
    TSimParticle* simp = fSimpBlock->Particle(i);
    int pdg_code       = simp->PDGCode();
    float mom          = simp->StartMom()->P();
    float time         = simp->StartPos()->T();
    if ((abs(pdg_code) == 13) or (abs(pdg_code) == 211)) {
      if ((simp->NStrawHits() >= fMinNStrawHits) and (time > 500)) {
	if (mom > 60) fNMuons[0] += 1;
	if (mom > 90) fNMuons[1] += 1;
	
      }
    }
    else if (pdg_code == 11) {
      if (simp->NStrawHits() >= fMinNStrawHits) {

	if ((mom > 90) and (mom < 130) and (time > 640)) {
	  fNElectrons += 1;
	}
      }
    }

    if (mom > 40) {
      if (pdg_code ==  11) fNEle40 += 1;
      if (pdg_code == -11) fNPos40 += 1;
    }
  }

  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TSimpAnaModule::Debug() {

  if (GetDebugBit(3) == 1) {
    if ((fNEle40 > 0) and (fNPos40 > 0)) GetHeaderBlock()->Print(Form("fNEle40, fNPos40 = %i %i",fNEle40,fNPos40));
  }

  // TStnTrack* trk;
  // TrackPar_t* tp;
  // int ntrk = fTrackBlock->NTracks();

}

//_____________________________________________________________________________
int TSimpAnaModule::EndJob() {
  printf("----- end job: ---- %s\n",GetName());
  return 0;
}

}
