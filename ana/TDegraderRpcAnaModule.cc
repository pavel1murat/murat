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
//  3  : events with good tracks P > 70 and T0 > 300
//  5  : events with N(tracks) > 1
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

//------------------------------------------------------------------------------
// Mu2e offline includes
//-----------------------------------------------------------------------------

#include "murat/ana/TDegraderRpcAnaModule.hh"

ClassImp(murat::TDegraderRpcAnaModule)

namespace murat {
//-----------------------------------------------------------------------------
TDegraderRpcAnaModule::TDegraderRpcAnaModule(const char* name, const char* title): 
  TAnaModule(name,title) { 
  fTrackBlockName         = "TrackBlock";
  fTrackStrawHitBlockName = "TrackStrawHitBlock";
  fTrackNumber.Set(100);

  //  fDiskCalorimeter = new TDiskCalorimeter();
  fCalorimeterType = 2;

  fMinT0           = 700; 
  //  fMinETrig        = 50.;               // MeV
					// track-cluster matching timing cut
  // fMinDtTcm        = -5.;
  // fMaxDtTcm        =  8.;
//-----------------------------------------------------------------------------
// initialize Track ID
// 0: SetC  1: TrkQual>0.1 2:TrkQual>0.4
// what about number of hits ? - 3: no cuts on the number of hits
//-----------------------------------------------------------------------------
  fNID             = 2;

  fTrackID[0]      = new TStnTrackID();

  int mask0 = TStnTrackID::kNActiveBit | TStnTrackID::kChi2DofBit | TStnTrackID::kMomErrBit ;
  fTrackID[0]->SetMaxChi2Dof(3.0 );
  fTrackID[0]->SetMinNActive(25  );
  fTrackID[0]->SetMaxMomErr (0.25);

  mask0    |= TStnTrackID::kTanDipBit;
  fTrackID[0]->SetMinTanDip (0.6);
  fTrackID[0]->SetMaxTanDip (1.0);

  mask0    |= TStnTrackID::kD0Bit;
  fTrackID[0]->SetMinD0 (-100);         // mm
  fTrackID[0]->SetMaxD0 ( 100);

  fTrackID[0]->SetUseMask(mask0);
//-----------------------------------------------------------------------------
// tight ID
//-----------------------------------------------------------------------------
  fTrackID[1]      = new TStnTrackID();

  int mask1 = TStnTrackID::kNActiveBit | TStnTrackID::kChi2DofBit | TStnTrackID::kMomErrBit ;
  fTrackID[1]->SetMaxChi2Dof(3.0 );
  fTrackID[1]->SetMinNActive(30  );     // 30 ... cut tighter on the number of hits
  fTrackID[1]->SetMaxMomErr (0.25);

  mask1    |= TStnTrackID::kTanDipBit;
  fTrackID[1]->SetMinTanDip (0.70);     // prev: 0.65 
  fTrackID[1]->SetMaxTanDip (1.20);

  mask1    |= TStnTrackID::kD0Bit;
  fTrackID[1]->SetMinD0 ( -50);         // mm
  fTrackID[1]->SetMaxD0 ( 100);

  fTrackID[1]->SetUseMask(mask1);

  fBestID           = 0;

  for (int i=0; i<20; i++) {
    TrackPar_t* tp  = fTrackPar+i;
    tp->fTrackID[0] = TAnaModule::fTrackID[0];
    tp->fTrackID[1] = TAnaModule::fTrackID[1];
  }

//-----------------------------------------------------------------------------
// MC truth: define which MC particle to consider as signal
// 2:conversionGun, 28:StoppedParticleReactionGun - see 
//-----------------------------------------------------------------------------
  fDirection        = 1;

  fDnMax            = 15;
//-----------------------------------------------------------------------------
// this is redefined by the dataset catalog anyway
//-----------------------------------------------------------------------------
  fPDGCode          = 22;
  fMCProcessCode    = 178;   // mu2e::ProcessCode::mu2eExternalRPC;
}

//-----------------------------------------------------------------------------
TDegraderRpcAnaModule::~TDegraderRpcAnaModule() {
  for (int i=0; i<fNID; i++) delete fTrackID[i];
}


//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TDegraderRpcAnaModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
  RegisterDataBlock(fTrackBlockName.Data()        ,"TStnTrackBlock"     ,&fTrackBlock  );
  // RegisterDataBlock(fTrackStrawHitBlockName.Data(),"TTrackStrawHitBlock",&fTrackStrawHitBlock);

  RegisterDataBlock("KSFBlock"        ,"TStnTrackSeedBlock"  ,&fTrackSeedBlock  );
  RegisterDataBlock("HelixBlock"      ,"TStnHelixBlock"      ,&fHelixBlock      );
  RegisterDataBlock("TimeClusterBlock","TStnTimeClusterBlock",&fTimeClusterBlock);

  RegisterDataBlock("ClusterBlock" ,"TStnClusterBlock" ,&fClusterBlock );
  RegisterDataBlock("CalDataBlock" ,"TCalDataBlock"    ,&fCalDataBlock );
  // RegisterDataBlock("StrawHitBlock","TStrawHitBlock"   ,&fStrawHitBlock);
  RegisterDataBlock("GenpBlock"    ,"TGenpBlock"       ,&fGenpBlock    );
  RegisterDataBlock("SimpBlock"    ,"TSimpBlock"       ,&fSimpBlock    );
  RegisterDataBlock("SpmcBlockVDet","TStepPointMCBlock",&fVDetBlock    );
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();

  return 0;
}

//-----------------------------------------------------------------------------
void TDegraderRpcAnaModule::BookDRpcHistograms(DRpcHist_t* Hist, const char* Folder) {
  HBook1F(Hist->fNHitsVD09 ,"nvd09"   ,Form("%s: N(hits) VD09"  ,Folder),  20, 0,   20,Folder);
  HBook1F(Hist->fSMomVD09[0]  ,"smvd09_0",Form("%s: sum mom VD09[0]"  ,Folder), 200, 50, 150,Folder);
  HBook1F(Hist->fSMomVD09[1]  ,"smvd09_1",Form("%s: sum mom VD09[1]"  ,Folder), 500, 90, 140,Folder);
  HBook1F(Hist->fNHitsVD10 ,"nvd10" ,Form("%s: N(hits) VD10"  ,Folder),  20, 0,   20,Folder);
  HBook1F(Hist->fSMomVD10[0]  ,"smvd10_0",Form("%s: sum mom VD10[0]"  ,Folder), 200, 50, 150,Folder);
  HBook1F(Hist->fSMomVD10[1]  ,"smvd10_1",Form("%s: sum mom VD10[1]"  ,Folder), 500, 90, 140,Folder);
  HBook1F(Hist->fNHitsVD13 ,"nvd13" ,Form("%s: N(hits) VD13"  ,Folder),  20, 0,   20,Folder);
  HBook1F(Hist->fSMomVD13[0]  ,"smvd13_0",Form("%s: sum mom VD13[0]"  ,Folder), 200, 50, 150,Folder);
  HBook1F(Hist->fSMomVD13[1]  ,"smvd13_1",Form("%s: sum mom VD13[1]"  ,Folder), 500, 90, 140,Folder);

  HBook2F(Hist->fSMomVD13VsSinTh,"smvd13_vs_sinth",Form("%s: sum mom VD13 vs sinth",Folder), 100,-0.1,0.1,200,50.,150.,Folder);

  HBook1F(Hist->fCPath           ,"cpath",Form("%s: cpath"  ,Folder), 200,  0,  2,Folder);
  HBook2F(Hist->fSMomVD13VsCPath ,"smvd13_vs_cpath",Form("%s: smvd13 vs conv path"  ,Folder), 200,  0,  2,200,50,150,Folder);
}

//-----------------------------------------------------------------------------
void TDegraderRpcAnaModule::BookTimeClusterHistograms(TimeClusterHist_t* Hist, const char* Folder) {
  HBook1F(Hist->fNsh  ,"nsh"  ,Form("%s: N(single straw hits)",Folder), 200, 0,  200,Folder);
  HBook1F(Hist->fNch  ,"nch"  ,Form("%s: N(combo hits)"       ,Folder), 200, 0,  200,Folder);
  HBook1F(Hist->fT0   ,"t0"   ,Form("%s: T0"                  ,Folder), 200, 0, 2000,Folder);
  HBook1F(Hist->fT0Err,"t0err",Form("%s: T0Err"               ,Folder), 200, 0, 2000,Folder);
}

//-----------------------------------------------------------------------------
void TDegraderRpcAnaModule::BookT2Histograms(T2Hist_t* Hist, const char* Folder) {
  HBook1F(Hist->fSMom[0] ,"smom_0",Form("%s: P1+P2"                ,Folder), 500,   0,  500,Folder);
  HBook1F(Hist->fSMom[1] ,"smom_1",Form("%s: P1+P2"                ,Folder), 500, 100,  150,Folder);
  HBook1F(Hist->fQ    ,"q"   ,Form("%s: Q1+Q2"                ,Folder),  10, -5,   5,Folder);
}

//_____________________________________________________________________________
void TDegraderRpcAnaModule::BookHistograms() {

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
  book_event_histset[ 1] = 1;	        // events with two simparticles and hits

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
// book helix histograms
//-----------------------------------------------------------------------------
  int book_helix_histset[kNHelixHistSets];
  for (int i=0; i<kNHelixHistSets; i++) book_helix_histset[i] = 0;

  book_helix_histset[ 0] = 1;		// all events

  for (int i=0; i<kNHelixHistSets; i++) {
    if (book_helix_histset[i] != 0) {
      sprintf(folder_name,"hel_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fHelix[i] = new HelixHist_t;
      BookHelixHistograms(fHist.fHelix[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book simp histograms
//-----------------------------------------------------------------------------
  int book_simp_histset[kNSimpHistSets];
  for (int i=0; i<kNSimpHistSets; i++) book_simp_histset[i] = 0;

  book_simp_histset[ 0] = 1;		// all events
  book_simp_histset[ 1] = 1;		// events with two reconstructable simparticles

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
// book track histograms
//-----------------------------------------------------------------------------
  int book_track_histset[kNTrackHistSets];
  for (int i=0; i<kNTrackHistSets; i++) book_track_histset[i] = 0;

  book_track_histset[  0] = 0;		// all tracks

  for (int i=0; i<kNTrackHistSets; i++) {
    if (book_track_histset[i] != 0) {
      sprintf(folder_name,"trk_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fTrack[i] = new TrackHist_t;
      BookTrackHistograms(fHist.fTrack[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book track histograms
//-----------------------------------------------------------------------------
  int book_track_id_histset[kNTrackIDHistSets];
  for (int i=0; i<kNTrackIDHistSets; i++) book_track_id_histset[i] = 0;

  book_track_id_histset[  0] = 1;		// all tracks

  for (int i=0; i<kNTrackIDHistSets; i++) {
    if (book_track_id_histset[i] == 0)                      continue;
    sprintf(folder_name,"tid_%i",i);
    fol = (TFolder*) hist_folder->FindObject(folder_name);
    if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
    fHist.fTrackID[i] = new TStnTrackID::Hist_t;
    BookTrackIDHistograms(fHist.fTrackID[i],Form("Hist/%s",folder_name));
  }
//-----------------------------------------------------------------------------
// book two-track histograms
//-----------------------------------------------------------------------------
  int book_t2_histset[kNT2HistSets];
  for (int i=0; i<kNT2HistSets; i++) book_t2_histset[i] = 0;

  book_t2_histset[  0] = 1;		// all 2-track combinations
  book_t2_histset[  1] = 1;		// Q = 0

  for (int i=0; i<kNT2HistSets; i++) {
    if (book_t2_histset[i] == 0)                      continue;
    sprintf(folder_name,"t2_%i",i);
    fol = (TFolder*) hist_folder->FindObject(folder_name);
    if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
    fHist.fT2[i] = new T2Hist_t;
    BookT2Histograms(fHist.fT2[i],Form("Hist/%s",folder_name));
  }
//-----------------------------------------------------------------------------
// book time cluster histograms
//-----------------------------------------------------------------------------
  int book_tc_histset[kNTimeClusterHistSets];
  for (int i=0; i<kNTimeClusterHistSets; i++) book_tc_histset[i] = 0;

  book_tc_histset[0] = 1;		// all clusters

  for (int i=0; i<kNTimeClusterHistSets; i++) {
    if (book_tc_histset[i] == 0)                            continue;
    sprintf(folder_name,"tc_%i",i);
    fol = (TFolder*) hist_folder->FindObject(folder_name);
    if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
    fHist.fTimeCluster[i] = new TimeClusterHist_t;
    BookTimeClusterHistograms(fHist.fTimeCluster[i],Form("Hist/%s",folder_name));
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
// book special histograms
//-----------------------------------------------------------------------------
  int book_drpc_histset[kNDRpcHistSets];
  for (int i=0; i<kNDRpcHistSets; i++) book_drpc_histset[i] = 0;

  book_drpc_histset[  0] = 1;		// all  events
  book_drpc_histset[  1] = 1;		//      events with two particles
  book_drpc_histset[  2] = 1;		//      events with two particles at vd13

  for (int i=0; i<kNDRpcHistSets; i++) {
    if (book_drpc_histset[i] == 0)                        continue;
    sprintf(folder_name,"drpc_%i",i);
    fol = (TFolder*) hist_folder->FindObject(folder_name);
    if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
    fHist.fDRpc[i] = new DRpcHist_t;
    BookDRpcHistograms(fHist.fDRpc[i],Form("Hist/%s",folder_name));
  }
}

//-----------------------------------------------------------------------------
void    TDegraderRpcAnaModule::FillDRpcHistograms    (DRpcHist_t*    Hist, 
                                                      double         Weight) {
  Hist->fNHitsVD09->Fill(fNHitsVD09,Weight);
  Hist->fSMomVD09[0]->Fill(fSMomVD09,Weight);
  Hist->fSMomVD09[1]->Fill(fSMomVD09,Weight);

  Hist->fNHitsVD10->Fill(fNHitsVD10,Weight);
  Hist->fSMomVD10[0]->Fill(fSMomVD10,Weight);
  Hist->fSMomVD10[1]->Fill(fSMomVD10,Weight);

  Hist->fNHitsVD13->Fill(fNHitsVD13,Weight);
  Hist->fSMomVD13[0]->Fill(fSMomVD13,Weight);
  Hist->fSMomVD13[1]->Fill(fSMomVD13,Weight);

  Hist->fCPath->Fill(fCPath,Weight);
  
  if (fEvtPar.fSimp != nullptr) {
    float sinth = fEvtPar.fSimp->StartMom()->Pz()/(fEvtPar.fSimp->StartMom()->P()+1.e-12);
    Hist->fSMomVD13VsSinTh->Fill(sinth,fSMomVD13,Weight);

    Hist->fSMomVD13VsCPath->Fill(fCPath,fSMomVD13,Weight);
  }
}

//-----------------------------------------------------------------------------
void    TDegraderRpcAnaModule::FillT2Histograms(T2Hist_t* Hist, T2Par_t* T2Par, double Weight) {
  Hist->fSMom[0]->Fill(T2Par->fSMom,Weight);
  Hist->fSMom[1]->Fill(T2Par->fSMom,Weight);
  Hist->fQ->Fill   (T2Par->fQ   ,Weight);
}

//-----------------------------------------------------------------------------
void    TDegraderRpcAnaModule::FillTimeClusterHistograms(TimeClusterHist_t* Hist  , 
                                                    TStnTimeCluster*   Tc    ,
                                                    double             Weight) {
  Hist->fNsh->Fill  (Tc->NHits()     ,Weight);
  Hist->fNch->Fill  (Tc->NComboHits(),Weight);
  Hist->fT0->Fill   (Tc->T0()        ,Weight);
  Hist->fT0Err->Fill(Tc->T0Err()     ,Weight);
}

//_____________________________________________________________________________
void TDegraderRpcAnaModule::FillHistograms() {

  //  double       cos_th (-2.); //,  cl_e(-1.);
  // int          disk_id(-1);
//-----------------------------------------------------------------------------
// event histograms
//-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0],&fEvtPar);

  if (fEvtPar.fNGoodParticles > 1) FillEventHistograms(fHist.fEvent[1],&fEvtPar);
//-----------------------------------------------------------------------------
// helix histograms
// HEL_0: all
// HEL_1: helices in events with reconstructed tracks
//-----------------------------------------------------------------------------
  int nhel = fHelixBlock->NHelices();
  for (int i=0; i<nhel; i++) {
    TStnHelix*  hel = fHelixBlock->Helix(i);
    HelixPar_t* hp  = fHelixPar+i;
    FillHelixHistograms(fHist.fHelix[0],hel,hp);
  }
//-----------------------------------------------------------------------------
// SIMP histograms
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// track histograms, fill them only for the downstream e- hypothesis
//-----------------------------------------------------------------------------
  // TStnTrack*   trk;
  // TrackPar_t*  tp;
//-----------------------------------------------------------------------------
// fill GENP histograms
// GEN_0: all particles
//-----------------------------------------------------------------------------
  TGenParticle* genp;
  for (int i=0; i<fNGenp; i++) {
    genp = fGenpBlock->Particle(i);
    FillGenpHistograms(fHist.fGenp[0],genp);
  }
//-----------------------------------------------------------------------------
// fill DRpc histograms
// DRPC_0: all events
// DRPC_1: events with NP = 2
// DRPC_2: events with NP = 2 at VD13 (tracker front)
//-----------------------------------------------------------------------------
  if (fQVD13 == 0) {
    FillDRpcHistograms(fHist.fDRpc[0]);
    if (fNHitsVD10 == 2) {
      FillDRpcHistograms(fHist.fDRpc[1]);
    }
    if (fNHitsVD13 == 2) {
      FillDRpcHistograms(fHist.fDRpc[2]);
    }
  }
//-----------------------------------------------------------------------------
// fill two-track histograms
//-----------------------------------------------------------------------------
  for (int i=0; i<fNT2Par; i++) {
    T2Par_t* t2par = &fT2Par[i];
    FillT2Histograms(fHist.fT2[0],t2par);
    if (t2par->fQ == 0) FillT2Histograms(fHist.fT2[1],t2par);
  }
}


//-----------------------------------------------------------------------------
// 2014-04-30: it looks that reading the straw hits takes a lot of time - 
//              turn off by default by commenting it out
//-----------------------------------------------------------------------------
int TDegraderRpcAnaModule::Event(int ientry) {

  //  TDiskCalorimeter::GeomData_t disk_geom;

  fTrackBlock->GetEntry(ientry);
  // fTrackSeedBlock->GetEntry(ientry);
  // fHelixBlock->GetEntry(ientry);
  // fTimeClusterBlock->GetEntry(ientry);
  //  fTrackStrawHitBlock->GetEntry(ientry);
  // fClusterBlock->GetEntry(ientry);
  // fStrawHitBlock->GetEntry(ientry);
  // fCalDataBlock->GetEntry(ientry);
  fGenpBlock->GetEntry(ientry);
  fSimpBlock->GetEntry(ientry);
  fVDetBlock->GetEntry(ientry);
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
  fEvtPar.fNGenp            = fGenpBlock->NParticles();
  fEvtPar.fNSimp            = fSimpBlock->NParticles();
  fEvtPar.fSimp             = nullptr;
  fEvtPar.fPartE            = -1.;
  fEvtPar.fPionSurvProb     = 1.;
  fEvtPar.fNHelices         = fHelixBlock->NHelices();
  fEvtPar.fNTimeClusters    = fTimeClusterBlock->NTimeClusters();
//-----------------------------------------------------------------------------
// luminosity weight
//-----------------------------------------------------------------------------
  fLumWt = GetHeaderBlock()->LumWeight();
//-----------------------------------------------------------------------------
// figure the MC truth
// pi+ --> e+ nu case : determine the event weight
//-----------------------------------------------------------------------------
  fCPath = -1;
  for (int i=0; i<fEvtPar.fNSimp; i++) {
    TSimParticle* simp = fSimpBlock->Particle(i);
    int pdg_code       = simp->PDGCode();
    int generator_id   = simp->GeneratorID();         // MC process code

    if ((pdg_code == fPDGCode) and (generator_id == fMCProcessCode)) {
      fEvtPar.fSimp  = simp;
      fEvtPar.fPartE = simp->StartMom()->Energy();
//-----------------------------------------------------------------------------
// also calculate the path
//-----------------------------------------------------------------------------
      double x  = simp->EndPos()->X()+3904.;
      double y  = simp->EndPos()->Y();
      double dr = 250-sqrt(x*x+y*y);
      fCPath    = dr*simp->StartMom()->P()/simp->StartMom()->Pt();
    }

    if ((abs(pdg_code) == 211) && (simp->GeneratorID() == 56)) {
//-----------------------------------------------------------------------------
// found the pion, survival probability
//-----------------------------------------------------------------------------
      fEvtPar.fPionSurvProb  = exp(-simp->fEndProperTime);
      // break;
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
//-----------------------------------------------------------------------------
// process virtual detectors - for fSimp need parameters at tracker entrance
// use the first fit
//-----------------------------------------------------------------------------
  fSimPar.fTFront   = NULL;
  fSimPar.fTMid     = NULL;
  fSimPar.fTBack    = NULL;

  fQVD09      = 0;
  fNHitsVD09  = 0;
  fSMomVD09   = 0;

  fQVD10      = 0;
  fNHitsVD10  = 0;
  fSMomVD10   = 0;

  fQVD13      = 0;
  fNHitsVD13  = 0;
  fSMomVD13   = 0;
  
  int nvdhits = fVDetBlock->NStepPoints();

  for (int i=0; i<nvdhits; i++) {
    TStepPointMC* vdhit = fVDetBlock->StepPointMC(i);
    
    int pdg_id = vdhit->PDGCode();
    int sim_id = vdhit->SimID();
    TSimParticle* simp = fSimpBlock->FindParticle(sim_id);
    float mom = vdhit->Mom()->Mag();
    
    if ((abs(pdg_id) == 11) and (simp->NStrawHits() >= 20) and (mom > 30)) {
      if (vdhit->VolumeID() ==  9) {
//-----------------------------------------------------------------------------
// count step pointMCs at VD09
//-----------------------------------------------------------------------------
        if (vdhit->Mom()->Pz() > 0) {
          fHitVD09[fNHitsVD09] = vdhit;
          fNHitsVD09          += 1;
          fSMomVD09           += mom;
          fQVD09              += pdg_id;
        }
      }
      else if (vdhit->VolumeID() == 10) {
//-----------------------------------------------------------------------------
// count step pointMCs at VD10
//-----------------------------------------------------------------------------
        if (vdhit->Mom()->Pz() > 0) {
          fHitVD10[fNHitsVD10] = vdhit;
          fNHitsVD10          += 1;
          fSMomVD10           += mom;
          fQVD10              += pdg_id;
        }
      }
      else if (vdhit->VolumeID() == 13) {
        if (vdhit->Mom()->Pz() > 0) {
//-----------------------------------------------------------------------------
// count step pointMCs in front of the tracker (VD13)
//-----------------------------------------------------------------------------
          fHitVD13[fNHitsVD13] = vdhit;
          fNHitsVD13          += 1;
          fSMomVD13           += mom;
          fQVD13              += pdg_id;
        }
      }
    }
  }

  if (fNHitsVD13 == 2) {
//-----------------------------------------------------------------------------
// calculate parameters of the two "circles" and the distance between them
// look at the MuHitDisplay code
//-----------------------------------------------------------------------------
  }
//-----------------------------------------------------------------------------
// init calorimeter
//-----------------------------------------------------------------------------
  fNClusters  = fClusterBlock->NClusters();
  fNCalHits   = fCalDataBlock->NHits();
//-----------------------------------------------------------------------------
// init calorimeter clusters - remember, the first one not necessarily is the 
// most energetic
//-----------------------------------------------------------------------------
  fNClusters = fClusterBlock->NClusters();
  if (fNClusters == 0) fCluster = 0;
  else                 fCluster = fClusterBlock->Cluster(0);
//-----------------------------------------------------------------------------
// loop over tracks and calculate needed parameters
//-----------------------------------------------------------------------------
  fNStrawHits          = GetHeaderBlock()->fNStrawHits;

  fEvtPar.fNTracks[0]  = fTrackBlock->NTracks();
  if (fEvtPar.fNTracks[0] == 0) fTrack = 0;
  else                          fTrack = fTrackBlock->Track(0);

  fEvtPar.fNGoodTracks[0] = 0;
  fEvtPar.fNGoodTracks[1] = 0;
  fNMatchedTracks         = 0;
//-----------------------------------------------------------------------------
// determine the number of CalPatRec tracks - this assumes that the list of 
// tracks has been created by MergePatRec
//-----------------------------------------------------------------------------
  int ntrk = fTrackBlock->NTracks();

  InitTrackPar(fTrackBlock,fClusterBlock,fTrackPar,&fSimPar);
//-----------------------------------------------------------------------------
// additional initializations - helices and time clusters
//-----------------------------------------------------------------------------
  for (int itrk=0; itrk<ntrk; itrk++) {
    TStnTrack*   trk = fTrackBlock->Track(itrk);
    TrackPar_t*  tp  = fTrackPar+itrk;

    tp->fHelix       = nullptr;
    tp->fTimeCluster = nullptr;

    int ih = trk->fHelixIndex;
    if (ih < fEvtPar.fNHelices) {
      tp->fHelix = fHelixBlock->Helix(ih);
    }
//-----------------------------------------------------------------------------
// currently, at least in the KK case, an TStnHelix doesn't have a helix index stored
// find a time cluster closest to teh track in time
//-----------------------------------------------------------------------------
    TStnTimeCluster* closest_tc(nullptr);
    float min_dt = 1.e6;
    for (int itc=0; itc<fEvtPar.fNTimeClusters; itc++) {
      TStnTimeCluster* tc = fTimeClusterBlock->TimeCluster(itc);
      float dt = trk->T0()-tc->T0();
      if (fabs(dt) < fabs(min_dt)) {
        closest_tc = tc;
        min_dt     = dt;
      }
    }

    tp->fTimeCluster = closest_tc;
    tp->fDtTrackTc   = min_dt;
  }

  fNT2Par                   = 0;
  for (int it1=0; it1<ntrk-1; it1++) {
    TStnTrack*   t1  = fTrackBlock->Track(it1);
    TrackPar_t*  tp1 = fTrackPar+it1;

    for (int it2=it1+1; it2<ntrk; it2++) {
      TStnTrack*   t2  = fTrackBlock->Track(it2);
      TrackPar_t*  tp2 = fTrackPar+it2;

      T2Par_t* t2par = &fT2Par[fNT2Par];
      t2par->fItrk1  = it1;
      t2par->fItrk2  = it2;
      t2par->fSMom   = tp1->fP+tp2->fP;
      t2par->fQ      = t1->Charge()+t2->Charge();
      fNT2Par++;
    }
  }


  fEventPassedSelections = 0;

  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TDegraderRpcAnaModule::Debug() {

  TStnTrack* trk;
  TrackPar_t* tp;
  int ntrk = fTrackBlock->NTracks();

  for (int itrk=0; itrk<ntrk; itrk++) {
    trk = fTrackBlock->Track(itrk);
    tp  = &fTrackPar[itrk];
//-----------------------------------------------------------------------------
// bit 0: all tracks
//-----------------------------------------------------------------------------
    if (GetDebugBit(0) == 1) {
      GetHeaderBlock()->Print(Form("bit_000: All p = %10.3f",trk->Momentum()->P()));
    }
//-----------------------------------------------------------------------------
// bit 1: IDWord =0 0 tracks
//-----------------------------------------------------------------------------
    if (GetDebugBit(1) == 1) {
      if (tp->fIDWord[fBestID] == 0) {
	GetHeaderBlock()->Print(Form("bit_001: IDWord=0 p = %10.3f",trk->Momentum()->P()));
      }
    }
  }
//-----------------------------------------------------------------------------
// bit 3,4: outliers
//-----------------------------------------------------------------------------
  if (GetDebugBit(3) == 1) {
    if ((fSMomVD10 > 130) or (fSMomVD13 > 130)) {
      GetHeaderBlock()->Print(Form("fSMomVD10: %f fSMomVD13: %f",fSMomVD10,fSMomVD13));
    }
  }

  if (GetDebugBit(4) == 1) {
    if ((fNHitsVD13 == 2) and (fQVD13 == 0)) {
      GetHeaderBlock()->Print(Form("fSMomVD10: %10.3f fSMomVD13: %10.3f fCPath: %10.3f",fSMomVD10,fSMomVD13,fCPath));
    }
  }

  if (GetDebugBit(5) == 1) {
    for (int i=0; i<fNT2Par; i++) {
      T2Par_t* t2par = &fT2Par[i];
      if ((t2par->fSMom >= 89) and (t2par->fSMom < 91)) {
        GetHeaderBlock()->Print(Form("t2par->fSMom, t2par->fQ: %10.3f %10.3f",t2par->fSMom,t2par->fQ));
      }
    }
  }
}

//_____________________________________________________________________________
int TDegraderRpcAnaModule::EndJob() {
  printf("----- end job: ---- %s\n",GetName());
  return 0;
}

//_____________________________________________________________________________
void TDegraderRpcAnaModule::Test001() {
}

}
