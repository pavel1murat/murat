///////////////////////////////////////////////////////////////////////////////
// 2017-05-24: analyse STNTUPLE produced by running on the output of Stage3 simulation
//             (output stream selected by tgtStopFilter)
// stopped muon is the last TSimParticle, its momentum in the end should be equal to zero
// 2019-10-04: old versions had "VDetBlock" named "VdetBlock"
// 
// use of tmp:
//
// Tmp(0) : 
// Tmp(1) : 
// 
// use of debug bits: bits 0-2 are reserved
// 0  : all events
// 1  : passed events
// 2  : rejected events
// 
// 3  : events with P(muon gp) > 8000
// 4  : events with electron or muon gparent=NULL
// 5  : events with electron or muon grandparent PDG code > 2500
// 6  : events with fSimp->EndPos()->Time > 900
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

#include "murat/ana/InitVirtualDetectors.hh"
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
  fVDetBlockName = "VDetBlock";
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

  HBook1F(hist->fStartMom   ,"mom"        ,Form("%s: start Mom"     ,Folder), 200,     0, 1000,Folder);
  HBook2F(hist->fYVsX       ,"y_vs_x"     ,Form("%s: yend vs Xend " ,Folder), 250,  -250, 250, 250, -250, 250,Folder);
  HBook2F(hist->fXEndVsZEnd ,"xe_vs_ze"   ,Form("%s: xend vs zend " ,Folder), 250,  -5000, 20000, 100, -5000, 5000,Folder);
  HBook2F(hist->fYVsX_2480  ,"y_vs_x_2480",Form("%s: Y vs X [2480]" ,Folder), 250,  -250, 250, 250, -250, 250,Folder);
  HBook2F(hist->fYVsX_2513  ,"y_vs_x_2513",Form("%s: Y vs X [2513]" ,Folder), 250,  -250, 250, 250, -250, 250,Folder);
}


//-----------------------------------------------------------------------------
void TMuonStopAnaModule::BookVDetHistograms(HistBase_t* Hist, const char* Folder) {

  VDetHist_t* hist = (VDetHist_t*) Hist;

  HBook1F(hist->fIndex   ,"index"   ,Form("%s: VD index"      ,Folder),1000, 0, 1000,Folder);
  HBook1F(hist->fPDGCode ,"pdg_code",Form("%s: PDG code"      ,Folder),2000,-1000, 1000,Folder);
  HBook1F(hist->fGenCode ,"gen_code",Form("%s: generator code",Folder), 100, -10, 90,Folder);
  HBook1F(hist->fMomentum,"mom"     ,Form("%s: Momentum"      ,Folder), 200, 0, 200,Folder);
  HBook1F(hist->fTime    ,"time"    ,Form("%s: Hit Time  "    ,Folder), 200, 0,2000,Folder);
  HBook2F(hist->fYVsX    ,"y_vs_x"  ,Form("%s: Y vs X (all)"  ,Folder), 250, -250, 250, 250, -250, 250,Folder);
  HBook2F(hist->fYVsZ    ,"y_vs_z"  ,Form("%s: Y vs Z (all)"  ,Folder), 250, -250, 250, 250, -250, 250,Folder);
  HBook1F(hist->fPt      ,"pt"      ,Form("%s: Pt"            ,Folder), 200, 0, 200,Folder);
  HBook1F(hist->fPp      ,"pp"      ,Form("%s: P(parallel)"   ,Folder), 200, 0, 200,Folder);
  HBook1F(hist->fTanTh   ,"tan_th"  ,Form("%s: tan(pitch ang)",Folder), 500, -1, 4,Folder);
}


//-----------------------------------------------------------------------------
void TMuonStopAnaModule::BookEventHistograms(HistBase_t* Hist, const char* Folder) {
  EventHist_t* hist = (EventHist_t*) Hist;

  HBook1F(hist->fEventNumber ,"evtnum" ,Form("%s: Event Number"  ,Folder), 1000, 0, 1.e4,Folder);
  HBook1F(hist->fRunNumber   ,"runnum" ,Form("%s: Run   Number"  ,Folder), 1000, 0, 1.e6,Folder);
  HBook1F(hist->fNVDetHits   ,"nvdh"   ,Form("%s: N VDET hits"   ,Folder),  100, 0,  100,Folder);
  HBook1F(hist->fNVDetHits_9 ,"nvdh_9" ,Form("%s: N VDET= 9 hits",Folder),  100, 0,  100,Folder);
  HBook1F(hist->fNVDetHits_13,"nvdh_13",Form("%s: N VDET=13 hits",Folder),  100, 0,  100,Folder);
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

  book_simp_histset[  0] = 1;		// all stopped muons
  book_simp_histset[  1] = 1;		// stopped muons p < 50 MeV/c
  book_simp_histset[  2] = 1;		// muons stopped in TS5

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
// virtual detector ID's are defined in DataProducts/inc/VirtualDetectorId.hh
// book VDet histograms - so far need only VDET=9 (in front of the stopping target)
//-----------------------------------------------------------------------------
  int book_vdet_histset[kNVDetHistSets];
  for (int i=0; i<kNVDetHistSets; i++) book_vdet_histset[i] = 0;

  book_vdet_histset[  9] = 1;		// all particles, VDET= 9 , ST_In
  book_vdet_histset[ 13] = 1;		// all particles, VDET=13 , TT_FrontHollow
  book_vdet_histset[ 14] = 1;		// all particles, VDET=14 , TT_FrontHollow, r > 40

  book_vdet_histset[109] = 1;		// e-  , VDET=9
  book_vdet_histset[209] = 1;		// e+  , VDET=9

  book_vdet_histset[301] = 1;		// all mu- , VDET= 1: Coll1_In
  book_vdet_histset[302] = 1;		// all mu- , VDET= 2: Coll1_Out
  book_vdet_histset[303] = 1;		// all mu- , VDET= 3: Coll31_In
  book_vdet_histset[304] = 1;		// all mu- , VDET= 4: Coll31_Out
  book_vdet_histset[305] = 1;		// all mu- , VDET= 5: Coll32_In 
  book_vdet_histset[306] = 1;		// all mu- , VDET= 6: Coll32_Out
  book_vdet_histset[307] = 1;		// all mu- , VDET= 7: Coll5_In
  book_vdet_histset[308] = 1;		// all mu- , VDET= 8: Coll5_Out
  book_vdet_histset[309] = 1;		// all mu- , VDET= 9: ST_In
  book_vdet_histset[310] = 1;		// all mu- , VDET=10: ST_Out
  book_vdet_histset[398] = 1;		// all mu- , VDET=98: mid-section TSu
  book_vdet_histset[399] = 1;		// all mu- , VDET=99: mid-section TSd

  book_vdet_histset[401] = 1;		// mu+ , VDET= 1: Coll1_In
  book_vdet_histset[402] = 1;		// mu+ , VDET= 2: Coll1_Out
  book_vdet_histset[403] = 1;		// mu+ , VDET= 3: Coll31_In
  book_vdet_histset[404] = 1;		// mu+ , VDET= 4: Coll31_Out
  book_vdet_histset[405] = 1;		// mu+ , VDET= 5: Coll32_In 
  book_vdet_histset[406] = 1;		// mu+ , VDET= 6: Coll32_Out
  book_vdet_histset[407] = 1;		// mu+ , VDET= 7: Coll5_In
  book_vdet_histset[408] = 1;		// mu+ , VDET= 8: Coll5_Out
  book_vdet_histset[409] = 1;		// mu+ , VDET= 9: ST_In
  book_vdet_histset[410] = 1;		// mu+ , VDET=10: ST_Out
  book_vdet_histset[498] = 1;		// mu+ , VDET=98: mid-section TSu
  book_vdet_histset[499] = 1;		// mu+ , VDET=99: mid-section TSd

  book_vdet_histset[501] = 1;		// p<50 MeV/c mu- , VDET= 1: Coll1_In
  book_vdet_histset[502] = 1;		// p<50 MeV/c mu- , VDET= 2: Coll1_Out
  book_vdet_histset[503] = 1;		// p<50 MeV/c mu- , VDET= 3: Coll31_In
  book_vdet_histset[504] = 1;		// p<50 MeV/c mu- , VDET= 4: Coll31_Out
  book_vdet_histset[505] = 1;		// p<50 MeV/c mu- , VDET= 5: Coll32_In 
  book_vdet_histset[506] = 1;		// p<50 MeV/c mu- , VDET= 6: Coll32_Out
  book_vdet_histset[507] = 1;		// p<50 MeV/c mu- , VDET= 7: Coll5_In
  book_vdet_histset[508] = 1;		// p<50 MeV/c mu- , VDET= 8: Coll5_Out
  book_vdet_histset[509] = 1;		// p<50 MeV/c mu- , VDET= 9: ST_In
  book_vdet_histset[510] = 1;		// p<50 MeV/c mu- , VDET=10: ST_Out
  book_vdet_histset[598] = 1;		// p<50 MeV/c mu- , VDET=98: mid-section TSu
  book_vdet_histset[599] = 1;		// p<50 MeV/c mu- , VDET=99: mid-section TSd

  book_vdet_histset[601] = 1;		// p>50 MeV/c mu- , VDET= 1: Coll1_In
  book_vdet_histset[602] = 1;		// p>50 MeV/c mu- , VDET= 2: Coll1_Out
  book_vdet_histset[603] = 1;		// p>50 MeV/c mu- , VDET= 3: Coll31_In
  book_vdet_histset[604] = 1;		// p>50 MeV/c mu- , VDET= 4: Coll31_Out
  book_vdet_histset[605] = 1;		// p>50 MeV/c mu- , VDET= 5: Coll32_In 
  book_vdet_histset[606] = 1;		// p>50 MeV/c mu- , VDET= 6: Coll32_Out
  book_vdet_histset[607] = 1;		// p>50 MeV/c mu- , VDET= 7: Coll5_In
  book_vdet_histset[608] = 1;		// p>50 MeV/c mu- , VDET= 8: Coll5_Out
  book_vdet_histset[609] = 1;		// p>50 MeV/c mu- , VDET= 9: ST_In
  book_vdet_histset[610] = 1;		// p>50 MeV/c mu- , VDET=10: ST_Out
  book_vdet_histset[698] = 1;		// p<50 MeV/c mu- , VDET=98: mid-section TSu
  book_vdet_histset[699] = 1;		// p<50 MeV/c mu- , VDET=99: mid-section TSd

  book_vdet_histset[2001] = 1;		// pbar , VDET= 1: Coll1_In
  book_vdet_histset[2002] = 1;		// pbar , VDET= 2: Coll1_Out
  book_vdet_histset[2003] = 1;		// pbar , VDET= 3: Coll31_In
  book_vdet_histset[2004] = 1;		// pbar , VDET= 4: Coll31_Out
  book_vdet_histset[2005] = 1;		// pbar , VDET= 5: Coll32_In 
  book_vdet_histset[2006] = 1;		// pbar , VDET= 6: Coll32_Out
  book_vdet_histset[2007] = 1;		// pbar , VDET= 7: Coll5_In
  book_vdet_histset[2008] = 1;		// pbar , VDET= 8: Coll5_Out
  book_vdet_histset[2009] = 1;		// pbar , VDET= 9: ST_In
  book_vdet_histset[2010] = 1;		// pbar , VDET=10: ST_Out
  book_vdet_histset[2091] = 1;		// pbar , VDET=91: before pbar window
  book_vdet_histset[2092] = 1;		// pbar , VDET=92: after pbar window
  book_vdet_histset[2098] = 1;		// pbar , VDET=98: mid-section TSu
  book_vdet_histset[2099] = 1;		// pbar , VDET=99: mid-section TSd

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
void TMuonStopAnaModule::FillSimpHistograms(HistBase_t* Hist, TSimParticle* Simp, SimpData_t* Sd) {

  SimpHist_t* hist = (SimpHist_t*) Hist;
  
  hist->fVolumeID->Fill(Simp->fEndVolumeIndex);
  hist->fGeneratorID->Fill(Simp->fGeneratorID);
  
  float xe = Simp->EndPos()->X()+3904.;
  float ye = Simp->EndPos()->Y();
  float ze = Simp->EndPos()->Z();
  float te = Simp->EndPos()->T();

  hist->fTime->Fill(te);
  hist->fParentMom->Fill(fParent->StartMom()->P());
  hist->fParentPDG->Fill(fParent->PDGCode());

  hist->fStartMom->Fill(Simp->StartMom()->P());
  hist->fYVsX->Fill(xe,ye);
  hist->fXEndVsZEnd->Fill(ze,xe);

  if (Simp->fEndVolumeIndex == 2480) hist->fYVsX_2480->Fill(xe,ye);
  if (Simp->fEndVolumeIndex == 2513) hist->fYVsX_2513->Fill(xe,ye);
}

//-----------------------------------------------------------------------------
VDetData_t*  TMuonStopAnaModule::GetVDetData(int ID) {
  VDetData_t* vd (NULL);

  if      (ID < 100) vd = fVDet+ID;

  return vd;
}


//-----------------------------------------------------------------------------
void TMuonStopAnaModule::FillVDetHistograms(HistBase_t* Hist, TStepPointMC* Step) {

  VDetHist_t* hist = (VDetHist_t*) Hist;

  int id = (Step->VolumeID());

  VDetData_t* vdd = GetVDetData(id) ; // fVDet+id;
  
  hist->fIndex   ->Fill(id);
  hist->fPDGCode ->Fill(Step->PDGCode());
  hist->fGenCode ->Fill(Step->GenIndex());
  hist->fMomentum->Fill(Step->Mom()->Mag());
  hist->fTime    ->Fill(Step->Time());

  // calculate local X and local Z

  double phi  = vdd->fPhiXZ*M_PI/180;

  double dx   = Step->Pos()->X()-vdd->fX;
  double dz   = Step->Pos()->Z()-vdd->fZ;

  double xloc =  dx*cos(phi) - dz*sin(phi);
  double zloc =  dx*sin(phi) + dz*cos(phi);

  hist->fYVsX    ->Fill(xloc,Step->Pos()->Y());
  hist->fYVsZ    ->Fill(zloc,Step->Pos()->Y());

  float py = Step->Mom()->Py();

  float px, pp, pt;
  if (vdd->fIZLocal == 3) {
    px = Step->Mom()->Px();
    pp = Step->Mom()->Pz();
  }
  else {
    px = Step->Mom()->Pz();
    pp = -Step->Mom()->Px();
  }

  pt = sqrt(px*px+py*py);

  hist->fPt->Fill(pt);
  hist->fPp->Fill(pp);
  hist->fTanTh->Fill(pt/pp);
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
  hist->fNVDetHits   ->Fill(fNVDetHits   );
  hist->fNVDetHits_9 ->Fill(fNVDetHits_9 );
  hist->fNVDetHits_13->Fill(fNVDetHits_13);
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
  float pmu = fMuon->StartMom()->P();
  if (pmu < 50) FillSimpHistograms(fHist.fSimp[1],fMuon,sd);

//-----------------------------------------------------------------------------
// SIMP[2]: muons stopped in TS5
//-----------------------------------------------------------------------------
  float ze = fMuon->EndPos()->Z();
  if (ze < 5000.) FillSimpHistograms(fHist.fSimp[2],fMuon,sd);
//-----------------------------------------------------------------------------
// VDET histograms
//-----------------------------------------------------------------------------
  for (int i=0; i<fNVDetHits; i++) {
    TStepPointMC* step = fVDetBlock->StepPointMC(i);

    if (step->VolumeID() ==  9) FillVDetHistograms(fHist.fVDet[ 9],step);
    if (step->VolumeID() == 13) {
      FillVDetHistograms(fHist.fVDet[13],step);
      float x = step->Pos()->X()+3904.;
      float y = step->Pos()->Y();
      float r = sqrt(x*x+y*y);
      if ((r >= 400) && (r < 800)) { 
	FillVDetHistograms(fHist.fVDet[14],step);
      }
      
    }

    if (step->PDGCode() == 11) {
      if (step->VolumeID() == 9) FillVDetHistograms(fHist.fVDet[109],step);
    }
    if (step->PDGCode() == -11) {
      if (step->VolumeID() == 9) FillVDetHistograms(fHist.fVDet[209],step);
    }
    if (step->PDGCode() == 13) {
      if (step->VolumeID() ==  1) FillVDetHistograms(fHist.fVDet[301],step);
      if (step->VolumeID() ==  2) FillVDetHistograms(fHist.fVDet[302],step);
      if (step->VolumeID() ==  3) FillVDetHistograms(fHist.fVDet[303],step);
      if (step->VolumeID() ==  4) FillVDetHistograms(fHist.fVDet[304],step);
      if (step->VolumeID() ==  5) FillVDetHistograms(fHist.fVDet[305],step);
      if (step->VolumeID() ==  6) FillVDetHistograms(fHist.fVDet[306],step);
      if (step->VolumeID() ==  7) FillVDetHistograms(fHist.fVDet[307],step);
      if (step->VolumeID() ==  8) FillVDetHistograms(fHist.fVDet[308],step);
      if (step->VolumeID() ==  9) FillVDetHistograms(fHist.fVDet[309],step);
      if (step->VolumeID() == 10) FillVDetHistograms(fHist.fVDet[310],step);
      if (step->VolumeID() == 98) FillVDetHistograms(fHist.fVDet[398],step);
      if (step->VolumeID() == 99) FillVDetHistograms(fHist.fVDet[399],step);

      float pmu = step->Mom()->Mag();
      if (pmu < 50) {
	if (step->VolumeID() ==  1) FillVDetHistograms(fHist.fVDet[501],step);
	if (step->VolumeID() ==  2) FillVDetHistograms(fHist.fVDet[502],step);
	if (step->VolumeID() ==  3) FillVDetHistograms(fHist.fVDet[503],step);
	if (step->VolumeID() ==  4) FillVDetHistograms(fHist.fVDet[504],step);
	if (step->VolumeID() ==  5) FillVDetHistograms(fHist.fVDet[505],step);
	if (step->VolumeID() ==  6) FillVDetHistograms(fHist.fVDet[506],step);
	if (step->VolumeID() ==  7) FillVDetHistograms(fHist.fVDet[507],step);
	if (step->VolumeID() ==  8) FillVDetHistograms(fHist.fVDet[508],step);
	if (step->VolumeID() ==  9) FillVDetHistograms(fHist.fVDet[509],step);
	if (step->VolumeID() == 10) FillVDetHistograms(fHist.fVDet[410],step);
	if (step->VolumeID() == 98) FillVDetHistograms(fHist.fVDet[598],step);
	if (step->VolumeID() == 99) FillVDetHistograms(fHist.fVDet[599],step);
      }
      else {
	if (step->VolumeID() ==  1) FillVDetHistograms(fHist.fVDet[601],step);
	if (step->VolumeID() ==  2) FillVDetHistograms(fHist.fVDet[602],step);
	if (step->VolumeID() ==  3) FillVDetHistograms(fHist.fVDet[603],step);
	if (step->VolumeID() ==  4) FillVDetHistograms(fHist.fVDet[604],step);
	if (step->VolumeID() ==  5) FillVDetHistograms(fHist.fVDet[605],step);
	if (step->VolumeID() ==  6) FillVDetHistograms(fHist.fVDet[606],step);
	if (step->VolumeID() ==  7) FillVDetHistograms(fHist.fVDet[607],step);
	if (step->VolumeID() ==  8) FillVDetHistograms(fHist.fVDet[608],step);
	if (step->VolumeID() ==  9) FillVDetHistograms(fHist.fVDet[609],step);
	if (step->VolumeID() == 10) FillVDetHistograms(fHist.fVDet[610],step);
	if (step->VolumeID() == 98) FillVDetHistograms(fHist.fVDet[698],step);
	if (step->VolumeID() == 99) FillVDetHistograms(fHist.fVDet[699],step);
      }
    }
    if (step->PDGCode() == -13) {
      if (step->VolumeID() ==  1) FillVDetHistograms(fHist.fVDet[401],step);
      if (step->VolumeID() ==  2) FillVDetHistograms(fHist.fVDet[402],step);
      if (step->VolumeID() ==  3) FillVDetHistograms(fHist.fVDet[403],step);
      if (step->VolumeID() ==  4) FillVDetHistograms(fHist.fVDet[404],step);
      if (step->VolumeID() ==  5) FillVDetHistograms(fHist.fVDet[405],step);
      if (step->VolumeID() ==  6) FillVDetHistograms(fHist.fVDet[406],step);
      if (step->VolumeID() ==  7) FillVDetHistograms(fHist.fVDet[407],step);
      if (step->VolumeID() ==  8) FillVDetHistograms(fHist.fVDet[408],step);
      if (step->VolumeID() ==  9) FillVDetHistograms(fHist.fVDet[409],step);
      if (step->VolumeID() == 10) FillVDetHistograms(fHist.fVDet[410],step);
      if (step->VolumeID() == 98) FillVDetHistograms(fHist.fVDet[498],step);
      if (step->VolumeID() == 99) FillVDetHistograms(fHist.fVDet[499],step);
    }
    if (step->PDGCode() == -2212) {
      if (step->VolumeID() ==  1) FillVDetHistograms(fHist.fVDet[2001],step);
      if (step->VolumeID() ==  2) FillVDetHistograms(fHist.fVDet[2002],step);
      if (step->VolumeID() ==  3) FillVDetHistograms(fHist.fVDet[2003],step);
      if (step->VolumeID() ==  4) FillVDetHistograms(fHist.fVDet[2004],step);
      if (step->VolumeID() ==  5) FillVDetHistograms(fHist.fVDet[2005],step);
      if (step->VolumeID() ==  6) FillVDetHistograms(fHist.fVDet[2006],step);
      if (step->VolumeID() ==  7) FillVDetHistograms(fHist.fVDet[2007],step);
      if (step->VolumeID() ==  8) FillVDetHistograms(fHist.fVDet[2008],step);
      if (step->VolumeID() ==  9) FillVDetHistograms(fHist.fVDet[2009],step);
      if (step->VolumeID() == 10) FillVDetHistograms(fHist.fVDet[2010],step);
      if (step->VolumeID() == 91) FillVDetHistograms(fHist.fVDet[2091],step);
      if (step->VolumeID() == 92) FillVDetHistograms(fHist.fVDet[2092],step);
      if (step->VolumeID() == 98) FillVDetHistograms(fHist.fVDet[2098],step);
      if (step->VolumeID() == 99) FillVDetHistograms(fHist.fVDet[2099],step);
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
  RegisterDataBlock("SimpBlock"          ,"TSimpBlock"       ,&fSimpBlock);
  RegisterDataBlock(fVDetBlockName.Data(),"TStepPointMCBlock",&fVDetBlock);
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();
//-----------------------------------------------------------------------------
// initialize virtual detector offsets - a convenience for histogram filling
//-----------------------------------------------------------------------------
  InitVirtualDetectors(fVDet,&fNVDet);

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
  fVDetBlock->GetEntry(ientry);

  fNVDetHits    = fVDetBlock->NStepPoints();
  fNVDetHits_9  = 0;
  fNVDetHits_13 = 0;
  fETot_13      = 0;

  for (int i=0; i<fNVDetHits; i++) {
    TStepPointMC* step = fVDetBlock->StepPointMC(i);
    if (step->VolumeID() ==  9) fNVDetHits_9 += 1;
    if (step->VolumeID() == 13) {
      float x = step->Pos()->X()+3904.;
      float y = step->Pos()->Y();
      float r = sqrt(x*x+y*y);
      if ((r >= 400) && (r < 800)) {
	fNVDetHits_13 += 1;
	fETot_13      += step->Mom()->Mag();
      }
    }
  }
//-----------------------------------------------------------------------------
// stopped muon - the last simulated particle in the list, proton - the first one
//-----------------------------------------------------------------------------
  int np  = fSimpBlock->NParticles();
  fProton = fSimpBlock->Particle(0);
  fMuon   = fSimpBlock->Particle(np-1);

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

  if (GetDebugBit(6) == 1) {
    if (fMuon->EndPos()->T() > 900) {
      GetHeaderBlock()->Print(Form("muon time = %12.5f",fMuon->EndPos()->T()));
    }
  }
}

//_____________________________________________________________________________
int TMuonStopAnaModule::EndJob() {
  printf("----- end job: ---- %s\n",GetName());
  return 0;
}

//_____________________________________________________________________________
void TMuonStopAnaModule::Test001() {
}

