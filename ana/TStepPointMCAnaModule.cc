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

#include "ana/InitVirtualDetectors.hh"

ClassImp(TStepPointMCAnaModule)
//-----------------------------------------------------------------------------
TStepPointMCAnaModule::TStepPointMCAnaModule(const char* name, const char* title):
  TStnModule(name,title)
{
  // fPdgCode       = 11;
  // fGeneratorCode = 28;
  fSpmcBlockName = "SpmcBlock";
  fVDetBlockName = "VdetBlock";

  fPdgDb = TDatabasePDG::Instance();

  for (int i=0; i<5000; i++) fParticleCache[i] = 0;

  SetParticleCache(   11,fPdgDb->GetParticle(   11)); // e
  SetParticleCache(  -11,fPdgDb->GetParticle(  -11)); // e
  SetParticleCache(   13,fPdgDb->GetParticle(   13)); // mu
  SetParticleCache(  -13,fPdgDb->GetParticle(  -13)); // mu
  SetParticleCache(   22,fPdgDb->GetParticle(   22)); // photon
  SetParticleCache(  211,fPdgDb->GetParticle(  211)); // pi^-
  SetParticleCache( -211,fPdgDb->GetParticle( -211)); // pi^-
  SetParticleCache( 2212,fPdgDb->GetParticle( 2212)); // proton
  SetParticleCache(-2212,fPdgDb->GetParticle(-2212)); // pbar
}

//-----------------------------------------------------------------------------
TStepPointMCAnaModule::~TStepPointMCAnaModule() {
}


//-----------------------------------------------------------------------------
void TStepPointMCAnaModule::BookEventHistograms(HistBase_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];
  EventHist_t* hist = (EventHist_t*) Hist;

  HBook1F(hist->fEventNumber,"evtnum",Form("%s: Event Number",Folder), 1000, 0,  1.e4,Folder);
  HBook1F(hist->fRunNumber  ,"runnum",Form("%s: Run   Number",Folder), 1000, 0,  1.e6,Folder);
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
  HBook1F(hist->fEKin           ,"ekin"    ,Form("%s: kinetic energy"  ,Folder), 400,   0,   100,Folder);

  HBook2F(hist->fYVsX           ,"y_vs_x"     ,Form("%s: Y vs X"       ,Folder), 100, -250,  250, 100, -250, 250, Folder);
  HBook2F(hist->fYVsZ           ,"y_vs_z"     ,Form("%s: Y vs Z"       ,Folder), 500, -250,  250, 500, -250, 250, Folder);
}

//-----------------------------------------------------------------------------
void TStepPointMCAnaModule::BookSimpHistograms(HistBase_t* Hist, const char* Folder) {
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
void TStepPointMCAnaModule::BookVDetHistograms(HistBase_t* Hist, const char* Folder) {

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
  HBook1F(hist->fEKin    ,"ekin"    ,Form("%s: kinetic energy",Folder), 400,  0, 100,Folder);

}

//-----------------------------------------------------------------------------
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

  book_spmc_histset[301] = 1;           // mu- p < 50 MeV/c
  book_spmc_histset[302] = 1;           // mu- p > 50 MeV/c

  book_spmc_histset[ 900] = 1;		// all particles , VD=9
  book_spmc_histset[ 913] = 1;		// mu- in VD=9

  book_spmc_histset[9100] = 1;		// all particles , VD=91
  book_spmc_histset[9113] = 1;		// mu- in VD=91

  for (int i=0; i<kNStepPointMCHistSets; i++) {
    if (book_spmc_histset[i] != 0) {
      sprintf(folder_name,"spmc_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fStepPointMC[i] = new StepPointMCHist_t;
      BookStepPointMCHistograms(fHist.fStepPointMC[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book SimParticle histograms ??  none so far
//-----------------------------------------------------------------------------
  int book_simp_histset[kNSimpHistSets];
  for (int i=0; i<kNSimpHistSets; i++) book_simp_histset[i] = 0;

  book_simp_histset[  0] = 0;		// all stopped muons
  book_simp_histset[  1] = 0;		// stopped muons p < 50 MeV/c

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
// book VDet histograms - so far need only VDET=9 (in front of the stopping target)
// virtual detector ID's are defined in DataProducts/inc/VirtualDetectorId.hh
//-----------------------------------------------------------------------------
  int book_vdet_histset[kNVDetHistSets];
  for (int i=0; i<kNVDetHistSets; i++) book_vdet_histset[i] = 0;

  book_vdet_histset[  9] = 1;		// all particles, VDET= 9 , ST_In
  book_vdet_histset[ 13] = 1;		// all particles, VDET=13 , TT_FrontHollow
  book_vdet_histset[ 14] = 1;		// all particles, VDET=13 , TT_FrontHollow, r > 40

  book_vdet_histset[109] = 1;		// e-  , VDET=9
  book_vdet_histset[209] = 1;		// e+  , VDET=9

  book_vdet_histset[301] = 1;		// all mu- , VDET=1: Coll1_In
  book_vdet_histset[302] = 1;		// all mu- , VDET=2: Coll1_Out
  book_vdet_histset[303] = 1;		// all mu- , VDET=3: Coll31_In
  book_vdet_histset[304] = 1;		// all mu- , VDET=4: Coll31_Out
  book_vdet_histset[305] = 1;		// all mu- , VDET=5: Coll32_In 
  book_vdet_histset[306] = 1;		// all mu- , VDET=6: Coll32_Out
  book_vdet_histset[307] = 1;		// all mu- , VDET=7: Coll5_In
  book_vdet_histset[308] = 1;		// all mu- , VDET=8: Coll5_Out
  book_vdet_histset[309] = 1;		// all mu- , VDET=9: ST_In

  book_vdet_histset[311] = 1;		// p<50 MeV/c mu- , VDET=1: Coll1_In
  book_vdet_histset[312] = 1;		// p<50 MeV/c mu- , VDET=2: Coll1_Out
  book_vdet_histset[313] = 1;		// p<50 MeV/c mu- , VDET=3: Coll31_In
  book_vdet_histset[314] = 1;		// p<50 MeV/c mu- , VDET=4: Coll31_Out
  book_vdet_histset[315] = 1;		// p<50 MeV/c mu- , VDET=5: Coll32_In 
  book_vdet_histset[316] = 1;		// p<50 MeV/c mu- , VDET=6: Coll32_Out
  book_vdet_histset[317] = 1;		// p<50 MeV/c mu- , VDET=7: Coll5_In
  book_vdet_histset[318] = 1;		// p<50 MeV/c mu- , VDET=8: Coll5_Out
  book_vdet_histset[319] = 1;		// p<50 MeV/c mu- , VDET=9: ST_In

  book_vdet_histset[321] = 1;		// p>50 MeV/c mu- , VDET=1: Coll1_In
  book_vdet_histset[322] = 1;		// p>50 MeV/c mu- , VDET=2: Coll1_Out
  book_vdet_histset[323] = 1;		// p>50 MeV/c mu- , VDET=3: Coll31_In
  book_vdet_histset[324] = 1;		// p>50 MeV/c mu- , VDET=4: Coll31_Out
  book_vdet_histset[325] = 1;		// p>50 MeV/c mu- , VDET=5: Coll32_In 
  book_vdet_histset[326] = 1;		// p>50 MeV/c mu- , VDET=6: Coll32_Out
  book_vdet_histset[327] = 1;		// p>50 MeV/c mu- , VDET=7: Coll5_In
  book_vdet_histset[328] = 1;		// p>50 MeV/c mu- , VDET=8: Coll5_Out
  book_vdet_histset[329] = 1;		// p>50 MeV/c mu- , VDET=9: ST_In

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
void TStepPointMCAnaModule::FillStepPointMCHistograms(HistBase_t* Hist, TStepPointMC* Step, SpmcData_t* SpmcData) {

  StepPointMCHist_t* hist = (StepPointMCHist_t*) Hist;
  
  int pdg_code = Step->PDGCode();

  hist->fVolumeID->Fill(Step->VolumeID());
  hist->fGenIndex->Fill(Step->GenIndex());
  hist->fSimID   ->Fill(Step->SimID());
  hist->fPDGCode[0] ->Fill(pdg_code);
  hist->fPDGCode[1] ->Fill(pdg_code);
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

  
  double m(0);
  if (SpmcData->fParticle) {
    m = SpmcData->fParticle->Mass()*1e3;   // convert GeV -> MeV
  }

  double ekin = sqrt(p*p+m*m)-m;
  hist->fEKin->Fill(ekin);

  float x = Step->Pos()->X();
  float y = Step->Pos()->Y();
  float z = Step->Pos()->Z();

  // hack for stage 2:

  if (fabs(z-2929.) < 0.1) x = x+3904;

  hist->fYVsX->Fill(x,y);		// useful for stage 2
  hist->fYVsZ->Fill(z,y);		// useful for stage 1
}

//-----------------------------------------------------------------------------
void TStepPointMCAnaModule::FillSimpHistograms(HistBase_t* Hist, TSimParticle* Simp, SimpData_t* Sd) {

  SimpHist_t* hist = (SimpHist_t*) Hist;
  
  hist->fVolumeID->Fill(Simp->fEndVolumeIndex);
  hist->fGeneratorID->Fill(Simp->fGeneratorID);
  
  float xe = Simp->EndPos()->X()+3904.;
  float ye = Simp->EndPos()->Y();
  float ze = Simp->EndPos()->Z();
  float te = Simp->EndPos()->T();

  hist->fTime->Fill(te);

  // hist->fParentMom->Fill(fParent->StartMom()->P());
  // hist->fParentPDG->Fill(fParent->PDGCode());

  hist->fStartMom->Fill(Simp->StartMom()->P());
  hist->fYVsX->Fill(xe,ye);
  hist->fXEndVsZEnd->Fill(ze,xe);

  if (Simp->fEndVolumeIndex == 2480) hist->fYVsX_2480->Fill(xe,ye);
  if (Simp->fEndVolumeIndex == 2513) hist->fYVsX_2513->Fill(xe,ye);
}

//-----------------------------------------------------------------------------
void TStepPointMCAnaModule::FillVDetHistograms(HistBase_t* Hist, TStepPointMC* Step) {

  VDetHist_t* hist = (VDetHist_t*) Hist;

  int id = (Step->VolumeID());

  VDetData_t* vdd = fVDet+id;
  
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
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TStepPointMCAnaModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks 'SpmcBlock' or 'StepPointMCBlock' (old)
//-----------------------------------------------------------------------------
  RegisterDataBlock(fSpmcBlockName.Data(),"TStepPointMCBlock",&fStepPointMCBlock);
  RegisterDataBlock("SimpBlock"          ,"TSimpBlock"       ,&fSimpBlock       );
  RegisterDataBlock(fVDetBlockName.Data(),"TStepPointMCBlock",&fVDetBlock       );
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


//-----------------------------------------------------------------------------
void TStepPointMCAnaModule::FillHistograms() {

//-----------------------------------------------------------------------------
// event histograms
//
// EVT_0: all events
//-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0]);
//-----------------------------------------------------------------------------
// StepPointMC histograms
// for beamline studies, there is only one 
//-----------------------------------------------------------------------------
  TStepPointMC* spmc;
  SpmcData_t    spmc_data;

  int nsteps = fStepPointMCBlock->NStepPoints();
  for (int i=0; i<nsteps; i++) {
    spmc         = fStepPointMCBlock->StepPointMC(i);
    float p      = spmc->Mom()->Mag();
    float t      = spmc->Time();
    int pdg_code = spmc->PDGCode();
    int abs_pdg_code = abs(pdg_code);
//-----------------------------------------------------------------------------
// particles of interest are electrons, pions, muons and photons,
// speed up the mass extraction
//-----------------------------------------------------------------------------
    if (abs_pdg_code < 2500) spmc_data.fParticle = GetParticleCache(pdg_code);
    else                     spmc_data.fParticle = NULL;

    if (spmc_data.fParticle == NULL) {
      printf(">>> WARNING: no particle with PDF code=%i in ROOT particle DB\n",pdg_code);
    }

    FillStepPointMCHistograms(fHist.fStepPointMC[0],spmc,&spmc_data);

    if      (pdg_code ==   11) {
      FillStepPointMCHistograms(fHist.fStepPointMC[1],spmc,&spmc_data);
      if (p >  60) FillStepPointMCHistograms(fHist.fStepPointMC[101],spmc,&spmc_data);
      if (t > 400) FillStepPointMCHistograms(fHist.fStepPointMC[102],spmc,&spmc_data);
      if (p <   2)              FillStepPointMCHistograms(fHist.fStepPointMC[103],spmc,&spmc_data);
      if ((p >=  2) && (p < 3)) FillStepPointMCHistograms(fHist.fStepPointMC[104],spmc,&spmc_data);
      if (p >   3)              FillStepPointMCHistograms(fHist.fStepPointMC[105],spmc,&spmc_data);
    }
    else if (pdg_code ==  -11) {
      FillStepPointMCHistograms(fHist.fStepPointMC[2],spmc,&spmc_data);
      if (p <   2)              FillStepPointMCHistograms(fHist.fStepPointMC[203],spmc,&spmc_data);
      if ((p >=  2) && (p < 3)) FillStepPointMCHistograms(fHist.fStepPointMC[204],spmc,&spmc_data);
      if (p >   3)              FillStepPointMCHistograms(fHist.fStepPointMC[205],spmc,&spmc_data);
    }
    else if (pdg_code ==   13) {
//-----------------------------------------------------------------------------
// mu-
//-----------------------------------------------------------------------------
      FillStepPointMCHistograms(fHist.fStepPointMC[3],spmc,&spmc_data);
      if   (p <  50)            FillStepPointMCHistograms(fHist.fStepPointMC[301],spmc,&spmc_data);
      else                      FillStepPointMCHistograms(fHist.fStepPointMC[302],spmc,&spmc_data);
    }
    else if (pdg_code ==  -13) FillStepPointMCHistograms(fHist.fStepPointMC[4],spmc,&spmc_data);
    else if (pdg_code ==   22) FillStepPointMCHistograms(fHist.fStepPointMC[5],spmc,&spmc_data);
    else if (pdg_code == -211) FillStepPointMCHistograms(fHist.fStepPointMC[6],spmc,&spmc_data);
    else if (pdg_code ==  211) FillStepPointMCHistograms(fHist.fStepPointMC[7],spmc,&spmc_data);
    else if (abs_pdg_code == 2212) FillStepPointMCHistograms(fHist.fStepPointMC[8],spmc,&spmc_data);  // protons+pbars
    else                           FillStepPointMCHistograms(fHist.fStepPointMC[9],spmc,&spmc_data);  // everything else

    if (spmc->VolumeID() == 9) {
      FillStepPointMCHistograms(fHist.fStepPointMC[900],spmc,&spmc_data);
      if (pdg_code == 13) FillStepPointMCHistograms(fHist.fStepPointMC[913],spmc,&spmc_data);
    }

    if (spmc->VolumeID() == 91) {
      FillStepPointMCHistograms(fHist.fStepPointMC[9100],spmc,&spmc_data);
      if (pdg_code == 13) FillStepPointMCHistograms(fHist.fStepPointMC[9113],spmc,&spmc_data);
    }

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
  fSimpBlock->GetEntry(ientry);
  fVDetBlock->GetEntry(ientry);
//-----------------------------------------------------------------------------
// assume electron in the first particle, otherwise the logic will need to 
// be changed
// if there are several hits, use the first one
//-----------------------------------------------------------------------------
  fNVDetHits = fVDetBlock->NStepPoints();

  //  int np     = fSimpBlock->NParticles();
  fProton    = fSimpBlock->Particle(0);
//   fMuon      = fSimpBlock->Particle(np-1);

//   TSimParticle* parent = fMuon;
// //-----------------------------------------------------------------------------
// // loop over "steps" in SimpBlock
// //-----------------------------------------------------------------------------
//   int loc = np-1;
//   while (loc >= 0) {
//     int parent_id = parent->ParentID();
//     parent = fSimpBlock->FindParticle(parent_id);
// //-----------------------------------------------------------------------------
// // sometimes history tree includes a scattered proton, which produces a pion,
// // decaying into a muon. In this case call pion a 'grandparent'
// //-----------------------------------------------------------------------------
//     if (parent && (parent->PDGCode() != 13)) break;
//   }

//  fParent = parent;
//-----------------------------------------------------------------------------
// everything is precalculated, fill histograms
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

