//////////////////////////////////////////////////////////////////////////////
// use of tmp:
//
// use of debug bits: bits 0-2 are reserved
//  0  : all events
//  1  : passed events
//  2  : rejected events
//  3  : N(events with 4+ CRVC around a track)
//  4  : events with a track [ nhits>=15] and a CRV stub candidate
// 
///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <format>

#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TEnv.h"
#include "TSystem.h"
#include "TGLViewer.h"

#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/alg/TStntuple.hh"
#include "Stntuple/val/stntuple_val_functions.hh"
//------------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
#include "Stntuple/geom/TCrvNumerology.hh"

#include "murat/gui/TEvdCrvNode.hh"
#include "murat/gui/TEvdCrvView.hh"

#include "murat/gui/TEvdCaloView.hh"
#include "murat/gui/TEvdCaloNode.hh"

#include "murat/gui/TEvdDetNode.hh"
#include "murat/gui/TEvdManager.hh"

#include "murat/ana/TCrvTcAnaModule.hh"


using namespace murat;

ClassImp(murat::TCrvTcAnaModule)

namespace murat {
//-----------------------------------------------------------------------------
TCrvTcAnaModule::TCrvTcAnaModule(const char* name, const char* title):
  TStnModule(name,title)
{
  TH1::AddDirectory(0);
//-----------------------------------------------------------------------------
// detector numerology
//-----------------------------------------------------------------------------
  fTpm            = TrkPanelMap::Instance();
  fCaloChannelMap = TCaloChannelMap::Instance();
  fCrvChannelMap  = TCrvChannelMap::Instance();
  
  fFrRef          = nullptr;
  fHist           = new Hist_t;
  
  // for now, no division by the disk...
  // fCrystals.resize(kNCrystals);
  // for (int i=0; i<kNCrystals; i++) {
  //   crystal_t* cr = &fCrystals[i];
  //   cr->fCid      = i;
  // }
  
  // fHitCrystals.clear();

  // fListOfTrkParam.reserve(100);
  // fListOfTrkParam.resize (100);

  // fListOfCrvbParam.reserve(100);
  // fListOfCrvbParam.resize (100);

  auto vm = TEvdManager::Instance();
  vm->InitGeometry(); // if not initialized
  
  fEvdCrv = (TEvdCrv*) vm->FindSubdetector("CRV");
  fCanvasWf = nullptr;
  //   fEvdCrv->InitGeometry("Offline/Mu2eG4/geom/geom_common_extracted.txt");
}

//-----------------------------------------------------------------------------
TCrvTcAnaModule::~TCrvTcAnaModule() {
}

//-----------------------------------------------------------------------------
int TCrvTcAnaModule::BookCrvcHistograms(CrvcHist_t* Hist, CrvIndex_t* Index, TFolder* Folder) {

  std::string prefix = std::format("run:{:06d} sel:{:02d}",fRunNumber,Index->sel);
  std::string name, title;

  name  = "crvc_tc_dt";
  title = std::format("{} : dt(CRVC-TC)",prefix);
  fBookHist->HBook1F(Hist->h_crvc_tc_dt,name.data(),title.data(),1000,-1000,1000,Folder);   // in us...

  return 0;
}

//-----------------------------------------------------------------------------
int TCrvTcAnaModule::BookCrvpHistograms(CrvpHist_t* Hist, CrvIndex_t* Index, TFolder* Folder) {

  std::string prefix = std::format("run:{:06d} sel:{} roc:{:02d} feb:{}",
                                   fRunNumber,Index->sel, Index->roc, Index->feb);
  std::string name, title;

  name  = "ph";
  title = std::format("{} : ph",prefix);
  fBookHist->HBook1F(Hist->h_ph,name.data(),title.data(),100,0,1000,Folder);   // in us...

  name  = "npes";
  title = std::format("{} : npes",prefix);
  fBookHist->HBook1F(Hist->h_npes,name.data(),title.data(),100,0,500,Folder);   // in us...

  name  = "time";
  title = std::format("{} : time",prefix);
  fBookHist->HBook1F(Hist->h_time,name.data(),title.data(),100,0,1.e5,Folder);   // in us...

  name  = "crvp_tc_dt";
  title = std::format("{} : dt(CRVP-TC)",prefix);
  fBookHist->HBook1F(Hist->h_crvp_tc_dt,name.data(),title.data(),1000,-1000,1000,Folder);   // in us...

  name  = "feb";
  title = std::format("{} : feb",prefix);
  fBookHist->HBook1F(Hist->h_feb,name.data(),title.data(),100,0,100,Folder);   // in us...

  name  = "ch";
  title = std::format("{} : ch",prefix);
  fBookHist->HBook1F(Hist->h_ch,name.data(),title.data(),2000,0,2000,Folder);   // in us...

  name  = "crvp_tc_dt_vs_feb";
  title = std::format("{} : dt(CRVP-TC) vs feb",prefix);
  fBookHist->HBook2F(Hist->h_crvp_tc_dt_vs_feb,name,title,100,0,100,1000,-1000,1000,Folder);   // in us...

  return 0;
}

//-----------------------------------------------------------------------------
int TCrvTcAnaModule::BookCrvsHistograms(CrvsHist_t* Hist, CrvIndex_t* Index, TFolder* Folder) {

  std::string prefix = std::format("run:{:06d} sel:{:02d}",fRunNumber,Index->sel);
  std::string name, title;

  name  = "crvc_tc_dt";
  title = std::format("{} : dt(CRVC-TC)",prefix);
  fBookHist->HBook1F(Hist->h_crvc_tc_dt,name,title,1000,-1000,1000,Folder);   // in us...

  name  = "crvp_tc_dt";
  title = std::format("{} : dt(CRVP-TC)",prefix);
  fBookHist->HBook1F(Hist->h_crvp_tc_dt,name,title,1000,-1000,1000,Folder);   // in us...

  return 0;
}

//-----------------------------------------------------------------------------
  int TCrvTcAnaModule::BookEventHistograms(EventHist_t*  Hist, EvtIndex_t*  Index, TFolder* Folder) {
    
    std::string prefix = std::format("run:{:06d} sel:{}",fRunNumber,Index->sel);
    std::string name, title;

    name  = "evn";
    title = std::format("{} : event number",prefix);
    fBookHist->HBook1F(Hist->h_evn,name,title,1000,0,1.e7,Folder);   // in us...

    return 0;
  }


//-----------------------------------------------------------------------------
int TCrvTcAnaModule::BookFebHistograms(FebHist_t* Hist, CrvIndex_t* Index, TFolder* Folder) {

  std::string prefix = std::format("run:{:06d} roc:{} feb:{:02d}",fRunNumber,Index->roc, Index->feb);
  std::string name, title;

  name  = "sbid";
  title = std::format("{} : SBID",prefix);
  fBookHist->HBook1F(Hist->h_sbid,name,title,1500,0,1500,Folder);

  name  = "ch_vs_crvp_tc_dt";
  title = std::format("{} : channel vs [T(pulse)-T(trk TC)]",prefix);
  fBookHist->HBook2F(Hist->h_ch_vs_crvp_tc_dt,name,title,400,-1000,1000,64,0,64,Folder);

//-----------------------------------------------------------------------------
// CRV reco pulses
//-----------------------------------------------------------------------------
  std::string folder_name = std::format("crvp");
  TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
  if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
  Hist->crvp = new CrvpHist_t();
  BookCrvpHistograms(Hist->crvp,Index,fol);
  
  return 0;
}

//-----------------------------------------------------------------------------
// by FEB: 32 scintillation bars per FEB 
//-----------------------------------------------------------------------------
int TCrvTcAnaModule::BookRocHistograms(RocHist_t* Hist, CrvIndex_t* Index, TFolder* Folder) {

  std::string prefix = std::format("run:{:06d} roc:{:02d}",fRunNumber,Index->roc);
  std::string name, title;

  name  = "sbid";
  title = std::format("{} : SBID",prefix);
  fBookHist->HBook1F(Hist->h_sbid,name,title,1500,0,1500,Folder);

  name  = "feb_vs_crvp_tc_dt";
  title = std::format("{} : FEB vs DT = T(CRVP)-T(TC)",prefix);
  fBookHist->HBook2F(Hist->h_feb_vs_crvp_tc_dt,name,title,1000,-1000,1000,30,0,30,Folder);

//-----------------------------------------------------------------------------
// CRV reco pulses
//-----------------------------------------------------------------------------
  std::string folder_name = std::format("crvp");
  TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
  if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
  Hist->crvp = new CrvpHist_t();
  BookCrvpHistograms(Hist->crvp,Index,fol);

//-----------------------------------------------------------------------------
// by FEB - 24 or 25 FEBs per ROC ?
//-----------------------------------------------------------------------------
  int n_feb_histsets(30);

  for (int i=0; i<n_feb_histsets; i++) {
    std::string folder_name = std::format("feb_{:02d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    Hist->feb[i] = new FebHist_t();
    Index->feb = i;
    BookFebHistograms(Hist->feb[i],Index,fol);
  }

  return 0;
}

//-----------------------------------------------------------------------------
int TCrvTcAnaModule::BookHistograms(Hist_t* Hist, TFolder* Folder) {

  // std::string prefix = std::format("");
  // std::string name, title;

  CrvIndex_t index;

  std::string prefix = std::format("run:{:06d}",fRunNumber);
  std::string name, title;

  name  = "sbid";
  title = std::format("{} : SBID",prefix);
  fBookHist->HBook1F(Hist->h_sbid,name,title,1500,0,1500,Folder);

  name  = "feb_vs_chid";
  title = std::format("{} : FEB vs CHID",prefix);
  fBookHist->HBook2F(Hist->h_feb_vs_chid,name,title,64,0,64,30,0,30,Folder);

  name  = "crvp_tc_dt_vs_chid";
  title = std::format("{} : dt(CRVP-TC) vs CHID",prefix);
  fBookHist->HBook2F(Hist->h_crvp_tc_dt_vs_chid,name,title,2400,0,2400,200,-500,500,Folder);

  name  = "feb_vs_sbid_0";
  title = std::format("{} : dt vs SBID ROC=1",prefix);
  fBookHist->HBook2F(Hist->h_feb_vs_sbid[0],name,title,600,0,600,30,0,30,Folder);

  name  = "feb_vs_sbid_1";
  title = std::format("{} : dt vs SBID ROC=2",prefix);
  fBookHist->HBook2F(Hist->h_feb_vs_sbid[1],name,title,600,0,600,30,0,30,Folder);

  name  = "chid";
  title = std::format("{} : offline channel ID",prefix);
  fBookHist->HBook1F(Hist->h_chid,name,title,2400,0,2400,Folder);

//-----------------------------------------------------------------------------
// CRV coincidence clusters
//-----------------------------------------------------------------------------
  int book_crvc_histset[kNCrvcHistSets];

  for (int i=0; i<kNCrvcHistSets; i++) { book_crvc_histset[i] = 0; }

  book_crvc_histset[0] = 1;

  for (int i=0; i<kNCrvcHistSets; i++) {
    if (book_crvc_histset[i] == 0) continue;
    std::string folder_name = std::format("crvc_{:02d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    Hist->crvc[i] = new CrvcHist_t();
    index.sel = i;
    BookCrvcHistograms(Hist->crvc[i],&index,fol);
  }

//-----------------------------------------------------------------------------
// CRVP histograms
//-----------------------------------------------------------------------------
  int book_crvp_histset[10];

  for (int i=0; i<kNCrvpHistSets; i++) { book_crvp_histset[i] = 0; }

  book_crvp_histset[0] = 1;             // all
  book_crvp_histset[1] = 1;             // ntc=1, dt_530 < 30
  book_crvp_histset[2] = 1;             // ntc=1, dt_590 < 30

  for (int i=0; i<kNCrvpHistSets; i++) {
    if (book_crvp_histset[i] == 0) continue;
    std::string folder_name = std::format("crvp_{:02d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    Hist->crvp[i] = new CrvpHist_t();
    index.sel = i;
    BookCrvpHistograms(Hist->crvp[i],&index,fol);
  }

//-----------------------------------------------------------------------------
// CRVS histograms - by sector
//-----------------------------------------------------------------------------
  int book_crvs_histset[10];

  for (int i=0; i<kNCrvsHistSets; i++) { book_crvs_histset[i] = 0; }

  book_crvs_histset[0] = 1;             // EX
  book_crvs_histset[1] = 1;             // T1
  book_crvs_histset[2] = 1;             // T2
  book_crvs_histset[3] = 1;             // all M's

  for (int i=0; i<kNCrvsHistSets; i++) {
    if (book_crvs_histset[i] == 0) continue;
    std::string folder_name = std::format("crvs_{:02d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    Hist->crvs[i] = new CrvsHist_t();
    index.sel = i;
    BookCrvsHistograms(Hist->crvs[i],&index,fol);
  }

//-----------------------------------------------------------------------------
// by ROC, links 0 and 3 --> rocs #1 and #4
//-----------------------------------------------------------------------------
  int book_roc_histset[kNRocHistSets];

  for (int i=0; i<kNRocHistSets; i++) { book_roc_histset[i] = 0; }

  book_roc_histset[1] = 1;
  book_roc_histset[2] = 1;
  book_roc_histset[4] = 1;

  for (int i=0; i<kNRocHistSets; i++) {
    if (book_roc_histset[i] == 0) continue;
    std::string folder_name = std::format("roc_{:02d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    Hist->roc[i] = new RocHist_t();
    index.roc = i;
    BookRocHistograms(Hist->roc[i],&index,fol);
  }

  return 0;
}

//-----------------------------------------------------------------------------
// need to optimize the filling time
//-----------------------------------------------------------------------------
  int TCrvTcAnaModule::FillCrvcHistograms(CrvcHist_t* Hist, TCrvCoincidenceCluster* Crvc, crvc_param_t* Crvcp) {
  // filling histograms: plot time differences between
  Hist->h_crvc_tc_dt->Fill(Crvcp->dt_min);

  return 0;
}

//-----------------------------------------------------------------------------
// need to optimize the filling time
//-----------------------------------------------------------------------------
  int TCrvTcAnaModule::FillCrvpHistograms(CrvpHist_t* Hist, TCrvRecoPulse* Crvp, crvp_param_t* Crvpp) {
  // filling histograms: plot time differences between
  Hist->h_ph->Fill(Crvp->fPed);
  Hist->h_npes->Fill(Crvp->fPes);
  Hist->h_time->Fill(Crvp->fTime);
  Hist->h_feb->Fill(Crvp->fFeb);
  Hist->h_ch->Fill(Crvp->Chid());

  Hist->h_crvp_tc_dt->Fill(Crvpp->dt_min);
  Hist->h_crvp_tc_dt_vs_feb->Fill(Crvp->fFeb,Crvpp->dt_min);

  return 0;
}

//-----------------------------------------------------------------------------
// need to optimize the filling time
//-----------------------------------------------------------------------------
  int TCrvTcAnaModule::FillCrvsHistograms() {

    for (int i=0; i<fNCrvp; i++) {
      //TCrvRecoPulse* crvp = fCrvpBlock->Pulse(i);
      crvp_param_t* crvpp = &fListOfCrvpParam[i];
      fHist->crvs[crvpp->sector]->h_crvp_tc_dt->Fill(crvpp->dt_min);
    }

    for (int i=0; i<fNCrvc; i++) {
      TCrvCoincidenceCluster* crvc = fCrvcBlock->Cluster(i);
      crvc_param_t*          crvcp = &fListOfCrvcParam[i];

      int sector = crvc->Sector();
                                        // all M's together
      if (sector >= 3) sector = 3;

      fHist->crvs[sector]->h_crvc_tc_dt->Fill(crvcp->dt_min);
    }

    return 0;
  }

//-----------------------------------------------------------------------------
  int TCrvTcAnaModule::FillFebHistograms(FebHist_t* Hist, TCrvRecoPulse* Crvp, crvp_param_t* Crvpp) {

    Hist->h_sbid->Fill(Crvp->fSbid);
    Hist->h_ch_vs_crvp_tc_dt->Fill(Crvpp->dt_min,Crvp->fFebCh);

    return 0;
  }

//-----------------------------------------------------------------------------
  int TCrvTcAnaModule::FillRocHistograms(RocHist_t* Hist, TCrvRecoPulse* Crvp, crvp_param_t* Crvpp) {

    Hist->h_sbid->Fill(Crvp->fSbid);
    Hist->h_feb_vs_crvp_tc_dt->Fill(Crvpp->dt_min,Crvp->fFeb);

    return 0;
  }

//-----------------------------------------------------------------------------
// need to optimize the filling time
//-----------------------------------------------------------------------------
int TCrvTcAnaModule::FillHistograms() {
  // filling histograms: plot time differences between

//-----------------------------------------------------------------------------
// missing parameters
//-----------------------------------------------------------------------------
  fListOfCrvcParam.resize(fNCrvc);

  for (int i=0; i<fNCrvc; i++) {
    TCrvCoincidenceCluster*  crvc  = fCrvcBlock->Cluster(i);
    crvc_param_t*   crvcp = &fListOfCrvcParam[i];
    crvcp->sector = crvc->Sector();
    crvcp->dt_min = 1.e12;
    for (int itc=0; itc<fNTc; itc++) {
      TStnTimeCluster* tc = fTcBlock->TimeCluster(itc);
      float dt    = crvc->fStartTime-tc->T0();
      if (fabs(dt) < fabs(crvcp->dt_min)) {
        crvcp->dt_min = dt;
      }
    }
  }

  fListOfCrvpParam.resize(fNCrvp);

  for (int i=0; i<fNCrvp; i++) {
    TCrvRecoPulse*  crvp  = fCrvpBlock->Pulse(i);
    crvp_param_t*   crvpp = &fListOfCrvpParam[i];

    TCrvChannelMap::Data_t* cdata = fCrvChannelMap->ch_data_by_offline(crvp->Chid());
    crvpp->sector = cdata->sector;
    crvpp->dt_min = 1.e12;
    for (int itc=0; itc<fNTc; itc++) {
      TStnTimeCluster* tc = fTcBlock->TimeCluster(itc);
      float dt    = crvp->fTime-tc->T0();
      if (fabs(dt) < fabs(crvpp->dt_min)) {
        crvpp->dt_min = dt;
      }
    }
  }

//-----------------------------------------------------------------------------
// CRVS histograms - by sector
//-----------------------------------------------------------------------------
  FillCrvsHistograms();

//-----------------------------------------------------------------------------
// CRV pulse-only histograms
//-----------------------------------------------------------------------------
  for (int i=0; i<fNCrvp; i++) {
    TCrvRecoPulse* crvp  = fCrvpBlock->Pulse(i);
    // crvp_param_t*  crvpp = &fListOfCrvpParam[i];

    // this is a global histogram
    fHist->h_feb_vs_chid->Fill(crvp->Chid(),crvp->fFeb);

    // have two histograms - one per roc , to color them
    // ROCs 1 and 2 --> hists 0 and 1
    fHist->h_feb_vs_sbid[crvp->fRoc-1]->Fill(crvp->fSbid,crvp->fFeb);
    
    fHist->h_sbid->Fill(crvp->fSbid);
    fHist->h_chid->Fill(crvp->Chid());
  }
//-----------------------------------------------------------------------------
// CRVC-TC histograms
//-----------------------------------------------------------------------------
  for (int i=0; i<fNCrvc; i++) {
    TCrvCoincidenceCluster*  crvc = fCrvcBlock->Cluster(i);
    crvc_param_t* crvcp           =&fListOfCrvcParam[i];
      
    FillCrvcHistograms(fHist->crvc[0],crvc,crvcp);
  }

//-----------------------------------------------------------------------------
// CRVP-TC histograms
//-----------------------------------------------------------------------------
  for (int i=0; i<fNCrvp; i++) {
    TCrvRecoPulse* crvp  = fCrvpBlock->Pulse(i);
    crvp_param_t*  crvpp = &fListOfCrvpParam[i];

    fHist->h_crvp_tc_dt_vs_chid->Fill(crvp->Chid(),crvpp->dt_min);

    RocHist_t* roc_hist = fHist->roc[crvp->fRoc];
    FillRocHistograms(roc_hist,crvp,crvpp);

    FebHist_t* feb_hist = fHist->roc[crvp->fRoc]->feb[crvp->fFeb];
    FillFebHistograms(feb_hist,crvp,crvpp);

    FillCrvpHistograms(fHist->crvp[0],crvp,crvpp);

    if (fNTc == 1) {
      if      (fabs(crvpp->dt_min - 80) < 30) {
//-----------------------------------------------------------------------------
// CRVP[1] : first peak
//-----------------------------------------------------------------------------
        FillCrvpHistograms(fHist->crvp[1],crvp,crvpp);
      }
      else if (fabs(crvpp->dt_min - 20) < 30) {
//-----------------------------------------------------------------------------
// CRVP[2] : second peak
//-----------------------------------------------------------------------------
        FillCrvpHistograms(fHist->crvp[2],crvp,crvpp);
      }
    }
  }
  
  return 0;
}


//-----------------------------------------------------------------------------
  int TCrvTcAnaModule::FillEventHistograms(EventHist_t* Hist) {
    return 0;
  }


//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TCrvTcAnaModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
  // RegisterDataBlock("CaloHitBlock"     ,"TCaloHitBlock"       ,&fCaloHitBlock     );
  // RegisterDataBlock("CaloRecoDigiBlock","TCaloRecoDigiBlock"  ,&fCaloRecoDigiBlock);
  // RegisterDataBlock("CaloClusterBlock" ,"TStnClusterBlock"    ,&fCaloClusterBlock );
  RegisterDataBlock("CrvcBlock"        ,"TCrvClusterBlock"    ,&fCrvcBlock        );
  RegisterDataBlock("CrvpBlock"        ,"TCrvPulseBlock"      ,&fCrvpBlock        );
  RegisterDataBlock("CrvdBlock"        ,"TCrvDigiBlock"       ,&fCrvdBlock        );
  RegisterDataBlock("TimeClusterBlock" ,"TStnTimeClusterBlock",&fTcBlock          );
  //  RegisterDataBlock("ComboHitBlock"    ,"TComboHitBlock"      ,&fChBlock          );
  //  RegisterDataBlock("TrackBlock"       ,"TStrTrackBlock"      ,&fTrackBlock       );
 
  return 0;
}


//_____________________________________________________________________________
int TCrvTcAnaModule::BeginRun() {
  fRunNumber = GetHeaderBlock()->RunNumber();
  TStntuple::Init(fRunNumber);
  fTpm->Init(fRunNumber);
  fCrvChannelMap->Init(fRunNumber);
//-----------------------------------------------------------------------------
// book histograms at begin run, there will be a flag telling whether we need
// to re-book them for every new run
// don't do that by default
//-----------------------------------------------------------------------------
  // relative positions of subdetectors wrt each other
  int print_level(0);
  fDetPos.GetFromRunDb(fRunNumber,print_level);
  
  BookHistograms(fHist,fFolder);
  return 0;
}


//-----------------------------------------------------------------------------
// 2026-09-09
//-----------------------------------------------------------------------------
int TCrvTcAnaModule::Event(int IEntry) {

  fTcBlock->GetEntry(IEntry);
  
  fCrvpBlock->GetEntry(IEntry);
  fCrvcBlock->GetEntry(IEntry);
  fCrvdBlock->GetEntry(IEntry);
  
  //  fTrackBlock->GetEntry(IEntry);
//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
  fNCrvc         = fCrvcBlock->NClusters();
  fNCrvd         = fCrvdBlock->NDigis();
  fNCrvp         = fCrvpBlock->NPulses();
  fNTc           = fTcBlock->NTimeClusters();

  //  CalculateMissingParameters();
  
  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TCrvTcAnaModule::Debug() {

  if (GetDebugBit(0) == 1) {
    GetHeaderBlock()->Print(Form("bit_00: all events"));
  }

  if (GetDebugBit(3) == 1) {
  }
  
}

//_____________________________________________________________________________
int TCrvTcAnaModule::EndJob() {
  return 0;
}


//-----------------------------------------------------------------------------
  int TCrvTcAnaModule::InitEvd() {
    int rc(0);

    auto vm = TEvdManager::Instance();
    // vm->InitGeometry(); // if not initialized
    
    vm->AddDataBlock("HeaderBlock"      ,GetHeaderBlock());
    // vm->AddDataBlock("CaloHitBlock"     ,fCaloHitBlock);
    // vm->AddDataBlock("CaloRecoDigiBlock",fCaloRecoDigiBlock);
    // vm->AddDataBlock("CaloClusterBlock" ,fCaloClusterBlock);
    
    vm->AddDataBlock("CrvpBlock"        ,fCrvpBlock);
    vm->AddDataBlock("CrvdBlock"        ,fCrvdBlock);
    vm->AddDataBlock("CrvcBlock"        ,fCrvcBlock);

//-----------------------------------------------------------------------------
// CRV-only view
//-----------------------------------------------------------------------------
    TEvdCrvNode* node = new TEvdCrvNode("crv_node");
    TEvdCrvView* view = new TEvdCrvView("crv_zy_view",TGLViewer::kCameraOrthoZOY);
    view->AddNode(node);
    vm->AddView(view);

    view = new TEvdCrvView("crv_xy_view",TGLViewer::kCameraOrthoXOY);
    view->AddNode(node);
    vm->AddView(view);
//-----------------------------------------------------------------------------
// 3D detector view
//-----------------------------------------------------------------------------
    // TEvdDetNode* det_node = new TEvdDetNode("det_node");
    // TEvdView*    det_view = new TEvdView   ("det_3D_view",TGLViewer::kCameraPerspXOY);
    // det_view->AddNode(det_node);
    // vm->AddView(det_view);

    return rc;
  }
  
//-----------------------------------------------------------------------------
  int TCrvTcAnaModule::DisplayEvent() {
    int rc(0);
    auto vm = TEvdManager::Instance();

    vm->DisplayEvent();

    return rc;
  }
  
//-----------------------------------------------------------------------------
void TCrvTcAnaModule::PlotWaveform(int ChannelID, bool NewCanvas) {

  if (fCanvasWf != nullptr) {
    if (NewCanvas) {
      // delete old canvas
      delete fCanvasWf;
      fCanvasWf = nullptr;
    }
  }

  if (fCanvasWf == nullptr) {
    fCanvasWf = new TCanvas("c","c",1200,800);
  }

  int nd = fCrvdBlock->NDigis();

  int nwf = 0;

  for (int i=0; i<nd; i++) {
    TCrvDigi* digi = fCrvdBlock->Digi(i);
    int chid = 4*digi->Sbid()+digi->Sipm();
    if (chid == ChannelID) {
      nwf++;
      int ns = digi->Ns();
      // assume that a single waveform is shorter than 1 us
      TH1F* hist = new TH1F(Form("h_wf_%04i_%i",ChannelID,i),
                            Form("event: %d:%d:%d waveform sipmid=_%04i",
                                 GetHeaderBlock()->RunNumber   (),
                                 GetHeaderBlock()->SubrunNumber(),
                                 GetHeaderBlock()->EventNumber (),
                                 ChannelID),20,0,20);
      hist->SetMinimum(2000);
      hist->SetMaximum(4200);

      for (int is=0; is<ns; is++) {
        int adc = digi->Adc()[is];
        hist->SetBinContent(is+1,adc);
        hist->SetBinError  (is+1,  0);
      }

      if (NewCanvas) {
        hist->Draw("hist");
      }
      else {
        hist->Draw("hist,sames");
      }
    }
  }

  fCanvasWf->Modified();
  fCanvasWf->Update();

  std::cout << std::format("plotted {} waveforms\n",nwf);
}
  

}
