//////////////////////////////////////////////////////////////////////////////
// use of tmp:
//
// use of debug bits: bits 0-2 are reserved
//  0  : all events
//  1  : passed events
//  2  : rejected events
//  3  : N(events with 4+ CRVC around a track)
//  4  : 
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

#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/alg/TStntuple.hh"
#include "Stntuple/geom/TDisk.hh"
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

#include "murat/ana/TCrvTimeAnaModule.hh"


using namespace murat;

ClassImp(murat::TCrvTimeAnaModule)

namespace murat {
//-----------------------------------------------------------------------------
TCrvTimeAnaModule::TCrvTimeAnaModule(const char* name, const char* title):
  TStnModule(name,title)
{
  TH1::AddDirectory(0);
//-----------------------------------------------------------------------------
// detector numerology
//-----------------------------------------------------------------------------
  fTpm            = TrkPanelMap::Instance();
  fCaloChannelMap = TCaloChannelMap::Instance();
  
  fFrRef          = nullptr;
  fHist           = new Hist_t;
  
  // for now, no division by the disk...
  fCrystals.resize(kNCrystals);
  for (int i=0; i<kNCrystals; i++) {
    crystal_t* cr = &fCrystals[i];
    cr->fCid      = i;
  }
  
  fHitCrystals.clear();

  fListOfTrkParam.reserve(100);
  fListOfTrkParam.resize (100);

  fListOfCrvbParam.reserve(100);
  fListOfCrvbParam.resize (100);

  auto vm = TEvdManager::Instance();
  vm->InitGeometry(); // if not initialized
  
  fEvdCrv = (TEvdCrv*) vm->FindSubdetector("CRV");
  fCanvasWf = nullptr;
  //   fEvdCrv->InitGeometry("Offline/Mu2eG4/geom/geom_common_extracted.txt");
}

//-----------------------------------------------------------------------------
TCrvTimeAnaModule::~TCrvTimeAnaModule() {
}

//-----------------------------------------------------------------------------
int TCrvTimeAnaModule::BookCrvsHistograms(CrvsHist_t* Hist, CrvIndex_t* Index, TFolder* Folder) {

  std::string prefix = std::format("run:{:06d} sel:{} sector:{}",fRunNumber,Index->sel,Index->sector);
  std::string name, title;

  name  = "dt";
  title = std::format("{} : dt",prefix);
  fBookHist->HBook1F(Hist->h_dt,name.data(),title.data(),100,-100,100,Folder);   // ns

  name  = "dt_vs_sbid";
  title = std::format("{} : dt vs SBID within the sector",prefix);
  fBookHist->HBook2F(Hist->h_dt_vs_sbid,name.data(),title.data(),100,-100,100,256,0,256,Folder);   // ns

  name  = "xc";
  title = std::format("{} : xc",prefix);
  fBookHist->HBook1F(Hist->h_xc,name.data(),title.data(),100,-10000,10000,Folder);   // mm

  name  = "dxc_0";
  title = std::format("{} : dxc[0], 1 fiber",prefix);
  fBookHist->HBook1F(Hist->h_dxc[0],name.data(),title.data(),100,-10000,10000,Folder);   // mm

  name  = "dxc_1";
  title = std::format("{} : dxc[1], both fibers",prefix);
  fBookHist->HBook1F(Hist->h_dxc[1],name.data(),title.data(),100,-10000,10000,Folder);   // mm

  name  = "sbid_vs_dxc";
  title = std::format("{} : SBID vs dxc within the sector",prefix);
  fBookHist->HBook2F(Hist->h_sbid_vs_dxc,name.data(),title.data(),100,-10000,10000,256,0,256,Folder);   // ns

  name  = "zc";
  title = std::format("{} : zc",prefix);
  fBookHist->HBook1F(Hist->h_zc,name.data(),title.data(),200,-10000,10000,Folder);   // mm

  name  = "dxloc_vs_dt10";
  title = std::format("{} : dxloc vs dt10",prefix);
  fBookHist->HBook2F(Hist->h_dxloc_vs_dt10,name.data(),title.data(),100,-50,50,100,-5000,5000,Folder);   // mm

  name  = "dxloc_vs_dt32";
  title = std::format("{} : dxloc vs dt32",prefix);
  fBookHist->HBook2F(Hist->h_dxloc_vs_dt32,name.data(),title.data(),100,-50,50,100,-5000,5000,Folder);   // mm

  return 0;
}

//-----------------------------------------------------------------------------
int TCrvTimeAnaModule::BookCrvHistograms(CrvHist_t* Hist, CrvIndex_t* Index, TFolder* Folder) {

  std::string prefix = std::format("run:{:06d} sel:{}",fRunNumber,Index->sel,Index->sector);
  std::string name, title;

  name  = "sbid";
  title = std::format("{} : sbid",prefix);
  fBookHist->HBook1F(Hist->h_sbid,name.data(),title.data(),800,0,800,Folder);   // ns

  name  = "chid";
  title = std::format("{} : chid",prefix);
  fBookHist->HBook1F(Hist->h_chid,name.data(),title.data(),3200,0,3200,Folder);   // ns

  name  = "sector";
  title = std::format("{} : sector",prefix);
  fBookHist->HBook1F(Hist->h_sector,name.data(),title.data(),30,0,30,Folder);   // ns

  name  = "ncrvc";
  title = std::format("{} : N(crvc)",prefix);
  fBookHist->HBook1F(Hist->h_ncrvc,name.data(),title.data(),100,0,100,Folder);   // ns

  name  = "crvc_mask";
  title = std::format("{} : crvc sector mask",prefix);
  fBookHist->HBook1F(Hist->h_crvc_mask,name.data(),title.data(),32,0,32,Folder);   // ns

  // CRV: assume already pre-calibrated
  name  = "dt20_vs_sbid";
  title = std::format("{} : dt20 vs SBID",prefix);
  fBookHist->HBook2F(Hist->h_dt20_vs_sbid,name.data(),title.data(),600,0,600,200,-50,50,Folder);
  
  // the other side
  name  = "dt31_vs_sbid";
  title = std::format("{} : dt31 vs SBID",prefix);
  fBookHist->HBook2F(Hist->h_dt31_vs_sbid,name.data(),title.data(),600,0,600,200,-50,50,Folder);
  
  // the other side
  name  = "dt1234_vs_sbid";
  title = std::format("{} : dt1234 vs SBID",prefix);
  fBookHist->HBook2F(Hist->h_dt1234_vs_sbid,name.data(),title.data(),600,0,600,100,-50,50,Folder);

  // charge vs chid
  name  = "pes_vs_chid";
  title = std::format("{} : PES vs ChID",prefix);
  fBookHist->HBook2F(Hist->h_pes_vs_chid,name.data(),title.data(),2400,0,2400,100,0,200,Folder);

//-----------------------------------------------------------------------------
// CRV histogram sets
//-----------------------------------------------------------------------------
  int book_crvs_histset[kNCrvsHistSets];

  for (int i=0; i<kNCrvsHistSets; i++) { book_crvs_histset[i] = 0; }

  book_crvs_histset[ 0] = 1;            // sector EX   , trk Nhits>=10
  book_crvs_histset[ 1] = 1;            // sector T1   , trk Nhits>=10
  book_crvs_histset[ 2] = 1;            // sector T2   , trk Nhits>=10
  book_crvs_histset[ 3] = 1;            // sector M1-M8, trk Nhits>=10

  book_crvs_histset[10] = 1;            // sector EX   , trk Nhits>=15
  book_crvs_histset[11] = 1;            // sector T1   , trk Nhits>=15
  book_crvs_histset[12] = 1;            // sector T2   , trk Nhits>=15
  book_crvs_histset[13] = 1;            // sector M1-M8, trk Nhits>=15
  
                                        // best resolution 
  book_crvs_histset[34] = 1;            // sector EX   , trk Nhits>=10, ro_group=4
  book_crvs_histset[37] = 1;            // sector EX   , trk Nhits>=10, ro_group=7
  
  for (int i=0; i<kNCrvsHistSets; i++) {
    // by sector
    if (book_crvs_histset[i] == 0) continue;
    std::string folder_name = std::format("sec_{:02d}",i);
    TFolder*    fol   = Folder->AddFolder(folder_name.data(),folder_name.data());
    Hist->sector[i]   = new CrvsHist_t();
    if      (i < 20) {
      Index->sector     = i % 10;
    }
    else if ((i >= 30) and (i < 40)) {
      Index->sector = 0;                // EX
    }
    
    BookCrvsHistograms(Hist->sector[i],Index,fol);
  }
  return 0;
}

//-----------------------------------------------------------------------------
int TCrvTimeAnaModule::BookTrkHistograms(TrkHist_t* Hist, TrkIndex_t* Index, TFolder* Folder) {

  std::string prefix = std::format("run:{:06d} sel:{}",fRunNumber,Index->sel);
  std::string name, title;

  name  = "nhits";
  title = std::format("{} : nhits",prefix);
  fBookHist->HBook1F(Hist->h_nhits,name.data(),title.data(),100,0,100,Folder);   // in MeV

  name  = "chi2d";
  title = std::format("{} : chi2/ndof",prefix);
  fBookHist->HBook1F(Hist->h_chi2d,name.data(),title.data(),100,0,20,Folder);   //

  name  = "t0";
  title = std::format("{} : t0",prefix);
  fBookHist->HBook1F(Hist->h_t0,name.data(),title.data(),500,0,2.5e6,Folder);   // in ns...

  name  = "dxdz";
  title = std::format("{} : dxdz",prefix);
  fBookHist->HBook1F(Hist->h_dxdz,name.data(),title.data(),200,-10,10,Folder);   //

  name  = "dydz";
  title = std::format("{} : dydz",prefix);
  fBookHist->HBook1F(Hist->h_dydz,name.data(),title.data(),200,-10,10,Folder);   // in ns...

  name  = "dt_tc";
  title = std::format("{} : dt TC",prefix);
  fBookHist->HBook1F(Hist->h_dt_tc,name.data(),title.data(),400,-100,100,Folder);   // in ns...

  name  = "ncrvc30";
  title = std::format("{} : n(CRVC) within 30 ns",prefix);
  fBookHist->HBook1F(Hist->h_ncrvc30,name.data(),title.data(),50,0,50,Folder);   // in MeV

  name  = "xcrv";
  title = std::format("{} : xcrv",prefix);
  fBookHist->HBook1F(Hist->h_xcrv,name.data(),title.data(),200,-10000,10000,Folder);   // in ns...

  name  = "zcrv";
  title = std::format("{} : zcrv",prefix);
  fBookHist->HBook1F(Hist->h_zcrv,name.data(),title.data(),200,-10000,10000,Folder);   // in ns...

  name  = "dt_crvc";
  title = std::format("{} : dt CRVC",prefix);
  fBookHist->HBook1F(Hist->h_dt_crvc,name.data(),title.data(),400,-100,100,Folder);   // in ns...

  name  = "dx_crvc";
  title = std::format("{} : dx_crvc",prefix);
  fBookHist->HBook1F(Hist->h_dx_crvc,name.data(),title.data(),200,-10000,10000,Folder);   //

  name  = "dz_crvc";
  title = std::format("{} : dx_crvc",prefix);
  fBookHist->HBook1F(Hist->h_dz_crvc,name.data(),title.data(),200,-10000,10000,Folder);   //

  name  = "dx_crvc_vs_dxdy";
  title = std::format("{} : dx_crvc vs dxdy",prefix);
  fBookHist->HBook2F(Hist->h_dx_crvc_vs_dxdy,name.data(),title.data(),100,-2,2,200,-1000,1000,Folder);

  name  = "dz_crvc_vs_dzdy";
  title = std::format("{} : dz_crvc vs dzdy",prefix);
  fBookHist->HBook2F(Hist->h_dz_crvc_vs_dzdy,name.data(),title.data(),100,-2,2,200,-1000,1000,Folder);

  return 0;
}

//-----------------------------------------------------------------------------
int TCrvTimeAnaModule::BookHistograms(Hist_t* Hist, TFolder* Folder) {

  // std::string prefix = std::format("");
  // std::string name, title;

  CaloIndex_t index;

  std::string prefix = std::format("run:{:06d}",fRunNumber);
  std::string name, title;
//-----------------------------------------------------------------------------
// Trk: tracks 
//-----------------------------------------------------------------------------
  name  = "ngtrk";
  title = std::format("{} : ngtrk",prefix);
  fBookHist->HBook1F(Hist->h_ngtrk,name.data(),title.data(),10,0,10,Folder);
  
  name  = "ntrk";
  title = std::format("{} : ntrk",prefix);
  fBookHist->HBook1F(Hist->h_ntrk,name.data(),title.data(),10,0,10,Folder);
  
  name  = "ntrk15";
  title = std::format("{} : ntrk 15",prefix);
  fBookHist->HBook1F(Hist->h_ntrk15,name.data(),title.data(),10,0,10,Folder);
  
  TrkIndex_t trk_index;

  // std::string prefix = std::format("run:{:06d}",fRunNumber);
  int book_trk_histset[kNTrkHistSets];

  for (int i=0; i<kNTrkHistSets; i++) { book_trk_histset[i] = 0; }

  book_trk_histset[0] = 1;             // all
  
  book_trk_histset[1] = 1;             // in-time tracks
  // book_trk_histset[2] = 1;             // in-time tracks DISK0
  // book_trk_histset[3] = 1;             // in-time tracks DISK1

  book_trk_histset[4] = 1;             // in-time CAL
  // book_trk_histset[5] = 1;             // in-time CAL DISK0
  // book_trk_histset[6] = 1;             // in-time CAL DISK1

  book_trk_histset[7] = 1;             // in-time trk_nhits>10 CAL
  // book_trk_histset[8] = 1;             // in-time trk_nhits>10 CAL DISK0
  // book_trk_histset[9] = 1;             // in-time trk_nhits>10 CAL DISK1

  book_trk_histset[10] = 1;             // in-time CRVC trk_nhits>10
  
  book_trk_histset[20] = 1;             // nhits>=20
  book_trk_histset[21] = 1;             // nhits>=20, at least one cluster in each disk

  for (int i=0; i<kNTrkHistSets; i++) {
    if (book_trk_histset[i] == 0) continue;
    std::string folder_name = std::format("trk_{:02d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    Hist->trk[i]   = new TrkHist_t();
    trk_index.sel  = i;
    BookTrkHistograms(Hist->trk[i],&trk_index,fol);
  }
//-----------------------------------------------------------------------------
// CRV
//-----------------------------------------------------------------------------
  CrvIndex_t crv_index;

  name  = "ncc";
  title = std::format("{} : N(close counters)",prefix);
  fBookHist->HBook1F(Hist->h_ncc,name.data(),title.data(),100,0,100,Folder);
  
  int book_crv_histset[kNTrkHistSets];

  for (int i=0; i<kNCrvHistSets; i++) { book_crv_histset[i] = 0; }

  book_crv_histset[0] = 1;             // all
  
  for (int i=0; i<kNCrvHistSets; i++) {
    if (book_crv_histset[i] == 0) continue;
    std::string folder_name = std::format("crv_{:02d}",i);
    TFolder* crv_fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! crv_fol) crv_fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    Hist->crv[i]   = new CrvHist_t();
    crv_index.sel = i;
    BookCrvHistograms(Hist->crv[i],&crv_index,crv_fol);
  }

  return 0;
}

//-----------------------------------------------------------------------------
  int TCrvTimeAnaModule::FillCrvsHistograms(CrvsHist_t* Hist, TEvdCrvCounter* Crvb, crvb_param_t* Crvbp) {
    // filling histograms: plot time differences between
    Hist->h_dt->Fill(Crvbp->dt); //}
    Hist->h_xc->Fill(Crvbp->trk_xc); //
    Hist->h_zc->Fill(Crvbp->trk_zc); //

    float dxc = Crvbp->xloc-Crvbp->trk_xc;

    //    int module_id   = Crvb->GetModuleID();
    int sector_id   = Crvb->GetSectorID();
    int sector_type = Crvb->GetSectorType();

    int first       = fEvdCrv->GetSector(Crvbp->hist_sid)->FirstCounter();
    int igr         = (Crvb->GetSbid()-first) / 32;
    
    if     (Crvbp->nxm == 1) {
      Hist->h_dxc[0]->Fill(dxc);
    }
    else if(Crvbp->nxm == 2) {
      Hist->h_dxc[1]->Fill(dxc);
      Hist->h_sbid_vs_dxc->Fill(dxc,Crvb->GetSbid()-first);
    }


    if ((sector_id > 0) and (igr > 3)) {
      GetHeaderBlock()->Print(Form("ERROR: sector_id:%i igr:%i should not be happening",
                                   sector_id,igr));
    }
    Hist->h_dt_vs_sbid->Fill(Crvbp->dt,Crvb->GetSbid()-first);
    
    if ((sector_type == 1) or (sector_type == 3)) {
      // EX or T2: long side of bars along X axis, the bars tell Z, time division - X
      Hist->h_dxloc_vs_dt10->Fill(Crvbp->dt10,dxc);
      Hist->h_dxloc_vs_dt32->Fill(Crvbp->dt32,dxc);
    }
    else {
      // T1 or M1-M8 : long side of bars along Z axis, the bars tell X, time division - Z
      // don't need these histograms, as the counters have one-sided readout
      Hist->h_dxloc_vs_dt10->Fill(Crvbp->trk_zc,Crvbp->dt10);
      Hist->h_dxloc_vs_dt32->Fill(Crvbp->trk_zc,Crvbp->dt32);
    }

    return 0;
  }
  
//-----------------------------------------------------------------------------
  int TCrvTimeAnaModule::FillCrvHistograms(CrvHist_t* Hist) {

    Hist->h_ncrvc->Fill(fNCrvc);
    Hist->h_crvc_mask->Fill(fCrvcMask);
    
    int nhc = fEvdCrv->GetNHitCounters();

    for (int ihc=0; ihc<nhc; ihc++) {
      TEvdCrvCounter* counter = fEvdCrv->GetHitCounter(ihc);
      // each counter has a parameter record
      crvb_param_t* crvbp = &fListOfCrvbParam[ihc];
       
      int sbid = counter->GetSbid();
      Hist->h_sbid->Fill(sbid);

      int sector_id = counter->GetSectorID();
      Hist->h_sector->Fill(sector_id);
      
      float dt20(1.e6), dt31(1.e6);
      int np = counter->GetNRecoPulses();
      for (int ip=0; ip<np; ip++) {
        TCrvRecoPulse* p = counter->GetRecoPulse(ip);
        int sipm = p->Sipm();
        int chid = 4*sbid+sipm;
        Hist->h_chid->Fill(chid);
        
        // sipms on the same size
        if (sipm < 2) {
          for (int ip2=0; ip2<np; ip2++) {
            TCrvRecoPulse* p2 = counter->GetRecoPulse(ip2);
            int sipm2 = p2->Sipm();
            if (sipm2 == sipm+2) {
              float dt = p2->Time()-p->Time();
              if (sipm == 0) {
                Hist->h_dt20_vs_sbid->Fill(sbid,dt);
                dt20 = dt;
              }
              else {
                Hist->h_dt31_vs_sbid->Fill(sbid,dt);
                dt31 = dt;
              }
            }
          }
        }
      }
      // here pulses are ordered by sipm
      for (int sipm=0; sipm<4; sipm++) {
        TCrvRecoPulse* p = crvbp->pulse[sipm]; 
        if (p != nullptr) {
          int chid = sbid*4+sipm;
          Hist->h_pes_vs_chid->Fill(chid,p->Pes());
        }
      }
      
      // it is a per-counter histogram
      if ((dt20 < 1.e6) and (dt31 < 1.e6)) {
        float dt1234 = dt31-dt20;
        Hist->h_dt1234_vs_sbid->Fill(sbid,dt1234);
      }
      
      // CRVS_00 - CRVS_03
      if ((fNTracks == 1) and (fNGoodTracks == 1)) {
        FillCrvsHistograms(Hist->sector[crvbp->hist_sid],counter,crvbp);
      }
      
      // CRVS10-CRVS13
      if ((fNTracks == 1) and (fNTrk15 == 1)) {
        FillCrvsHistograms(Hist->sector[crvbp->hist_sid+10],counter,crvbp);
      }

      // CRVS_34, CRVS_37: best resolution RO groups in EX
      if ((fNTracks == 1) and (fNGoodTracks == 1)) {
        if (crvbp->sector_id == 0) {
          if      (crvbp->ro_group == 4) FillCrvsHistograms(Hist->sector[34],counter,crvbp);
          else if (crvbp->ro_group == 7) FillCrvsHistograms(Hist->sector[37],counter,crvbp);
        }
      }
      
    }
    return 0;
  }
  
//-----------------------------------------------------------------------------
int TCrvTimeAnaModule::FillTrkHistograms(TrkHist_t* Hist, TStrTrack* Trk, trk_param_t* Tp) {
  // filling histograms: plot time differences between
  
  float dxdz = Trk->fNx/Trk->fNz;
  float dydz = Trk->fNy/Trk->fNz;
  
  float dxdy = dxdz/dydz;
  float dzdy = 1./dydz;
  
  Hist->h_nhits->Fill(Trk->fNHits);
  Hist->h_chi2d->Fill(Trk->fChi2/Trk->fNDof);
  Hist->h_t0->Fill(Trk->fT0);
  Hist->h_dt_tc->Fill(Tp->dtmin_tc);

  Hist->h_dxdz->Fill(dxdz);
  Hist->h_dydz->Fill(dydz);

  Hist->h_ncrvc30->Fill(Tp->ncrvc30);
  Hist->h_xcrv->Fill(Tp->xcrv);
  Hist->h_zcrv->Fill(Tp->zcrv);

  Hist->h_dt_crvc->Fill(Tp->dtmin_crvc);

  Hist->h_dx_crvc->Fill(Tp->dx_crvc);
  Hist->h_dz_crvc->Fill(Tp->dz_crvc);

  if (Tp->crvc) {
    int sector_type = Tp->crvc->SectorType();
    if ((sector_type == 1) or (sector_type == 3)) {
      // EX or T2: long side of bars along X axis, measuring Z
      Hist->h_dz_crvc_vs_dzdy->Fill(dzdy,Tp->dz_crvc);
    }
    else {
      // T1 or M1-M8 : long side of bars along Z axis, the bars tell X
      Hist->h_dx_crvc_vs_dxdy->Fill(dxdy,Tp->dx_crvc);
    }
  }
  
  return 0;
}

//-----------------------------------------------------------------------------
// need to optimize the filling time
//-----------------------------------------------------------------------------
int TCrvTimeAnaModule::FillHistograms() {
//-----------------------------------------------------------------------------
// 1. fill track histograms
//-----------------------------------------------------------------------------
  fHist->h_ntrk->Fill(fNTracks);
  fHist->h_ngtrk->Fill(fNGoodTracks);
  fHist->h_ntrk15->Fill(fNTrk15);

  fHist->h_ncc->Fill(fNCloseCounters);

  for (int i=0; i<fNTracks; i++) {
    TStrTrack* trk = fTrackBlock->Track(i);
    trk_param_t* tp = &fListOfTrkParam[i];
    FillTrkHistograms(fHist->trk[0],trk,tp);
    if (tp->intime) {
      FillTrkHistograms(fHist->trk[1],trk,tp);
    }
    
    if (tp->intime_calc) {
      FillTrkHistograms(fHist->trk[4],trk,tp);
    }

    if ((trk->NHits() > 10) and (tp->intime_calc)) {
      FillTrkHistograms(fHist->trk[7],trk,tp);
    }

    if ((trk->NHits() > 10) and (tp->intime_crvc)) {
      FillTrkHistograms(fHist->trk[10],trk,tp);
    }    

    if (trk->NHits() >= 20) {
      FillTrkHistograms(fHist->trk[20],trk,tp);
    }    
  }
//-----------------------------------------------------------------------------
// fill CRV DT20 and DT31 histograms
//-----------------------------------------------------------------------------
  FillCrvHistograms(fHist->crv[0]);

  return 0;
}

//-----------------------------------------------------------------------------
int TCrvTimeAnaModule::CalculateMissingTrkParameters() {
//-----------------------------------------------------------------------------
// extra track parameters
//-----------------------------------------------------------------------------
  fNGoodTracks = 0;
  fNTrk15      = 0;
  for (int i1=0; i1<fNTracks; i1++) {
    TStrTrack*  trk = fTrackBlock->Track(i1);
    trk_param_t* tp = &fListOfTrkParam[i1];
    tp->good_track = 0;
    if (trk->fNHits >= 10) {
      tp->good_track = 1;
      fNGoodTracks++;
      if (trk->fNHits >= 15) {
        fNTrk15++;
      }
    }
//-----------------------------------------------------------------------------
// initialize parameter record
//-----------------------------------------------------------------------------
    float trk_dxdy = trk->fNx/trk->fNy;
    float trk_dzdy = trk->fNz/trk->fNy;

    float trk_dxdz = trk->fNx/trk->fNz;
    float trk_dydz = trk->fNy/trk->fNz;
    
    tp->dtmin_tc    = 1.e6;
    tp->tc          = nullptr;
    for (int k=0; k<2; k++) {
      float dz = fDetPos.calo_pos[k][2]+fDetPos.calo_off[k][2]-(trk->fZ0+fDetPos.trk_pos[2]+fDetPos.trk_off[2]);

      tp->dtmin_calc[k] = 1.e6;
      tp->x_disk[k]     = trk->fX0 + (fDetPos.trk_pos [0]+fDetPos.trk_off [0]) + trk_dxdz*dz;
      tp->y_disk[k]     = trk->fY0 + (fDetPos.trk_pos [1]+fDetPos.trk_off [1]) + trk_dydz*dz;
      tp->calc[k]       = nullptr;
      tp->dx_calc[k]    = 1.e6;
      tp->dy_calc[k]    = 1.e6;
    }
    tp->crvc        = nullptr;
    tp->dtmin_crvc  = 1.e6;
    tp->intime      = 0;
    tp->ncrvc30     = 0;
    tp->intime_calc = 0;
    tp->intime_crvc = 0;
    tp->xcrv        = 1.e6;
    tp->zcrv        = 1.e6;
    tp->dx_crvc     = 1.e6;
    tp->dz_crvc     = 1.e6;
//-----------------------------------------------------------------------------    
// find the closest time cluster - should always be there
//-----------------------------------------------------------------------------    
    for (int i2=0; i2<fNTc; i2++) {
      TStnTimeCluster* tc = fTcBlock->TimeCluster(i2);
      float dt = trk->fT0-tc->fT0;
      if (fabs(dt) < fabs(tp->dtmin_tc)) {
        tp->dtmin_tc = dt;
        tp->tc       = tc;
      }
    }
//-----------------------------------------------------------------------------
// determine the closest CRV coincidence
// as the calibration used time clusters, look at the time cluster T0
//-----------------------------------------------------------------------------
    for (int i2=0; i2<fNCrvc; i2++) {
      TCrvCoincidenceCluster* crvc = fCrvcBlock->Cluster(i2);
      float crvc_y = crvc->Position()->Y();
      float dy     = crvc_y - trk->fY0;
      float path   = dy*sqrt(1.+trk_dxdy*trk_dxdy+trk_dzdy*trk_dzdy);
      
      float dt     = trk->fT0-path/2.997e2-crvc->StartTime();

      if (fabs(dt) < 30) {
        tp->ncrvc30 += 1;
      }
      
      if (fabs(dt) < fabs(tp->dtmin_crvc)) {
        tp->dtmin_crvc = dt;
        tp->crvc       = crvc;
      }
    }
    
    if (tp->crvc) {
      if (fabs(tp->dtmin_crvc) <  30) {
        tp->intime_crvc = 1;
      }
//-----------------------------------------------------------------------------
// transform track coordinates to global coordinate system
// constants from Offline/Mu2eG4/geom/geom_common_extracted_v04.txt
//-----------------------------------------------------------------------------
      float y_crvc = tp->crvc->Position()->Y() + fDetPos.crv_off[1];  // -80.;//-4280; // for kicks
      float dy     = y_crvc-trk->fY0;

      tp->xcrv     = trk->fX0+(trk->fNx/trk->fNy)*dy; // track coordinates at Y_CRVC
      tp->zcrv     = trk->fZ0+(trk->fNz/trk->fNy)*dy;
      
      // X(CRVC) is defined in the global coordinate system, transform that
      // to the local reference frame of the tracker

      float x_crvc = tp->crvc->Position()->X()-(fDetPos.trk_pos[0]+fDetPos.trk_off[0]);
      tp->dx_crvc  = tp->xcrv-x_crvc;

      float z_crvc = tp->crvc->Position()->Z()-(fDetPos.trk_pos[2]+fDetPos.trk_off[2]);
      tp->dz_crvc  = tp->zcrv-z_crvc;
    }
    
    if ((fabs(tp->dtmin_tc  ) <  30) and
        (fabs(tp->dtmin_crvc) <  30) and
        (trk->fNHits          >= 10)     ) {
      tp->intime = 1;
    }
  }
  
  return 0;
}

//-----------------------------------------------------------------------------
int TCrvTimeAnaModule::CalculateMissingParameters() {

  CalculateMissingTrkParameters();

//-----------------------------------------------------------------------------
// deal with the CRV
//-----------------------------------------------------------------------------
  float dxdt_prop[4] = {72., 0, 78., 0};
  
  fEvdCrv->Clear();
  fListOfCrvbParam.resize(fNCrvp);

  for (int i=0; i<fNCrvp; i++) {
    TCrvRecoPulse* p = fCrvpBlock->Pulse(i);
    fEvdCrv->AddRecoPulse(p);
  }

  fCrvcMask = 0;
  for (int i=0; i<fNCrvc; i++) {
    // TCrvCoincidenceCluster* crvc = fCrvcBlock->Cluster(i);
    int ind = fCrvcBlock->ClusterPulseIndex(i,0);
    TCrvRecoPulse* p = fCrvpBlock->Pulse(ind);            // each cluster chould have at least one
    int sector = TCrvNumerology::Instance()->GetSector(p->Sbid());
    fCrvcMask |= (1<<sector);
  }

  fNCloseCounters = 0;
  int nhc = fEvdCrv->GetNHitCounters();
  // for each bar, extrapolate track to its Y

  for (int ihc=0; ihc<nhc; ihc++) {
    TEvdCrvCounter* counter = fEvdCrv->GetHitCounter(ihc);
    crvb_param_t* crvbp     = &fListOfCrvbParam[ihc];
    crvbp->sector_id        = counter->GetSectorID();
    crvbp->sector_type      = counter->GetSectorType();
    // hist_sid is used only to merge M1-M8 into one set
    crvbp->hist_sid         = crvbp->sector_id; // default
    crvbp->ro_group         = counter->GetSbid() / 32;
      
    if (crvbp->sector_id > 2) {                 // all M's go into group 3
      crvbp->hist_sid = 3;
    }
    crvbp->t                = 1.e6;
    crvbp->dt               = 1.e6;
    crvbp->trk_xc           = 1.e6;
    crvbp->trk_zc           = 1.e6;
    crvbp->dt10             = 1.e6;
    crvbp->dt32             = 1.e6;
    
    for (int sipm=0; sipm<4; sipm++) {
      crvbp->pulse[sipm] = nullptr;
    }

    const double* pos = counter->GetGeoMatrix()->GetTranslation();
    float crv_y  = pos[1]+fDetPos.crv_off[1];
    
    int nh = counter->GetNRecoPulses();

    // assume only one good track per event

    for (int it=0; it<fNTracks; it++) {
      TStrTrack*   trk = fTrackBlock->Track(it);
      trk_param_t* tp  = &fListOfTrkParam[it];
      if (not tp->good_track)  continue;
//-----------------------------------------------------------------------------
// coordinates
//-----------------------------------------------------------------------------
      float trk_y  = trk->fY0-fDetPos.trk_pos[1]+fDetPos.trk_off[1];
      float dy     = crv_y-trk_y;
      
      // track time on the counter assume all cosmic particls come from the top - not necessarily true
      float dxdy   = trk->fNx/trk->fNy;
      float dzdy   = trk->fNz/trk->fNy;
      float path   = dy*sqrt(1+dxdy*dxdy+dzdy*dzdy);
      crvbp->trk_t = trk->fT0-path/2.997e2; // dy is in mm, so the speed of lignt - in mm/ns

      // track coordinates in global RS at Y_CRVC
      crvbp->trk_x = trk->fX0+fDetPos.trk_pos[0]+fDetPos.trk_off[0]+(trk->fNx/trk->fNy)*dy;
      crvbp->trk_z = trk->fZ0+fDetPos.trk_pos[2]+fDetPos.trk_off[2]+(trk->fNz/trk->fNy)*dy;

      crvbp->trk_xc = crvbp->trk_x-(pos[0]+fDetPos.crv_off[0]);
      crvbp->trk_zc = crvbp->trk_z-(pos[2]+fDetPos.crv_off[2]);
//-----------------------------------------------------------------------------
// timing: figure out the time closest to the track
//-----------------------------------------------------------------------------
      TCrvRecoPulse* pp[4] = {};
      
      for (int ih=0; ih<nh; ih++) {
        TCrvRecoPulse* p = counter->GetRecoPulse(ih);
        if (fabs(p->Time()-crvbp->trk_t) > 80.) continue;
        int sipm = p->Sipm();
        // if there is more than one pulse, use the highest one
        if (pp[sipm] == nullptr) {
          pp[sipm] = p;
        }
        else {
          if (p->Pes() > pp[sipm]->Pes()) {
            pp[sipm] = p;
          }
        }
      }

      float t(0);
      int nap(0);                       // number of "active" (or good) pulses
      for (int sipm=0; sipm<4; sipm++) {
        crvbp->pulse[sipm] = pp[sipm];
        if (pp[sipm]) {
          t   += pp[sipm]->Time();
          nap += 1;
        }
      }

      if (nap > 0) {
        crvbp->t  = t/nap;
        crvbp->dt = crvbp->t-crvbp->trk_t;
      }
//-----------------------------------------------------------------------------
// calculate coordinate along the counter
//-----------------------------------------------------------------------------
      float xloc = 0;
      crvbp->nxm = 0;
 
      float vprop    = dxdt_prop[crvbp->sector_id];
      if (pp[0] and pp[1]) {
        crvbp->dt10 = pp[1]->Time()-pp[0]->Time();
        xloc += vprop*crvbp->dt10;  // 79 = 75+3
        crvbp->nxm++;
      }
      if (pp[2] and pp[3]) {
        crvbp->dt32 = pp[3]->Time()-pp[2]->Time();
        xloc += vprop*crvbp->dt32;
        crvbp->nxm++;
      }
                                        // swap the sign
      if (crvbp->nxm > 0) {
        crvbp->xloc = -xloc/crvbp->nxm;
      }

      if (nap == 4) {
        fNCloseCounters++;
      }
    }
  }
  
  return 0;
}

//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TCrvTimeAnaModule::BeginJob() {
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
  RegisterDataBlock("ComboHitBlock"    ,"TComboHitBlock"      ,&fChBlock          );
  RegisterDataBlock("TrackBlock"       ,"TStrTrackBlock"      ,&fTrackBlock       );
 
  return 0;
}


//_____________________________________________________________________________
int TCrvTimeAnaModule::BeginRun() {
  fRunNumber = GetHeaderBlock()->RunNumber();
  TStntuple::Init(fRunNumber);
  fTpm->Init(fRunNumber);
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
int TCrvTimeAnaModule::Event(int IEntry) {

  fTcBlock->GetEntry(IEntry);
  
  fCrvpBlock->GetEntry(IEntry);
  fCrvcBlock->GetEntry(IEntry);
  fCrvdBlock->GetEntry(IEntry);
  
  fTrackBlock->GetEntry(IEntry);
//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
  fNCrvc         = fCrvcBlock->NClusters();
  fNCrvp         = fCrvpBlock->NPulses();
  fNTc           = fTcBlock->NTimeClusters();
  //  fNCaloClusters = fCaloClusterBlock->NClusters();
  fNTracks          = fTrackBlock->NTracks();
  //  fNCalh         = fCaloHitBlock->NHits();

  CalculateMissingParameters();
  
  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TCrvTimeAnaModule::Debug() {

  if (GetDebugBit(0) == 1) {
    GetHeaderBlock()->Print(Form("bit_00: all events"));
  }

  if (GetDebugBit(3) == 1) {
    int found(0);
    for (int i=0; i<fNTracks; i++) {
      trk_param_t* tp = &fListOfTrkParam[i];
      if (tp->ncrvc30 >=4) {
        found = 1;
        break;
      }
    }
    if (found) {
      GetHeaderBlock()->Print(Form("bit_03"));
      fTrackBlock->Print();
     fCrvcBlock->Print();
    }
  }
}

//_____________________________________________________________________________
int TCrvTimeAnaModule::EndJob() {
  return 0;
}


//-----------------------------------------------------------------------------
// do that for all ROCs and all FEBs
//-----------------------------------------------------------------------------
/*
# roc  feb   fit(crvp.time-tc.t0)  chi2dof
   1    1           530.724         3.004
   1    2           530.498         2.150
*/
int TCrvTimeAnaModule::FitCaloTimeOffsets(float TMin, float TMax) {
  /*
  TH2F* h2 = nullptr; // fHist->h_dt_vs_sipmid;

  int nch = kNCaloChannels;
  for (int i=0; i<nch; i++) {
    
    std::string hpname = std::format("hpx_{:04d}",i);
    TH1D* hp = h2->ProjectionY(hpname.data(),i+1,i+1);
    fit_result_t* fr = &fFr[i];
    FitHistogram(hp,fr,TMin,TMax,10);
  }
  
  // done fitting, write results

  std::ofstream os("CaloTimeCalib_corr.txt");

  os << std::format("TABLE CalTimeCalib {}\n",fRunNumber);
  os << std::format("# corrections to  CalTimeCalib cid= XXX, to be subtracted \n");
  os << std::format("# sipmid     dT      sigma      chi2 \n");
  
  for (int i=0; i<nch; i++) {
    
    // to do no harm if no fit, dt should be initialized to zero 
    float dt{0}; 
    if (fFr[i].chi2dof > 0) {
      dt = fFr[i].p[1];
    }

    os << std::format("{:5}   {:8.3f}  {:8.3f}  {:8.3f}\n",i,dt,fFr[i].p[2],fFr[i].chi2dof);
  }
  os.close();
  */
  return 0;
}

//-----------------------------------------------------------------------------
int TCrvTimeAnaModule::PrintTracks() {
  std::cout << std::format(" i      T0         Z0       Nx     Ny    Nz    Chi2D     Xc[0]      Yc[0]      Xc[1]      Yc[1]    DxCal[0]    DyCal[0]    DxCal[1]   DyCal[1]\n");
  for (int i=0; i<fNTracks; i++) {
    TStrTrack*   trk = fTrackBlock->Track(i);
    trk_param_t* tp  = &fListOfTrkParam.at(i);
    std::cout << std::format("{:2d} {:10.3f} {:10.3f} {:6.3f} {:6.3f} {:6.3f} {:7.2f}",
                             i,trk->fT0,trk->fZ0,trk->fNx,trk->fNy,trk->fNz,trk->fChi2/trk->fNDof);
    
    std::cout << std::format(" {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f}\n",
                             tp->x_disk [0],tp->y_disk [0],tp->x_disk[1],tp->y_disk[1],
                             tp->dx_calc[0],tp->dy_calc[0],tp->dx_calc[1],tp->dy_calc[1]);
  }
  return 0;
}


  //-----------------------------------------------------------------------------
  int TCrvTimeAnaModule::PrintCrvBars() {

    int ns = fEvdCrv->GetNSectors();

    bool banner_printed = false;
    
    for (int is=0; is<ns; is++) {
      TEvdCrvSector* sector = fEvdCrv->GetSector(is);
      if (sector->GetNHits() == 0) continue;
      int nm = sector->GetNModules();
      for (int im=0; im<nm; im++) {
        TEvdCrvModule* module = sector->GetModule(im);
        if (module->GetNHits() == 0) continue;
        int nl = module->GetNLayers();
        for (int il=0; il<nl; il++) {
          TEvdCrvLayer* layer = module->GetLayer(il);
          if (layer->GetNHits() == 0) continue;
          int nc = layer->GetNCounters();
          for (int ic=0; ic<nc; ic++) {
            TEvdCrvCounter* counter = layer->GetCounter(il);
            int nh = counter->GetNRecoPulses();
            if (nh == 0) continue;
            for (int ih=0; ih<nh; ih++) {
              TCrvRecoPulse* p = counter->GetRecoPulse(ih);
              if (not banner_printed) {
                p->Print("banner");
                banner_printed = true;
              }
              p->Print("data");
            }
          }
        }
      }
    }
            
    return 0;
  }


//-----------------------------------------------------------------------------
  int TCrvTimeAnaModule::InitEvd() {
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

    vm->AddDataBlock("TrackBlock"       ,fTrackBlock);
//-----------------------------------------------------------------------------
// CRV-only view
//-----------------------------------------------------------------------------
    TEvdCrvNode* node = new TEvdCrvNode("crv_node");
    TEvdCrvView* view = new TEvdCrvView("crv_view");
    view->AddNode(node);
    vm->AddView(view);
//-----------------------------------------------------------------------------
// 3D detector view
//-----------------------------------------------------------------------------
    TEvdDetNode* det_node = new TEvdDetNode("det_node");
    TEvdView*    det_view = new TEvdView   ("det_view");
    det_view->AddNode(det_node);
    vm->AddView(det_view);

    return rc;
  }
  
//-----------------------------------------------------------------------------
  int TCrvTimeAnaModule::DisplayEvent() {
    int rc(0);
    auto vm = TEvdManager::Instance();

    vm->DisplayEvent();

    return rc;
  }
  
//-----------------------------------------------------------------------------
void TCrvTimeAnaModule::PlotWaveform(int ChannelID, bool NewCanvas) {

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
