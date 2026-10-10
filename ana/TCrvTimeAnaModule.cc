//////////////////////////////////////////////////////////////////////////////
// use of tmp:
//
// use of debug bits: bits 0-2 are reserved
//  0  : all events
//  1  : passed events
//  2  : rejected events
//  3  : N(Calc  disk 0 ) > 0 and N(calc disk1 > 0)
//  4  : N(Calc  disk 0 ) > 0 and N(calc disk1 > 0) and a track with N>=20 hits
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
#include "murat/gui/TEvdCaloView.hh"
#include "murat/gui/TEvdCaloNode.hh"
#include "murat/gui/TEvdDetNode.hh"
#include "murat/gui/TEvdManager.hh"

#include "murat/ana/TCrvTimeAnaModule.hh"


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

  fEvdCrv = new TEvdCrv();
  fEvdCrv->InitGeometry("Offline/Mu2eG4/geom/geom_common_extracted.txt");
}

//-----------------------------------------------------------------------------
TCrvTimeAnaModule::~TCrvTimeAnaModule() {
}

//-----------------------------------------------------------------------------
int TCrvTimeAnaModule::BookTrkHistograms(TrkHist_t* Hist, TrkIndex_t* Index, TFolder* Folder) {

  std::string prefix = std::format("run:{:06d} sel:{}",fRunNumber,Index->sel);
  std::string name, title;

  name  = "nhits";
  title = std::format("{} : nhits",prefix);
  fBookHist->HBook1F(Hist->h_nhits,name.data(),title.data(),100,0,100,Folder);   // in MeV

  name  = "chi2d";
  title = std::format("{} : chi2/ndod",prefix);
  fBookHist->HBook1F(Hist->h_chi2d,name.data(),title.data(),100,0,20,Folder);   // in ns...

  name  = "t0";
  title = std::format("{} : t0",prefix);
  fBookHist->HBook1F(Hist->h_t0,name.data(),title.data(),500,0,2.5e6,Folder);   // in ns...

  name  = "dxdz";
  title = std::format("{} : dxdz",prefix);
  fBookHist->HBook1F(Hist->h_dxdz,name.data(),title.data(),200,-10,10,Folder);   // in ns...

  name  = "dydz";
  title = std::format("{} : dydz",prefix);
  fBookHist->HBook1F(Hist->h_dydz,name.data(),title.data(),200,-10,10,Folder);   // in ns...

  name  = "dt_tc";
  title = std::format("{} : dt TC",prefix);
  fBookHist->HBook1F(Hist->h_dt_tc,name.data(),title.data(),400,-100,100,Folder);   // in ns...

  // for (int k=0; k<2; k++) {
  //   name  = std::format("dt_calc_{}",k);
  //   title = std::format("{} : dt CALC disk:{}",prefix,k);
  //   fBookHist->HBook1F(Hist->h_dt_calc[k],name.data(),title.data(),400,-100,100,Folder);   // in ns...

  //   name  = std::format("dx_calc_{}",k);
  //   title = std::format("{} : dx_calc = x(trk)-x(calc) disk:{}",prefix,k);
  //   fBookHist->HBook1F(Hist->h_dx_calc[k],name.data(),title.data(),200,-1000,1000,Folder);   // in ns...

  //   name  = std::format("dy_calc_{}",k);
  //   title = std::format("{} : dy_calc = y(trk)-y(calc) disk:{}",prefix,k);
  //   fBookHist->HBook1F(Hist->h_dy_calc[k],name.data(),title.data(),200,-1000,1000,Folder);   // in ns...

  //   name  = std::format("dx_calc_vs_dxdz_{}",k);
  //   title = std::format("{} : dx_calc vs dxdz disk:{}",prefix,k);
  //   fBookHist->HBook2F(Hist->h_dx_calc_vs_dxdz[k],name.data(),title.data(),200,-1,1,200,-1000,1000,Folder);

  //   name  = std::format("dy_calc_vs_dydz_{}",k);
  //   title = std::format("{} : dy_calc vs dydz disk:{}",prefix,k);
  //   fBookHist->HBook2F(Hist->h_dy_calc_vs_dydz[k],name.data(),title.data(),200,-1.0,1.0,200,-1000,1000,Folder);
  // }
  
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

  // name  = "dt_vs_sipmid";
  // title = std::format("{} : dt vs SIPM ID",prefix);
  // fBookHist->HBook2F(Hist->h_dt_vs_sipmid,name.data(),title.data(),2800,0,2800,1000,-1000,1000,Folder);

  // name  = "nsipms_vs_cid";
  // title = std::format("{} : nsipms vs crystal ID",prefix);
  // fBookHist->HBook2F(Hist->h_nsipms_vs_cid,name.data(),title.data(),1400,0,1400,3,0,3,Folder);

  // name  = "sipmid";
  // title = std::format("{} : Sipm ID",prefix);
  // fBookHist->HBook1F(Hist->h_sipmid,name.data(),title.data(),2800,0,2800,Folder);

  // name  = "n2_vs_n1";
  // title = std::format("{} : N1:N1 calh E>10",prefix);
  // fBookHist->HBook2F(Hist->h_n2_vs_n1,name.data(),title.data(),20,0,20,20,0,20,Folder);

  // name  = "ntrk_0";
  // title = std::format("{} : ntrk[0]",prefix);
  // fBookHist->HBook1F(Hist->h_ntrk[0],name.data(),title.data(),10,0,10,Folder);

  // name  = "ntrk_1";
  // title = std::format("{} : ntrk[1]",prefix);
  // fBookHist->HBook1F(Hist->h_ntrk[1],name.data(),title.data(),10,0,10,Folder);

  // CRV: assume already pre-calibrated
  name  = "dt20_vs_sbid";
  title = std::format("{} : dt20 vs SBID",prefix);
  fBookHist->HBook2F(Hist->h_dt20_vs_sbid,name.data(),title.data(),1000,0,800,200,-50,50,Folder);
  // the other side
  name  = "dt31_vs_sbid";
  title = std::format("{} : dt31 vs SBID",prefix);
  fBookHist->HBook2F(Hist->h_dt31_vs_sbid,name.data(),title.data(),1000,0,800,200,-50,50,Folder);
  
//-----------------------------------------------------------------------------
// CALC: calorimeter clusters 
//-----------------------------------------------------------------------------
  // int book_calc_histset[kNCalcHistSets];
  // for (int i=0; i<kNCalcHistSets; i++) { book_calc_histset[i] = 0; }

  // book_calc_histset[0] = 1;             // all
//-----------------------------------------------------------------------------
// Trk: tracks 
//-----------------------------------------------------------------------------
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

  // for (int k=0; k<2; k++) {
  //   if (Tp->calc[k]) {
  //     Hist->h_dt_calc[k]->Fill(Tp->dtmin_calc[k]);
  //     Hist->h_dx_calc[k]->Fill(Tp->dx_calc[k]);
  //     Hist->h_dy_calc[k]->Fill(Tp->dy_calc[k]);

  //     Hist->h_dx_calc_vs_dxdz[k]->Fill(dxdz,Tp->dx_calc[k]);
  //     Hist->h_dy_calc_vs_dydz[k]->Fill(dydz,Tp->dy_calc[k]);
  //   }
  // }

  Hist->h_xcrv->Fill(Tp->xcrv);
  Hist->h_zcrv->Fill(Tp->zcrv);

  Hist->h_dt_crvc->Fill(Tp->dtmin_crvc);

  Hist->h_dx_crvc->Fill(Tp->dx_crvc);
  Hist->h_dz_crvc->Fill(Tp->dz_crvc);

  if (Tp->crvc) {
    int s_type = Tp->crvc->SectorType();
    if ((s_type == 1) or (s_type == 3)) {
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
  // filling histograms: plot time differences between

//-----------------------------------------------------------------------------
// fill track histograms
//-----------------------------------------------------------------------------
  fHist->h_ntrk[0]->Fill(fNTrk);

  for (int i=0; i<fNTrk; i++) {
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

    // if ((fNCcDisk[0] >= 1) and (fNCcDisk[1] >= 1) and (trk->NHits() >= 20)) {
    //   FillTrkHistograms(fHist->trk[21],trk,tp);
    // }    
  }
  
//-----------------------------------------------------------------------------
// fill CRV DT20 and DT31 histograms
//-----------------------------------------------------------------------------
  int ns = fEvdCrv->GetNSectors();

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
          TEvdCrvCounter* counter = layer->GetCounter(ic);
          int nh = counter->GetNRecoPulses();
          if (nh == 0) continue;

          
          // hits in the counter - plot dt_20 and dt_31
          std::vector<TCrvRecoPulse*> hits[4];
          
          for (int ih=0; ih<nh; ih++) {
            TCrvRecoPulse* p = counter->GetRecoPulse(ih);
            int sipm = p->Sipm();
            hits[sipm].push_back(p);
          }

          int nhh[4];
          for (int i=0; i<4; i++) {
            nhh[i] = hits[i].size();
          }
          // 0-2 side
          if ((nhh[0] > 0) and (nhh[2] > 0)) {
            for (int i0=0; i0<nhh[0]; i0++) {
              TCrvRecoPulse* p0 = hits[0][i0];
              int sbid = p0->Sbid();
              if (p0->Pes() < 10) continue;
              for (int i2=0; i2<nhh[2]; i2++) {
                TCrvRecoPulse* p2 = hits[2][i2];
                if (p2->Pes() < 10) continue;
                float dt = p2->Time()-p0->Time();
                fHist->h_dt20_vs_sbid->Fill(sbid,dt);
              }
            }
          }
          // 1-3 side
          if ((nhh[1] > 0) and (nhh[3] > 0)) {
            for (int i1=0; i1<nhh[1]; i1++) {
              TCrvRecoPulse* p1 = hits[1][i1];
              int sbid = p1->Sbid();
              if (p1->Pes() < 10) continue;
              for (int i3=0; i3<nhh[3]; i3++) {
                TCrvRecoPulse* p3 = hits[3][i3];
                if (p3->Pes() < 10) continue;
                float dt = p3->Time()-p1->Time();
                fHist->h_dt31_vs_sbid->Fill(sbid,dt);
              }
            }
          }
        }
      }
    }
  }

  return 0;
}

//-----------------------------------------------------------------------------
int TCrvTimeAnaModule::CalculateMissingTrkParameters() {
//-----------------------------------------------------------------------------
// global positioning constants
// tracker.z0: 22935.0
// tracker.
//-----------------------------------------------------------------------------
    float crv_time_offset = 0.; // 21; // today

    double crv_off[3] = {    0., -145.0, 0.0 };
      
    double trk_pos[3] = {-3904., 0., 22935.0 }; // nominal
    //    double trk_off[3] = {    0., 0., -1171.0 }; // offset to be added
    double trk_off[3] = {    0., 0.,     0.0 }; // offset to be added

    double calo_pos[2][3] = {           // nominal
      -3904., 0., 23000.0+2383.-64.,        // 25842;
      -3904., 0., 23000.0+3517.-64.
    };
    
    double calo_off[2][3] = {           // offset to be added
      0.    , 0.,     0.,
      0.    , 0.,     0.
    };
//-----------------------------------------------------------------------------
// extra track parameters
//-----------------------------------------------------------------------------
  for (int i1=0; i1<fNTrk; i1++) {
    TStrTrack*  trk = fTrackBlock->Track(i1);
    trk_param_t* tp = &fListOfTrkParam[i1];
//-----------------------------------------------------------------------------
// initialize parameter record
//-----------------------------------------------------------------------------
    tp->dtmin_tc    = 1.e6;
    tp->tc          = nullptr;
    for (int k=0; k<2; k++) {
      float dz          = calo_pos[k][2]+calo_off[k][2]-(trk->fZ0+trk_pos[2]+trk_off[2]);

      tp->dtmin_calc[k] = 1.e6;
      tp->x_disk[k]     = trk->fX0     + (trk_pos [0]+trk_off [0]) + (trk->fNx/trk->fNz)*dz;
      tp->y_disk[k]     = trk->fY0     + (trk_pos [1]+trk_off [1]) + (trk->fNy/trk->fNz)*dz;
      tp->calc[k]       = nullptr;
      tp->dx_calc[k]    = 1.e6;
      tp->dy_calc[k]    = 1.e6;
    }
    tp->crvc        = nullptr;
    tp->dtmin_crvc  = 1.e6;
    tp->intime      = 0;
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
      float dt = tp->tc->T0()-(crvc->StartTime()-crv_time_offset);
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
      float y_crvc = tp->crvc->Position()->Y() + crv_off[1];  // -80.;//-4280; // for kicks
      float dy     = y_crvc-trk->fY0;

      tp->xcrv     = trk->fX0+(trk->fNx/trk->fNy)*dy; // track coordinates at Y_CRVC
      tp->zcrv     = trk->fZ0+(trk->fNz/trk->fNy)*dy;
      
      // X(CRVC) is defined in the global coordinate system, transform that
      // to the local reference frame of the tracker

      float x_crvc = tp->crvc->Position()->X()-(trk_pos[0]+trk_off[0]);
      tp->dx_crvc  = tp->xcrv-x_crvc;

      float z_crvc = tp->crvc->Position()->Z()-(trk_pos[2]+trk_off[2]);
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
  fEvdCrv->Clear();

  for (int i=0; i<fNCrvp; i++) {
    TCrvRecoPulse* p = fCrvpBlock->Pulse(i);
    fEvdCrv->AddRecoPulse(p);
  }

  // for each bar, extrapolate track to its Y

  for (int it=0; it<fNTracks; it++) {
    TStrTrack*   trk = fTrackBlock->Track(it);
    trk_param_t* tp  = &fListOfTrkParam[it];
    
    int nhc = fEvdCrv->GetNHitCounters();
    for (int ihc=0; ihc<nhc; ihc++) {
      TEvdCrvCounter*    counter = fEvdCrv->GetHitCounter(ihc);
      crv_counter_par_t* crvbp   = &fListOfCrvbParam[ihc];
      // extrapolate track
      double* pos = counter->GetTranslation();
      float dy     = p[1]-trk->fY0;

      crvbp->xtrk = trk->fX0+(trk->fNx/trk->fNy)*dy; // track coordinates at Y_CRVC
      crvbp->ztrk = trk->fZ0+(trk->fNz/trk->fNy)*dy;

      if ((s_type == 1) or (s_type == 3)) {
        // EX or T2: long side of bars along X axis, measuring Z
        crvbp->dz = pos[2]-crvbp->ztrk;
        // for dx, need to do the time/charge division (assuming two-sided readout)
      }
      else {
        // T1 or M1-M8 : long side of bars along Z axis, the bars tell X
        crvbp->dx = pos[0]-crvbp->xtrk;
        // for dx, need to do the time/charge division (assuming two-sided readout)
      }

      int nh = counter->GetNRecoPulses();
      // two types of counter orientation - along Z and along X
      crvbp->d
    
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
  
  fTrackBlock->GetEntry(IEntry);
//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
  fNCrvc         = fCrvcBlock->NClusters();
  fNCrvp         = fCrvpBlock->NPulses();
  fNTc           = fTcBlock->NTimeClusters();
  //  fNCaloClusters = fCaloClusterBlock->NClusters();
  fNTrk          = fTrackBlock->NTracks();
  //  fNCalh         = fCaloHitBlock->NHits();

  CalculateMissingParameters();
  
  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TCrvTimeAnaModule::Debug() {

  // if (GetDebugBit(3) == 1) {
  //   if ((fNCcDisk[0] > 0) and (fNCcDisk[1] > 0)) {
  //     GetHeaderBlock()->Print(Form("NCcDisk[0]:%2i fNCcDisk[1]:%2i",
  //                                  fNCcDisk[0],fNCcDisk[1]));
  //     fCaloClusterBlock->Print();
  //   }
  // }

  // if (GetDebugBit(4) == 1) {
  //   if ((fNCcDisk[0] > 0) and (fNCcDisk[1] > 0)) {
  //     if (fNTrk >= 1) {
  //       int nh = fTrackBlock->Track(0)->NHits();
  //       if (nh >= 20) {
  //         GetHeaderBlock()->Print(Form("NCcDisk[0]:%2i fNCcDisk[1]:%2i N_trk_hits = %3d",
  //                                      fNCcDisk[0],fNCcDisk[1],nh));
  //         fCaloClusterBlock->Print();
  //       }
  //     }
  //   }
  // }
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
  
  return 0;
}

//-----------------------------------------------------------------------------
int TCrvTimeAnaModule::PrintTracks() {
  std::cout << std::format(" i      T0         Z0       Nx     Ny    Nz    Chi2D     Xc[0]      Yc[0]      Xc[1]      Yc[1]    DxCal[0]    DyCal[0]    DxCal[1]   DyCal[1]\n");
  for (int i=0; i<fNTrk; i++) {
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
  int TCrvTimeAnaModule::InitEvd() {
    int rc(0);

    auto vm = TEvdManager::Instance();
    vm->InitGeometry(); // if not initialized
    
    vm->AddDataBlock("HeaderBlock"      ,GetHeaderBlock());
    vm->AddDataBlock("CaloHitBlock"     ,fCaloHitBlock);
    vm->AddDataBlock("CaloRecoDigiBlock",fCaloRecoDigiBlock);
    vm->AddDataBlock("CaloCLusterBlock" ,fCaloClusterBlock);
    
    vm->AddDataBlock("CrvpBlock"        ,fCrvpBlock);
    vm->AddDataBlock("CrvcBlock"        ,fCrvcBlock);

    vm->AddDataBlock("TrackBlock"       ,fTrackBlock);
//-----------------------------------------------------------------------------
// calorimeter view
//-----------------------------------------------------------------------------
    TEvdCaloNode* node = new TEvdCaloNode("calo_node");
    TEvdCaloView*    view = new TEvdCaloView("calo_view");
    view->AddNode(node);
    vm->AddView(view);
//-----------------------------------------------------------------------------
// 3D detector view
//-----------------------------------------------------------------------------
    TEvdDetNode* det_node = new TEvdDetNode("det_node");
    TEvdView* det_view       = new TEvdView("det_view");
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
  
}
