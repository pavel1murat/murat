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
// Mu2e offline includes
//-----------------------------------------------------------------------------
#include "murat/ana/TCaloTimeAnaModule.hh"


ClassImp(murat::TCaloTimeAnaModule)

namespace murat {
//-----------------------------------------------------------------------------
TCaloTimeAnaModule::TCaloTimeAnaModule(const char* name, const char* title):
  TStnModule(name,title)
{
  TH1::AddDirectory(0);
//-----------------------------------------------------------------------------
// detector numerology
//-----------------------------------------------------------------------------
//  fTpm            = TrkPanelMap::Instance();
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
}

//-----------------------------------------------------------------------------
TCaloTimeAnaModule::~TCaloTimeAnaModule() {
}

//-----------------------------------------------------------------------------
int TCaloTimeAnaModule::BookCalcHistograms(CalcHist_t* Hist, CaloIndex_t* Index, TFolder* Folder) {

  std::string prefix = std::format("run:{:06d} sel:{} disk:{:02d} crate:{}",
                                   fRunNumber,Index->sel, Index->disk, Index->crate);
  std::string name, title;

  name  = "edep";
  title = std::format("{} : edep",prefix);
  fBookHist->HBook1F(Hist->h_edep,name.data(),title.data(),100,0,1000,Folder);   // in MeV

  return 0;
}

//-----------------------------------------------------------------------------
int TCaloTimeAnaModule::BookCalhHistograms(CalhHist_t* Hist, CaloIndex_t* Index, TFolder* Folder) {

  std::string prefix = std::format("run:{:06d} sel:{} disk:{:02d} crate:{}",
                                   fRunNumber,Index->sel, Index->disk, Index->crate);
  std::string name, title;

  // name  = "ph";
  // title = std::format("{} : ph",prefix);
  // fBookHist->HBook1F(Hist->h_ph,name.data(),title.data(),100,0,1000,Folder);   // in us...



  return 0;
}


//-----------------------------------------------------------------------------
int TCaloTimeAnaModule::BookHistograms(Hist_t* Hist, TFolder* Folder) {

  // std::string prefix = std::format("");
  // std::string name, title;

  CaloIndex_t index;

  std::string prefix = std::format("run:{:06d}",fRunNumber);
  std::string name, title;

  name  = "dt_vs_sipmid";
  title = std::format("{} : dt vs SIPM ID",prefix);
  fBookHist->HBook2F(Hist->h_dt_vs_sipmid,name.data(),title.data(),2800,0,2800,1000,-1000,1000,Folder);

  name  = "dt10_vs_crystal";
  title = std::format("{} : dt10 vs crystal ID",prefix);
  fBookHist->HBook2F(Hist->h_dt10_vs_crystal,name.data(),title.data(),1400,0,1400,1000,-500,500,Folder);

  name  = "nsipms_vs_cid";
  title = std::format("{} : nsipms vs crystal ID",prefix);
  fBookHist->HBook2F(Hist->h_nsipms_vs_cid,name.data(),title.data(),1400,0,1400,3,0,3,Folder);

  name  = "sipmid";
  title = std::format("{} : Sipm ID",prefix);
  fBookHist->HBook1F(Hist->h_sipmid,name.data(),title.data(),2800,0,2800,Folder);

  name  = "n2_vs_n1";
  title = std::format("{} : N1:N1 calh E>10",prefix);
  fBookHist->HBook2F(Hist->h_n2_vs_n1,name.data(),title.data(),20,0,20,20,0,20,Folder);
//-----------------------------------------------------------------------------
// calorimeter hits 
//-----------------------------------------------------------------------------
  int book_calh_histset[kNCalhHistSets];

  for (int i=0; i<kNCalhHistSets; i++) { book_calh_histset[i] = 0; }

  book_calh_histset[0] = 1;             // all

  for (int i=0; i<kNCalhHistSets; i++) {
    if (book_calh_histset[i] == 0) continue;
    std::string folder_name = std::format("calh_{:02d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    Hist->calh[i] = new CalhHist_t();
    index.sel = i;
    BookCalhHistograms(Hist->calh[i],&index,fol);
  }

//-----------------------------------------------------------------------------
// CALC: calorimeter clusters 
//-----------------------------------------------------------------------------
  int book_calc_histset[kNCalcHistSets];
  for (int i=0; i<kNCalcHistSets; i++) { book_calc_histset[i] = 0; }

  book_calc_histset[0] = 1;             // all
  book_calc_histset[1] = 1;             // all DISK 0
  book_calc_histset[2] = 1;             // all DISK 1

  int sel_disk[100];
  sel_disk[0] = -1;
  sel_disk[1] =  0;
  sel_disk[2] =  1;
  sel_disk[3] = -1;
  sel_disk[4] =  0;
  sel_disk[5] =  1;
  sel_disk[6] = -1;
  sel_disk[7] =  0;
  sel_disk[8] =  1;

  for (int i=0; i<kNCalcHistSets; i++) {
    if (book_calc_histset[i] == 0) continue;
    std::string folder_name = std::format("calc_{:02d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    Hist->calc[i] = new CalcHist_t();
    index.sel  = i;
    index.disk = sel_disk[i];
    BookCalcHistograms(Hist->calc[i],&index,fol);
  }

  return 0;
}

//-----------------------------------------------------------------------------
// need to optimize the filling time
//-----------------------------------------------------------------------------
int TCaloTimeAnaModule::FillCalcHistograms(CalcHist_t* Hist, TStnCluster* Calc, calc_param_t* Cp) {
  // filling histograms: plot time differences between
  
  Hist->h_edep->Fill(Calc->Energy());
  return 0;
}

//-----------------------------------------------------------------------------
// need to optimize the filling time
//-----------------------------------------------------------------------------
int TCaloTimeAnaModule::FillDiskHistograms(DiskHist_t* Hist, TCaloRecoDigi* Calrd) {
  // filling histograms: plot time differences between
  // Hist->h_ph->Fill(Crvp->ph);
  // Hist->h_npes->Fill(Crvp->npes);
  // Hist->h_time->Fill(Crvp->time);
  // Hist->h_feb->Fill(Crvp->feb);
  // Hist->h_ch->Fill(Crvp->ch);
  return 0;
}

//-----------------------------------------------------------------------------
// need to optimize the filling time
//-----------------------------------------------------------------------------
int TCaloTimeAnaModule::FillHistograms() {
  // filling histograms: plot time differences between

  //  Index_t index;
  
  for (int i=0; i<kNCrystals; ++i) {
    crystal_t* cr = &fCrystals[i];
    int n0 = cr->fCrd[0].size();
    int n1 = cr->fCrd[1].size();
    for (int i0=0; i0<n0; ++i0) {
      TCaloRecoDigi* crd0 = cr->fCrd[0].at(i0);
      // look for closest in time
      float min_dt(1.e6);
      for (int i1=0; i1<n1; ++i1) {
        TCaloRecoDigi* crd1 = cr->fCrd[1].at(i1);
        float dt = crd1->fTime-crd0->fTime;
        if (fabs(dt) < fabs(min_dt)) {
          min_dt = dt;
        }
      }
      // now histogram
      fHist->h_dt10_vs_crystal->Fill(cr->Cid(),min_dt);
    }
  }
  
  for (int i=0; i<fNCalh; i++) {
    TCaloHit*  calh = fCaloHitBlock->Hit(i);
    fHist->h_nsipms_vs_cid->Fill(calh->Cid(),calh->NSipms());
  }
//-----------------------------------------------------------------------------
// double-nested loops start here
//-----------------------------------------------------------------------------
  for (int i1=0; i1<fNCalord; i1++) {
    TCaloRecoDigi*  calrd = fCaloRecoDigiBlock->CaloRecoDigi(i1);

    fHist->h_sipmid->Fill(calrd->SipmID());
  }
//-----------------------------------------------------------------------------
// fill cluster histograms
//-----------------------------------------------------------------------------
//  int n_good_tc = 0;
  
  for (int i=0; i<fNCaloClusters; i++) {
    TStnCluster* calc        = fCaloClusterBlock->Cluster(i);
    calc_param_t*   calc_par = &fListOfCalcParam[i];
    FillCalcHistograms(fHist->calc[0],calc,calc_par);
    if (calc->DiskID() == 0) FillCalcHistograms(fHist->calc[1],calc,calc_par);
    else                     FillCalcHistograms(fHist->calc[2],calc,calc_par);
  }
  
  return 0;
}


//-----------------------------------------------------------------------------
int TCaloTimeAnaModule::CalculateMissingParameters() {

  fNCalh10[0] = 0;
  fNCalh10[1] = 0;

  fListOfCalcParam.clear();
  if (fNCaloClusters > 0) fListOfCalcParam.resize(fNCaloClusters);
//-----------------------------------------------------------------------------
// brute force clear - can this be done more intelligently ?
//-----------------------------------------------------------------------------
  for (int i=0; i<kNCrystals; i++) {
    crystal_t* cr = &fCrystals[i];
    cr->fCrd[0].clear();
    cr->fCrd[1].clear();
  }
//-----------------------------------------------------------------------------
// calo reco digis
//-----------------------------------------------------------------------------
  for (int i1=0; i1<fNCalord; i1++) {
    TCaloRecoDigi* crd = fCaloRecoDigiBlock->CaloRecoDigi(i1);
    int sipmid         = crd->SipmID();
    int cid            = sipmid / 2;
    int sipm           = sipmid % 2;

    crystal_t* cr      = &fCrystals[cid];
    cr->fCrd[sipm].push_back(crd);
  }
//-----------------------------------------------------------------------------
// calorimeter hits
// 1. clear hit channels
//-----------------------------------------------------------------------------
  int nhc = fHitCrystals.size();
  for (int i=0; i<nhc; ++i) {
    crystal_t* cr = fHitCrystals[i];
    cr->fHits.clear();
  }
  fHitCrystals.clear();
  fMaxHitsPerCrystal = 0;
  
  for (int i1=0; i1<fNCalh; i1++) {
    TCaloHit*  calh = fCaloHitBlock->Hit(i1);
    crystal_t* cr   = &fCrystals[calh->fCid];

    cr->fHits.push_back(calh);
    int nh = cr->NHits();
    if (nh == 1) {
      // first hit: add crystal to the list of crystals with hits
      fHitCrystals.push_back(cr);
    }
    
    if (nh > fMaxHitsPerCrystal) {
      fMaxHitsPerCrystal = nh;
    }
       
    if (calh->fEDep > 10) {
      int disk = calh->Disk();
      fNCalh10[disk] += 1;
    }
  }
//-----------------------------------------------------------------------------
// extra parameters of the calorimeter clusters
//-----------------------------------------------------------------------------
  fNCcDisk[0] = 0;
  fNCcDisk[1] = 0;
  
  for (int i1=0; i1<fNCaloClusters; i1++) {
    TStnCluster*  calc = fCaloClusterBlock->Cluster(i1);
    fNCcDisk[calc->DiskID()] += 1;

    calc_param_t& cp = fListOfCalcParam[i1];
//-----------------------------------------------------------------------------
// find the closest time cluster
//-----------------------------------------------------------------------------
    cp.dtmin_tc = 1.e6;
    cp.tc       = nullptr;
//-----------------------------------------------------------------------------
// find the closest CRV coincidence
//-----------------------------------------------------------------------------
    cp.dtmin_crvc = 1.e6;
    //    cp.crvc       = nullptr;
//------------------------------;-----------------------------------------------
// find the closest track
//-----------------------------------------------------------------------------
    cp.dtmin_trk = 1.e6;
    // cp.trk       = nullptr;
  }

  return 0;
}

//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TCaloTimeAnaModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
  RegisterDataBlock("CaloHitBlock"     ,"TCaloHitBlock"       ,&fCaloHitBlock     );
  RegisterDataBlock("CaloRecoDigiBlock","TCaloRecoDigiBlock"  ,&fCaloRecoDigiBlock);
  RegisterDataBlock("CaloClusterBlock" ,"TStnClusterBlock"    ,&fCaloClusterBlock );
 
  return 0;
}

//_____________________________________________________________________________
int TCaloTimeAnaModule::BeginRun() {
  fRunNumber = GetHeaderBlock()->RunNumber();
  TStntuple::Init(fRunNumber);
  //  fTpm->Init(fRunNumber);
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
int TCaloTimeAnaModule::Event(int IEntry) {

  //  TLorentzVector        mom;

  fCaloHitBlock->GetEntry(IEntry);
  fCaloRecoDigiBlock->GetEntry(IEntry);
  fCaloClusterBlock->GetEntry(IEntry);
//-----------------------------------------------------------------------------
// channel_# = 2*crystal_# + sipm_#
//-----------------------------------------------------------------------------
  fNCaloClusters = fCaloClusterBlock->NClusters();
  fNCalh         = fCaloHitBlock->NHits();
  fNCalord       = fCaloRecoDigiBlock->NDigis();

  CalculateMissingParameters();
  
  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TCaloTimeAnaModule::Debug() {

  if (GetDebugBit(3) == 1) {
    if ((fNCcDisk[0] > 0) and (fNCcDisk[1] > 0)) {
      GetHeaderBlock()->Print(Form("NCcDisk[0]:%2i fNCcDisk[1]:%2i",
                                   fNCcDisk[0],fNCcDisk[1]));
      fCaloClusterBlock->Print();
    }
  }

  if (GetDebugBit(5) == 1) {
    if (fMaxHitsPerCrystal > 1) {
      int n2plus = 0;
      int n_hit_crystals = fHitCrystals.size();
      for (int i=0; i<n_hit_crystals; i++) {
        crystal_t* cr = fHitCrystals.at(i);
        // figure the number of good hits in the crystal
        int nh = cr->NHits();
        int n_good_hits = 0;
        for (int ih=0; ih<nh; ih++) {
          TCaloHit* calh = cr->Hit(ih);
          if ((calh->fNSipms == 2) and (calh->fEDep > 10.)) {
            n_good_hits += 1;
          }
        }
        if (n_good_hits >= 2) {
          n2plus += 1;
          GetHeaderBlock()->Print(Form("CID:%4i n_good_hits:%2i",
                                       cr->Cid(),n_good_hits));
        }
      }
      if (n2plus > 0) {
          fCaloHitBlock->Print();
      }
    }
  }
}

//_____________________________________________________________________________
int TCaloTimeAnaModule::EndJob() {
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
int TCaloTimeAnaModule::FitCaloTimeOffsets(float TMin, float TMax) {

  TH2F* h2 = fHist->h_dt_vs_sipmid;

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

// //-----------------------------------------------------------------------------
// int TCaloTimeAnaModule::PrintTracks() {
//   std::cout << std::format(" i      T0         Z0       Nx     Ny    Nz    Chi2D     Xc[0]      Yc[0]      Xc[1]      Yc[1]    DxCal[0]    DyCal[0]    DxCal[1]   DyCal[1]\n");
//   for (int i=0; i<fNTrk; i++) {
//     TStrTrack*   trk = fTrackBlock->Track(i);
//     trk_param_t* tp  = &fListOfTrkParam.at(i);
//     std::cout << std::format("{:2d} {:10.3f} {:10.3f} {:6.3f} {:6.3f} {:6.3f} {:7.2f}",
//                              i,trk->fT0,trk->fZ0,trk->fNx,trk->fNy,trk->fNz,trk->fChi2/trk->fNDof);
    
//     std::cout << std::format(" {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f}\n",
//                              tp->x_disk [0],tp->y_disk [0],tp->x_disk[1],tp->y_disk[1],
//                              tp->dx_calc[0],tp->dy_calc[0],tp->dx_calc[1],tp->dy_calc[1]);
//   }
//   return 0;
// }
  
}
