//////////////////////////////////////////////////////////////////////////////
// use of tmp:
//
// use of debug bits: bits 0-2 are reserved
//  0  : all events
//  1  : passed events
//  2  : rejected events
//  3  : print additional parameters of the time cluster
//  4  : for max_n_hits_panel >= 6, print time cluster with hits
// 
///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <format>

#include "TROOT.h"
#include "TStyle.h"
#include "TPaveStats.h"
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
#include "murat/gui/TEvdDetNode.hh"
#include "murat/gui/TEvdManager.hh"

#include "Stntuple/loop/TStnInputModule.hh"
#include "murat/ana/TPanelDtAnaModule.hh"

using namespace murat;

ClassImp(murat::TPanelDtAnaModule)

namespace murat {
//-----------------------------------------------------------------------------
TPanelDtAnaModule::TPanelDtAnaModule(const char* name, const char* title):
  TStnModule(name,title)
{
  TH1::AddDirectory(0);
//-----------------------------------------------------------------------------
// detector numerology
//-----------------------------------------------------------------------------
  fTpm            = TrkPanelMap::Instance();
  
  //  fFrRef          = nullptr;
  fHist           = new Hist_t;
  
  fListOfTcParam.reserve(100);
  fListOfTcParam.resize (100);

  auto vm = TEvdManager::Instance();
  vm->InitGeometry(); // if not initialized

  fNEvents  = -1;
  fMaxEvent = -1;

  // fCanvasWf = nullptr;
  //   fEvdCrv->InitGeometry("Offline/Mu2eG4/geom/geom_common_extracted.txt");
}

//-----------------------------------------------------------------------------
TPanelDtAnaModule::~TPanelDtAnaModule() {
}

//-----------------------------------------------------------------------------
int TPanelDtAnaModule::BookChannelHistograms(ChannelHist_t* Hist, TrkIndex_t* Index, TFolder* Folder) {

  std::string prefix = std::format("run {:6d} slot:plane:panel:{:02d}:{:2d}:{:02d} MNID:{} channel:{:02d}",
                                   fRunNumber,Index->slot,Index->plane,Index->panel,Index->mnid,Index->ch);
  std::string name, title;

  // name  = "ph";
  // title = std::format("{} : pulse height",prefix);
  // fBookHist->HBook1F(Hist->h_ph,name.data(),title.data(),200,0,200,Folder);

  // name  = "bl";
  // title = std::format("{} : baseline",prefix);
  // fBookHist->HBook1F(Hist->h_bl,name.data(),title.data(),500,0,500,Folder);

  // name  = "fs";
  // title = std::format("{} : first sample",prefix);
  // fBookHist->HBook1F(Hist->h_fs,name.data(),title.data(),50,0,50,Folder);

  name  = "tdc0";
  title = std::format("{} : tdc0, us",prefix);
  fBookHist->HBook1F(Hist->h_tdc0,name.data(),title.data(),1500,0,150,Folder);

  name  = "dt10";
  title = std::format("{} : dt10=tdc1-tdc0, ns",prefix);
  fBookHist->HBook1F(Hist->h_dt10,name.data(),title.data(),1000,-50,50,Folder);

  name  = "edep";
  title = std::format("{} : edep, keV",prefix);
  fBookHist->HBook1F(Hist->h_edep,name.data(),title.data(),200,-0.001,0.009,Folder);

  return 0;
}

//-----------------------------------------------------------------------------
int TPanelDtAnaModule::BookEventHistograms(EventHist_t* Hist, TrkIndex_t* Index, TFolder* Folder) {

  std::string prefix = std::format("run:{:06d}",fRunNumber);
  std::string name, title;

  name  = "nch";
  title = std::format("{} : nch",prefix);
  fBookHist->HBook1F(Hist->h_nch,name.data(),title.data(),200,0,2000,Folder);

  name  = "nsh";
  title = std::format("{} : nsh",prefix);
  fBookHist->HBook1F(Hist->h_nsh,name.data(),title.data(),200,0,2000,Folder);

  name  = "ntc";
  title = std::format("{} : ntc",prefix);
  fBookHist->HBook1F(Hist->h_ntc,name.data(),title.data(),100,0,100,Folder);

  name  = "ntrk";
  title = std::format("{} : ntrk",prefix);
  fBookHist->HBook1F(Hist->h_ntrk,name.data(),title.data(),10,0,10,Folder);

  return 0;
}

//-----------------------------------------------------------------------------
int TPanelDtAnaModule::BookPanelHistograms(PanelHist_t* Hist, TrkIndex_t* Index, TFolder* Folder) {

  std::string prefix = std::format("run:{:06d} slot:{:02d} plane:{:02d} panel:{:03d} mnid:MN{:03d}",
                                   fRunNumber,Index->slot,Index->plane,Index->panel,Index->mnid);
  std::string name, title;

  name  = "occup";
  title = std::format("{} : occupancy",prefix);
  fBookHist->HBook1F(Hist->h_occup,name.data(),title.data(),100,0,100,Folder);

  name  = "edep";
  title = std::format("{} : edep, keV",prefix);
  fBookHist->HBook1F(Hist->h_edep,name.data(),title.data(),200,-0.001,0.009,Folder);

  for (int i=0; i<96; i++) {
    std::string folder_name = std::format("ch_{:02d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    Hist->ch[i] = new ChannelHist_t;
    Index->ch   = i;
    BookChannelHistograms(Hist->ch[i],Index,fol);
  }
  return 0;
}


//-----------------------------------------------------------------------------
int TPanelDtAnaModule::BookSlotHistograms(SlotHist_t* Hist, TrkIndex_t* Index, TFolder* Folder) {

  std::string prefix = std::format("run:{:06d} slot:{:02d}",fRunNumber,Index->slot);
  std::string name, title;

  name  = "edep";
  title = std::format("{} : edep",prefix);
  fBookHist->HBook1F(Hist->h_edep,name.data(),title.data(),200,-0.001,0.009,Folder);

  for (int i=0; i<12; i++) {
    std::string folder_name = std::format("pnl_{:03d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    Hist->panel[i] = new PanelHist_t;
    
    Index->plane = 2*Index->slot + i/6;
    Index->panel = i % 6;
    Index->pnl12 = i;
    //    TLOG(TLVL_DEBUG+1) << std::format("Index->plane:{:02d} Index->panel:{}",Index->plane,Index->panel);

    TrkPanelMap::Data_t* tpmd = fTpm->panel_data_by_offline(Index->plane,Index->panel);
    Index->mnid    = tpmd->mnid;
    
    BookPanelHistograms(Hist->panel[i],Index,fol);
  }
  return 0;
}


//-----------------------------------------------------------------------------
int TPanelDtAnaModule::BookSelHistograms(SelHist_t* Hist, TrkIndex_t* Index, TFolder* Folder) {

  std::string prefix = std::format("run:{:06d} selection:{:02d}",fRunNumber,Index->sel);
  std::string name, title;

  name  = "occup_2d";
  title = std::format("{} : occupancy, panel vs ch",prefix);
  fBookHist->HBook2F(Hist->h_occup_2d,name.data(),title.data(),100,0,100,216,0,216,Folder);

  name  = "occup";
  title = std::format("{} : per panel occupancy",prefix);
  fBookHist->HBook1F(Hist->h_occup,name.data(),title.data(),216,0,216,Folder);

  TrkIndex_t index;

  // for (int i=0; i<36; i++) {
  //   name  = std::format("dt05_%{}",i);
  //   title = std::format("{} : dt05 ({:02} - {:02})",prefix,i,i-1);
  //   fBookHist->HBook1F(Hist->h_dt05[i],name.data(),title.data(),1000,-5000,5000,Folder);
  // }
  
  for (int i=0; i<18; i++) {
    std::string folder_name = std::format("slot_{:02d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    Hist->slot[i] = new SlotHist_t;
    index.slot  =  i;
    index.panel = -1;
    index.mnid  = -1;
    index.ch    = -1;
    BookSlotHistograms(Hist->slot[i],&index,fol);
  }
  return 0;
}

//-----------------------------------------------------------------------------
int TPanelDtAnaModule::BookTcHistograms(TcHist_t* Hist, TrkIndex_t* Index, TFolder* Folder) {

  std::string prefix = std::format("run:{:06d} selection:{:02d}",fRunNumber,Index->sel);
  std::string name, title;

  name  = "nch";
  title = std::format("{} : N(combo hits)",prefix);
  fBookHist->HBook1F(Hist->h_nch,name.data(),title.data(),200,0,200,Folder);

  name  = "ngh";
  title = std::format("{} : N(good hits)",prefix);
  fBookHist->HBook1F(Hist->h_ngh,name.data(),title.data(),200,0,200,Folder);

  name  = "nplanes";
  title = std::format("{} : N(planes)",prefix);
  fBookHist->HBook1F(Hist->h_nplanes,name.data(),title.data(),36,0,36,Folder);

  name  = "npanels";
  title = std::format("{} : N(panels)",prefix);
  fBookHist->HBook1F(Hist->h_npanels,name.data(),title.data(),250,0,250,Folder);

  name  = "max_nh_plane";
  title = std::format("{} : Max NH per plane)",prefix);
  fBookHist->HBook1F(Hist->h_max_nh_plane,name.data(),title.data(),250,0,250,Folder);

  name  = "max_nh_panel";
  title = std::format("{} : Max NH per panel)",prefix);
  fBookHist->HBook1F(Hist->h_max_nh_panel,name.data(),title.data(),100,0,100,Folder);

  name  = "edep_vs_pnl";
  title = std::format("{} : EDep vs panel",prefix);
  fBookHist->HBook2F(Hist->h_edep_vs_pnl,name.data(),title.data(),216,0,216,200,-0.1e-3,9.9e-3,Folder);

  name  = "dt_vs_pnl";
  title = std::format("{} : T(panel)-T0(TC) vs panel, ns",prefix);
  fBookHist->HBook2F(Hist->h_dt_vs_pnl,name.data(),title.data(),216,0,216,400,-100,100,Folder);

  for (int i1=0; i1<36; i1++) {
    for (int i2=i1+1; i2<36; i2++) {
      name  = std::format("dt_pln_{:02}{:02}",i2,i1);
      title = std::format("{}: time cluster  T[{:02}]-T[{:02}]",prefix,i2,i1);
      fBookHist->HBook1F(Hist->h_dt_pln[i2][i1],name.data(),title.data(),1000,-500,500,Folder);
    }
  }

  return 0;
}

//-----------------------------------------------------------------------------
  int TPanelDtAnaModule::BookHistograms(Hist_t* Hist, TFolder* Folder) {

  std::string prefix = std::format("");
  std::string name, title;

  TrkIndex_t index;
//-----------------------------------------------------------------------------
// events
//-----------------------------------------------------------------------------
  int book_event_histset[kNEventHistSets];
  for (int i=0; i<kNEventHistSets; i++) { book_event_histset[i] = 0; }

  book_event_histset[0] = 1;                  // all events

  for (int i=0; i<kNEventHistSets; i++) {
    if (book_event_histset[i] == 0) continue;
    std::string folder_name = std::format("evt_{:02d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    Hist->event[i] = new EventHist_t;
    index.sel   =  i;
    BookEventHistograms(Hist->event[i],&index,fol);
  }
//-----------------------------------------------------------------------------
// straw hits , two selections
//-----------------------------------------------------------------------------
  int book_sel_histset[kNSelHistSets];
  for (int i=0; i<kNSelHistSets; i++) { book_sel_histset[i] = 0; }

  book_sel_histset[0] = 1;                  // all events
  book_sel_histset[1] = 1;                  // 

  for (int i=0; i<kNSelHistSets; i++) {
    if (book_sel_histset[i] == 0) continue;
    std::string folder_name = std::format("sel_{:02d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    Hist->sel[i] = new SelHist_t;
    index.sel   =  i;
    index.slot  = -1;
    index.panel = -1;
    index.mnid  = -1;
    index.ch    = -1;
    BookSelHistograms(Hist->sel[i],&index,fol);
  }
//-----------------------------------------------------------------------------
// time clusters
//-----------------------------------------------------------------------------
  int book_tc_histset[kNTcHistSets];
  for (int i=0; i<kNTcHistSets; i++) { book_tc_histset[i] = 0; }

  book_tc_histset[0] = 1;                  // all 
  book_tc_histset[1] = 1;                  // N(hits) >= 8 , ngh >= 2

  for (int i=0; i<kNTcHistSets; i++) {
    if (book_tc_histset[i] == 0) continue;
    std::string folder_name = std::format("tc_{:02d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    Hist->tc[i]  = new TcHist_t;
    index.sel   =  i;
    BookTcHistograms(Hist->tc[i],&index,fol);
  }
  return 0;
}



//-----------------------------------------------------------------------------
int TPanelDtAnaModule::FillChannelHistograms(ChannelHist_t* Hist, TrkIndex_t* Index, TStrawHit* Sh) {
  // float tdc0_ns = Sd->tdc0*5./256;
  // float tdc1_ns = Sd->tdc1*5./256;
  float dt10    = Sh->Time(1)-Sh->Time(0);
  
  Hist->h_tdc0->Fill(Sh->Time(0));
  Hist->h_dt10->Fill(dt10);

  Hist->h_edep->Fill(Sh->EDep());
                                        // no edep for the moment
  return 0;
}

//-----------------------------------------------------------------------------
int TPanelDtAnaModule::FillEventHistograms(EventHist_t* Hist, TrkIndex_t* Index) {
  Hist->h_nsh->Fill(fNSh);
  Hist->h_nch->Fill(fNCh);
  Hist->h_ntc->Fill(fNTc);
  Hist->h_ntrk->Fill(fNTracks);
  return 0;
}

//-----------------------------------------------------------------------------
int TPanelDtAnaModule::FillPanelHistograms(PanelHist_t* Hist, TrkIndex_t* Index, TStrawHit* Sh) {
  Hist->h_occup->Fill(Index->ch);
  Hist->h_edep->Fill(Sh->EDep());
  return 0;
}

//-----------------------------------------------------------------------------
int TPanelDtAnaModule::FillSlotHistograms(SlotHist_t* Hist, TrkIndex_t* Index, TStrawHit* Sh) {
  
  Hist->h_edep->Fill(Sh->EDep());
  return 0;
}

//-----------------------------------------------------------------------------
int TPanelDtAnaModule::FillSelHistograms(SelHist_t* Hist, TrkIndex_t* Index, TStrawHit* Sh) {
  
  int offline_panel = Index->slot*12+Index->pnl12;
  Hist->h_occup->Fill(offline_panel);
  Hist->h_occup_2d->Fill(Index->ch,offline_panel);
  
  SlotHist_t* slot_hist_set = Hist->slot[Index->slot];
  FillSlotHistograms(slot_hist_set,Index,Sh);
  
  PanelHist_t* panel_hist_set = slot_hist_set->panel[Index->pnl12];
  FillPanelHistograms(panel_hist_set,Index,Sh);

  ChannelHist_t* ch_hist_set = panel_hist_set->ch[Index->ch];
  FillChannelHistograms(ch_hist_set,Index,Sh);

  return 0;
}

//-----------------------------------------------------------------------------
  int TPanelDtAnaModule::FillTcHistograms(TcHist_t*        Hist ,
                                          TrkIndex_t*      Index,
                                          TStnTimeCluster* Tc   ,
                                          tc_param_t*      Tcp  ) {

    Hist->h_nch->Fill(Tc->NComboHits());
    Hist->h_ngh->Fill(Tcp->ngh);
    Hist->h_nplanes->Fill(Tcp->nplanes);
    Hist->h_npanels->Fill(Tcp->npanels);
    Hist->h_max_nh_panel->Fill(Tcp->max_nh_panel);
    Hist->h_max_nh_plane->Fill(Tcp->max_nh_plane);
    
    for (int pln1=0; pln1<36; pln1++) {
      if (Tcp->nh_plane[pln1] < 2) continue;
      for (int pln2=pln1+1; pln2<36; pln2++) {
        if (Tcp->nh_plane[pln2] < 2) continue;
        // two planes with 2 or more hits in this time cluster
        float dt = Tcp->time_plane[pln2]-Tcp->time_plane[pln1];
        Hist->h_dt_pln[pln2][pln1]->Fill(dt);
      }
    }
    // mean EDep
    for (int p216=0; p216<216; p216++) {
      if (Tcp->nh_panel[p216] < 2) continue;
      Hist->h_edep_vs_pnl->Fill(p216,Tcp->edep_panel[p216]);
    }

    if (Tcp->npanels > 10) {
      for (int p216=0; p216<216; p216++) {
        if (Tcp->nh_panel[p216] < 2) continue;
        float dt = Tcp->time_panel[p216]-Tc->T0();
        Hist->h_dt_vs_pnl->Fill(p216,dt);
      }
    }

    return 0;
  }

//-----------------------------------------------------------------------------
// all histograms here are for straw hits only,
// what about time clusters ?
//-----------------------------------------------------------------------------
int TPanelDtAnaModule::FillHistograms() {
  TrkIndex_t index;
//-----------------------------------------------------------------------------
// event
//-----------------------------------------------------------------------------
  FillEventHistograms(fHist->event[0],&index);
//-----------------------------------------------------------------------------
// straw hit selections
//-----------------------------------------------------------------------------
  for (int i=0; i<fNSh; i++) {
    TStrawHit*  sh = fShBlock->Hit(i);
    // std::cout << std::format("i:{:5d} sd_mnid[i]:{:03d}\n",i,sd->mnid);
    mu2e::StrawId sid = sh->Sid();
    index.plane = sid.plane();
    index.panel = sid.panel();
    index.ch    = sid.straw();
    index.slot  = index.plane/2;
    index.pnl12 = 6*(index.plane%2)+index.panel;
    index.mnid  = sh->Mnid();

    index.sel   = 0;
    FillSelHistograms(fHist->sel[0],&index,sh);
    
    if (sh->EDep() > 0.0005) {
      index.sel   = 1;
      FillSelHistograms(fHist->sel[1],&index,sh);
    }
  }
//-----------------------------------------------------------------------------
// time clusters
//-----------------------------------------------------------------------------
  for (int i=0; i<fNTc; i++) {
    TStnTimeCluster* tc  = fTcBlock->TimeCluster(i);
    tc_param_t*      tcp = &fListOfTcParam[i];
    index.sel = 0;
    FillTcHistograms(fHist->tc[0],&index,tc,tcp);
    if ((tc->NComboHits() >= 8) and (tcp->ngh >= 2)) {
      index.sel = 1;
      FillTcHistograms(fHist->tc[1],&index,tc,tcp);
    }
  }
  return 0;
}


 
//-----------------------------------------------------------------------------
int TPanelDtAnaModule::CalculateMissingParameters() {

  for (int itc=0; itc<fNTc; itc++) {
    TStnTimeCluster* tc = fTcBlock->TimeCluster(itc);
    tc_param_t* tcp     = &fListOfTcParam[itc];
    // initialize tcp
    tcp->npanels      = 0;
    tcp->nplanes      = 0;
    tcp->ngh          = 0;
    tcp->max_nh_plane = 0;
    tcp->max_nh_panel = 0;
    tcp->max_nh_face  = 0;
    tcp->max_edep     = 0.;
    tcp->tmin         = 0;
    tcp->tmax         = 0;
    
    // can do 1D, by panel
    for (int p216=0; p216<216; p216++) {
      tcp->nh_panel  [p216] = 0;
      tcp->ngh_panel [p216] = 0;
      tcp->time_panel[p216] = 0;
      tcp->edep_panel[p216] = 0;
    }
    
      // can do 1D, just by plane 
    for (int pln=0; pln<36; pln++) {
      tcp->nh_plane  [pln] = 0;
      tcp->ngh_plane [pln] = 0;
      tcp->time_plane[pln] = 0;
    }

    for (int stn=0; stn<18; stn++) {
      for (int face=0; face<4; face++) {
        tcp->nh_face[stn][face] = 0;
        tcp->np_face[stn][face] = 0;
      }
    }
    
    TStnLinkBlock* ch_links = fTcBlock->ListOfChLinks();

    int nch = tc->NComboHits();
    for (int ih=0; ih<nch; ih++) {
      // index of the combohit in the combohit block
      int hit_id = ch_links->Index(itc,ih);
      
      TComboHit* ch = fChBlock->GetHitByID(hit_id);
      if (ch == nullptr) {
        GetHeaderBlock()->Print(Form("ERROR: CH id=%i is missing. SKIP",hit_id));
        continue;
      }

      int plane = ch->Sid().plane();
      int stn   = plane / 2;
      // int pln2  = plane % 2;            // 0 or 1
      int panel = ch->Sid().panel();
      const TrkPanelMap::Data_t* tpmd = fTpm->panel_data_by_offline(plane,panel);

      int zface = tpmd->zface;          // ordered in Z
      //int p12   = 6*pln2+panel;
      int p216  = plane*6+panel;        // offline plane index
      
      if (tcp->nh_panel[p216] == 0) tcp->npanels++;
      tcp->nh_panel    [p216] += 1;
      tcp->time_panel  [p216] += ch->Time(); // this one is corrected
      tcp->edep_panel  [p216] += ch->EDep();

      if (ch->EDep() > tcp->max_edep) tcp->max_edep = ch->EDep();

      if (ch->EDep() > 0.0005) {
        tcp->ngh              += 1;
        tcp->ngh_panel[p216 ] += 1;
        tcp->ngh_plane[plane] += 1;
      }

      if (tcp->nh_plane[plane] == 0) tcp->nplanes++;
      tcp->nh_plane    [plane] += 1;
      tcp->time_plane  [plane] += ch->Time();

      if (tcp->nh_face[stn][zface] == 0) tcp->nfaces++;
      tcp->nh_face  [stn][zface] += 1;
      tcp->time_face[stn][zface] += ch->Time();

      if (ch->Time() < tcp->tmin) tcp->tmin = ch->Time();
      if (ch->Time() > tcp->tmax) tcp->tmax = ch->Time();
    }
//-----------------------------------------------------------------------------
// average time of the hits in a given plane
//-----------------------------------------------------------------------------
    for (int pln=0; pln<36; pln++) {
      if (tcp->nh_plane[pln] > 1) {
        tcp->time_plane[pln] = tcp->time_plane[pln]/tcp->nh_plane[pln];
      }
      if (tcp->nh_plane[pln] > tcp->max_nh_plane) {
        tcp->max_nh_plane = tcp->nh_plane[pln];
      }
    }
    
    for (int stn=0; stn<18; stn++) {
//-----------------------------------------------------------------------------
// average time of the hits in a given face
//-----------------------------------------------------------------------------
      for (int face=0; face<4; face++) {
        if (tcp->nh_face[stn][face] > 1) {
          tcp->time_face[stn][face] = tcp->time_face[stn][face]/tcp->nh_face[stn][face];
        }
        if (tcp->nh_face[stn][face] > tcp->max_nh_face) {
          tcp->max_nh_face = tcp->nh_face[stn][face];
        }
      }
    }
//-----------------------------------------------------------------------------
// average time and charge of the hits in a given panel
//-----------------------------------------------------------------------------
    for (int p216=0; p216<216; p216++) {
        // do not need division if p12=0,1
      if (tcp->nh_panel[p216] > 1) {
        tcp->time_panel[p216] = tcp->time_panel[p216]/tcp->nh_panel[p216];
        tcp->edep_panel[p216] = tcp->edep_panel[p216]/tcp->nh_panel[p216];
      }

      if (tcp->nh_panel[p216] > tcp->max_nh_panel) {
        tcp->max_nh_panel = tcp->nh_panel[p216];
      }
    }
  }
  return 0;
}

//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TPanelDtAnaModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
  RegisterDataBlock("TimeClusterBlock" ,"TStnTimeClusterBlock",&fTcBlock          );
  RegisterDataBlock("ComboHitBlock"    ,"TComboHitBlock"      ,&fChBlock          );
  RegisterDataBlock("StrawHitBlock"    ,"TStrawHitBlock"      ,&fShBlock          );
  RegisterDataBlock("TrackBlock"       ,"TStrTrackBlock"      ,&fTrackBlock       );

  fNEvents = GetAna()->GetInputModule()->GetEntries();
  return 0;
}

//_____________________________________________________________________________
int TPanelDtAnaModule::BeginRun() {
  fRunNumber = GetHeaderBlock()->RunNumber();
  TStntuple::Init(fRunNumber);
  fTpm->Init(fRunNumber);
//-----------------------------------------------------------------------------
// book histograms at begin run, there will be a flag telling whether we need
// to re-book them for every new run
// don't do that by default
//-----------------------------------------------------------------------------
  // relative positions of subdetectors wrt each other
  //  int print_level(0);
  //  fDetPos.GetFromRunDb(fRunNumber,print_level);
  
  BookHistograms(fHist,fFolder);
  return 0;
}

//-----------------------------------------------------------------------------
// 2026-09-09
//-----------------------------------------------------------------------------
int TPanelDtAnaModule::Event(int IEntry) {

  fChBlock->GetEntry(IEntry);
  fShBlock->GetEntry(IEntry);
  fTcBlock->GetEntry(IEntry);
  fTrackBlock->GetEntry(IEntry);
//-----------------------------------------------------------------------------
// expect ChBlock and ShBlock to have only "close" hits, i.e. using compressed
// ntuples
//-----------------------------------------------------------------------------
  fNTc           = fTcBlock->NTimeClusters();
  fNTracks       = fTrackBlock->NTracks();
  fNCh           = fChBlock->NHits();
  fNSh           = fShBlock->NHits();

  int evn = GetHeaderBlock()->EventNumber();
  if (evn > fMaxEvent) {
    fMaxEvent = evn;
  }

  CalculateMissingParameters();
  
  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TPanelDtAnaModule::PrintComboHits() {
  for (int i=0; i<fNTc; i++) {
    TStnTimeCluster* tc  = fTcBlock->TimeCluster(i);
    tc_param_t*      tcp = &fListOfTcParam[i];
      
    bool banner_printed(false);
    if (tcp->max_nh_panel >= 6) {
        
      if (not banner_printed) {
        tc->Print("banner");
      }
      tc->Print("data");

      int nch = tc->NComboHits();

      bool ch_banner_printed(false);
      for (int ich=0; ich<nch; ich++) {
        int hit_id = fTcBlock->ListOfChLinks()->Index(i,ich);
        TComboHit* ch = fChBlock->GetHitByID(hit_id);
        if (ch == nullptr) {
          std::cout << std::format("ERROR: Time cluster CH id={:4} is missing. SKIP\n",hit_id);
          continue;
        }
        if (not ch_banner_printed) {
          ch->Print("banner");
          ch_banner_printed = true;
        }
        ch->Print("data");
      }
    }
  }
}

//-----------------------------------------------------------------------------
void TPanelDtAnaModule::Debug() {

  if (GetDebugBit(0) == 1) {
    GetHeaderBlock()->Print(Form("bit_00: all events"));
  }

  if (GetDebugBit(3) == 1) {
    // additional time cluster parameters 
    GetHeaderBlock()->Print("bit_003:");
    for (int i=0; i<fNTc; i++) {
      TStnTimeCluster* tc  = fTcBlock->TimeCluster(i);
      tc_param_t*      tcp = &fListOfTcParam[i];
      
      printf("%5i %5i %12.2f %12.2f %12.2f\n",
             tc->NStrawHits(),
             tc->NComboHits(),
             tc->T0(),
             tcp->tmin,
             tcp->tmax);
    }
  }

  if (GetDebugBit(4) == 1) {
    // additional time cluster parameters 
    GetHeaderBlock()->Print("bit_004:");
    PrintComboHits();
  }
}

//_____________________________________________________________________________
int TPanelDtAnaModule::EndJob() {
  return 0;
}


//-----------------------------------------------------------------------------
  int TPanelDtAnaModule::FitPlaneTimeOffsets(float TMin, float TMax, int DebugLevel) {
  
  TStnModule::fit_result_t frr[35];

  frr[0].p[1]    = 0;
  frr[0].e[1]    = 0;
  frr[0].p[2]    = 0;
  frr[0].e[2]    = 0;
  frr[0].chi2dof = -1.;
  
  for (int i=1; i<36; i++) {
    int i2=i-1;
    // run 124155: station 7 (planes 14 and 15) was off
    // if (i == 16) i2 = 13;
    TH1F* h = fHist->tc[0]->h_dt_pln[i][i2];
    
    fit_result_t* fr = &frr[i];
    fr->chi2dof = -1;
    FitHistogram(h,fr,10, TMin,TMax,DebugLevel);
  }
  
  // done fitting, print results

  // std::ofstream os("CaloTimeCalib_corr.txt");

  // os << std::format("TABLE CalTimeCalib {}\n",fRunNumber);
  // os << std::format("# corrections to  CalTimeCalib cid= XXX, to be subtracted \n");
  // os << std::format("# sipmid     dT      sigma      chi2 \n");
  
  std::cout << std::format("|   i  | correction |   Mean |  sig(mean) |  width | sig(width) |  chi2d |  \n");

  float vsum = 0;
  for (int i=0; i<36; i++) {
    vsum += frr[i].p[1];
    for (int j=0; j<6; j++) {
      int pnl = 6*i+j;
      // for plane-level corrections, the correction is the same for all panels in a given plane
      std::cout << std::format("|{:4} | {:10.3f} | {:10.3f} | {:10.3f} | {:10.3} | {:10.3f} |  {:10.3f} |\n",
                               pnl,vsum,frr[i].p[1],frr[i].e[1],frr[i].p[2],frr[i].e[2],frr[i].chi2dof);
    }
  }
  // plot histogram
  TH1F* h_m = new TH1F("h_m","Plane #Dt fit summary: mean values",36,0,36);
  for (int i=0; i<36; i++) {
    h_m->SetBinContent(i+1,frr[i].p[1]);
    h_m->SetBinError  (i+1,frr[i].e[1]);
  }

  TH1F* h_w = new TH1F("h_w","Plane #Dt fit summary: widths",36,0,36);
  for (int i=0; i<36; i++) {
    h_w->SetBinContent(i+1,frr[i].p[2]);
    h_w->SetBinError  (i+1,frr[i].e[2]);
  }

  TCanvas* c2 = new TCanvas("c2","c2",1200,600);
  c2->Divide(2,1);
  c2->Draw();
  c2->cd(1);
  h_m->Draw();
  c2->cd(2);
  h_w->Draw();
  
  return 0;
}

//-----------------------------------------------------------------------------
// one occupancy canvas per station
//-----------------------------------------------------------------------------
void TPanelDtAnaModule::FitPanelTimeOffsets(TH2F* H_DtVsPanel, float Dt, int PrintLevel) {

  // gROOT->SetBatch(kTRUE);   // no GUI windows

  // std::string fn = std::format("run_{:6d}_n002_tc.pdf",fRunNumber);

  TH2F* h2 =  H_DtVsPanel;
  if (h2 == nullptr) {
    h2 = fHist->tc[0]->h_dt_vs_pnl;
  }
  
  gStyle->SetStatW(0.30);   // wider (NDC)

  std::string prefix = std::format("run {:6d}",fRunNumber);
  std::string name, title;

  name  = "dt_vs_panel_fit";
  title = std::format("{} : dt vs panel fit",prefix);
  fHist->h_dt_vs_panel_fit = new TH1F(name.data(),title.data(),216,0,216);
  
  name  = "pdt";
  title = std::format("{} : panel Dt",prefix);
  fHist->h_pdt             = new TH1F(name.data(),title.data(),500,-50,50);
  
  for (int ip=0; ip<216; ip++) {
    int binx = ip+1;
    TH1* hpy = h2->ProjectionY("py",binx,binx);
    FitHistogram(hpy,&fPanelFr[ip],ip,0);
    fit_result_t* fr = &fPanelFr[ip];
    // finx maximum, fit within +/- Dt from it
    FitHistogram(hpy,fr,Dt, 1,-1,PrintLevel);

    float y  = fPanelFr[ip].p[1];
    float ey = fPanelFr[ip].e[1];
    fHist->h_dt_vs_panel_fit->SetBinContent(ip+1,y );
    fHist->h_dt_vs_panel_fit->SetBinError  (ip+1,ey);
    if (fPanelFr[ip].chi2dof > 0) {
      fHist->h_pdt->Fill(y);
    }
  }

  fHist->h_dt_vs_panel_fit->Draw();

  // finally, print the resulting table
  // not sure yet how to name the file

  std::cout << std::format("# i    correction      err      fit_chi2\n");
  for (int ip=0; ip<216; ip++) {
    //    corr[ip] = -1;
    fit_result_t* fr = &fPanelFr[ip];

    // inverse signs of the step and total offsets (to add to calibrations)
    std::cout  << std::format(" {:3d} {:10.4f} {:10.4f} {:10.4f} \n",ip,fr->p[1],fr->e[1],fr->chi2dof);
  }
  
  return;
}


// //-----------------------------------------------------------------------------
// int TPanelDtAnaModule::PrintTracks() {
//   std::cout << std::format(" i      T0         Z0       Nx     Ny    Nz    Chi2D     Xc[0]      Yc[0]      Xc[1]      Yc[1]    DxCal[0]    DyCal[0]    DxCal[1]   DyCal[1]\n");
//   for (int i=0; i<fNTracks; i++) {
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


//-----------------------------------------------------------------------------
  int TPanelDtAnaModule::InitEvd() {
    int rc(0);

    auto vm = TEvdManager::Instance();
    
    vm->AddDataBlock("HeaderBlock"      ,GetHeaderBlock());
    vm->AddDataBlock("TrackBlock"       ,fTrackBlock);
//-----------------------------------------------------------------------------
// 3D detector view
//-----------------------------------------------------------------------------
    TEvdDetNode* det_node = new TEvdDetNode("det_node");
    TEvdView*    det_view = new TEvdView   ("det_view",TGLViewer::kCameraPerspXOY);
    det_view->AddNode(det_node);
    vm->AddView(det_view);

    return rc;
  }
  
//-----------------------------------------------------------------------------
  int TPanelDtAnaModule::DisplayEvent() {
    int rc(0);
    auto vm = TEvdManager::Instance();

    vm->DisplayEvent();

    return rc;
  }
  
//-----------------------------------------------------------------------------
// one occupancy canvas per station
//-----------------------------------------------------------------------------
int TPanelDtAnaModule::PrintHistograms(int ISet) {

  gROOT->SetBatch(kTRUE);   // no GUI windows

  std::string fn = std::format("run_{:6d}_set_{:02}_occup.pdf",fRunNumber,ISet);
  
  gStyle->SetStatW(0.30);   // wider (NDC)

 for (int is=0; is<18; is++) {
    TCanvas c(Form("c_%02i",is),Form("c_%02i",is),1600,1800);
    c.Divide(3,4);
    for (int ip=0; ip<12; ip++) {
      TH1F* h = fHist->sel[ISet]->slot[is]->panel[ip]->h_occup;
      c.cd(ip+1);
      float hmax = (int(fMaxEvent/1.e6)+1)*1e6;
      // normalization to the rate :

      float input_rate = 1.e4;   // 10 kHz 
      float scale = input_rate/(fNEvents+1.e-12);
      h->Scale(scale);
      h->SetMaximum(hmax);
      gPad->SetLogy(kTRUE);
      h->Draw("hist");
      // make statbox transparent
      gPad->Update();           // create stats box

      auto st = (TPaveStats*)h->FindObject("stats");
      if (st) {
        st->SetFillStyle(0);    // transparent
        st->SetBorderSize(1);   // optional
      }
      gPad->Modified();
      gPad->Update();
    }

    if (is == 0) {
      c.Print(Form("%s(",fn.data()));     // or .pdf, .root, ...
    }
    else if (is == 17) {
      c.Print(Form("%s)",fn.data()));     // or .pdf, .root, ...
    }
    else {
      c.Print(fn.data());     // or .pdf, .root, ...
    }
  }

  return 0;
}

//-----------------------------------------------------------------------------
// print channels which contain more than Percentage of the total, default - 1%
//-----------------------------------------------------------------------------
int TPanelDtAnaModule::PrintNoisyChannels(float Percentage) {

  int   nbx  = fHist->sel[1]->h_occup_2d->GetNbinsX();
  int   nby  = fHist->sel[1]->h_occup_2d->GetNbinsY();
  float qtot = fHist->sel[1]->h_occup_2d->GetEntries();

  std::cout << std::format("qtot:{}\n",qtot);

  std::cout << std::format("   ix    iy  plane  panel   channel   DTC   link   MNID   nxy   percentage\n");
  for (int ix=0; ix<nbx; ix++) {
    for (int iy=0; iy<nby; iy++) {
      float nxy = fHist->sel[1]->h_occup_2d->GetBinContent(ix+1,iy+1);
      if (nxy > qtot*Percentage) {
        // ix is the channel number
        int plane = iy / 6;
        int panel = iy % 6;
        TrkPanelMap::Data_t* tpmd = fTpm->panel_data_by_offline(plane,panel);
        
        std::cout << std::format("{:4} {:4}   {:02d}   {}    {:02d} {}  {} MN{:03d} {:8} {:8.4f}\n",
                                 ix,iy,plane,panel,ix,tpmd->dtc_id,tpmd->link,tpmd->mnid,
                                 nxy,float(nxy)/qtot);
      }
    }
  }

  return 0;
}
  

}
