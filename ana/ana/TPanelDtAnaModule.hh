///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef __murat_ana_TPanelDtAnaModule_hh__
#define __murat_ana_TPanelDtAnaModule_hh__

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TCanvas.h"

#include "Stntuple/loop/TStnModule.hh"

#include "Stntuple/obj/TStrawHitBlock.hh"
#include "Stntuple/obj/TComboHitBlock.hh"
#include "Stntuple/obj/TStnTimeClusterBlock.hh"
#include "Stntuple/obj/TStrTrackBlock.hh"

#include "Stntuple/geom/TrkPanelMap.hh"

// subdetector-level alignment table
#include "murat/alg/DetPos_t.hh"

namespace murat {
class TPanelDtAnaModule: public TStnModule {
  
public:
  
  enum {
    kNStations         = 18,
    kNPlanes           = 36,
    kNPanelsPerStation = 12,
  };


  enum {
    kNHistSets         = 100,
    kNSelHistSets      = 100,
    kNTcHistSets       = 100,
  };
//-----------------------------------------------------------------------------
// indices
//-----------------------------------------------------------------------------
  struct TrkIndex_t {
    int sel;
    int slot;                 // 0-17
    int plane;                // offline
    int panel;                // offline
    int pnl12;                // panel index within the station (0-11)
    int mnid;
    int ch;
  };
//-----------------------------------------------------------------------------
// histograms
//-----------------------------------------------------------------------------
  struct ChannelHist_t {
    TH1F*     h_tdc0;
    TH1F*     h_dt10;                   // HV-CAL
    TH1F*     h_edep;                   // 
  };
  
  struct PanelHist_t {
    ChannelHist_t* ch[96];
    TH1F*          h_occup;             // pulse height
    TH1F*          h_edep;              // 
  } ;
  
  struct SlotHist_t {
    PanelHist_t* panel[12];
    TH1F*        h_edep;                // 
  };

  struct SelHist_t {
    TH1F*        h_occup;               // occupancy per panel, offline indexing
    TH2F*        h_occup_2d;            // occupancy
    TH1F*        h_dt05[36];
    SlotHist_t*  slot[18];
  };

  struct TcHist_t {
    TH1F*        h_nch;                 // same as nsh
    TH2F*        h_edep_vs_pnl;
    TH1F*        h_dt_pln[36][36];
    TH2F*        h_dt_vs_pnl;             // offline index : plane*6+pln6
  };
  
  struct Hist_t {
    SelHist_t* sel[kNSelHistSets];      // select straw hits
    TcHist_t*  tc [kNTcHistSets ];
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
  struct tc_param_t {
    int   npanels;
    int   ngh;                          // N "good" hits (EDep > 0.5 keV)
    int   nh_panel  [216];
    int   ngh_panel [216];
    float time_panel[216];
    float edep_panel[216];
    int   max_nh_panel;
    float max_edep;
    int   nplanes;
    int   nh_plane  [36];                 // N hits per plane
    int   ngh_plane [36];
    float time_plane[36];
    int   max_nh_plane;
    int   nfaces;
    int   nh_face  [18][4];                 // n(hits) per face
    int   np_face  [18][4];                 // n(panels) / face
    int   time_face[18][4];                 // time / face
    //    int   max_np_face;                  // max number of panels/face
    int   max_nh_face;                  // max number of hits/face
    float tmin;                         // min hit time
    float tmax;                         // max hit time
  };
  
public:
					// pointers to the data blocks used
					// 0: TPR, 1: CPR
  TStnTimeClusterBlock*     fTcBlock;
  TStrawHitBlock*           fShBlock;
  TComboHitBlock*           fChBlock;
  TStrTrackBlock*           fTrackBlock;

  TrkPanelMap*              fTpm;
  // TCaloChannelMap*          fCaloChannelMap;
  // TCrvChannelMap*           fCrvChannelMap;
  
  int                       fRunNumber;

  int                       fMaxEvent;   // for X-axis truncation
  int                       fNEvents;
  int                       fNTc;
  int                       fNCh;
  int                       fNSh;
  int                       fNTracks;
 					
  Hist_t*                   fHist;       // histograms to be filled

  std::vector<tc_param_t>   fListOfTcParam;

  float           dt05[36];            // ???
  float           t05 [36];
  float           dtpl[36];
  int             n05 [36];
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TPanelDtAnaModule(const char* name="PanelDtAna", const char* title="Stntuple PanelDtAna");
  ~TPanelDtAnaModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// overloaded methods of TStnModule
//-----------------------------------------------------------------------------
  virtual int     BeginJob()           override;
  virtual int     BeginRun()           override;
  virtual int     Event   (int ientry) override;
  virtual int     EndJob  ()           override;
//-----------------------------------------------------------------------------
// other methods
//-----------------------------------------------------------------------------
  int              BookChannelHistograms(ChannelHist_t* Hist, TrkIndex_t* Index, TFolder* Folder);
  int              BookPanelHistograms  (PanelHist_t*   Hist, TrkIndex_t* Index, TFolder* Folder);
  int              BookSlotHistograms   (SlotHist_t*    Hist, TrkIndex_t* Index, TFolder* Folder);
  int              BookSelHistograms    (SelHist_t*     Hist, TrkIndex_t* Index, TFolder* Folder);
  int              BookTcHistograms     (TcHist_t*      Hist, TrkIndex_t* Index, TFolder* Folder);
  int              BookHistograms       (Hist_t*        Hist, TFolder* Folder);

  int              FillChannelHistograms(ChannelHist_t* Hist, TrkIndex_t* Index, TStrawHit* Sh);
  int              FillPanelHistograms  (PanelHist_t*   Hist, TrkIndex_t* Index, TStrawHit* Sh);
  int              FillSlotHistograms   (SlotHist_t*    Hist, TrkIndex_t* Index, TStrawHit* Sh);
  int              FillSelHistograms    (SelHist_t*     Hist, TrkIndex_t* Index, TStrawHit* Sh);
  
  int              FillTcHistograms     (TcHist_t*      Hist, TrkIndex_t* Index, TStnTimeCluster* Tc,
                                         tc_param_t* Tcp);

  int              FillHistograms       ();

  int              CalculateMissingParameters   ();
  
  void             Debug();

  int              DisplayEvent();
  int              InitEvd();

  //  int              PrintTracks();
  
  int              PrintHistograms   (int ISet);
  int              PrintNoisyChannels(float Percentage = 0.01);

  ClassDefOverride(murat::TPanelDtAnaModule,0)
};
}
#endif
