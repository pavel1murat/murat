///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef __murat_ana_TCrvTimeAnaModule_hh__
#define __murat_ana_TCrvTimeAnaModule_hh__

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

#include "Stntuple/loop/TStnModule.hh"

#include "Stntuple/obj/TComboHitBlock.hh"

#include "Stntuple/obj/TCaloHitBlock.hh"
#include "Stntuple/obj/TCaloRecoDigiBlock.hh"
#include "Stntuple/obj/TStnClusterBlock.hh"

#include "Stntuple/obj/TCrvClusterBlock.hh"
#include "Stntuple/obj/TCrvPulseBlock.hh"

#include "Stntuple/obj/TStnTimeClusterBlock.hh"
#include "Stntuple/obj/TStrTrackBlock.hh"

#include "Stntuple/geom/TCrvChannelMap.hh"
#include "Stntuple/geom/TCaloChannelMap.hh"
#include "Stntuple/geom/TrkPanelMap.hh"

#include "murat/gui/TEvdCrv.hh"

namespace murat {
class TCrvTimeAnaModule: public TStnModule {
  
public:
  
  enum {
    kNStations         = 18,
    kNPlanes           = 36,
    kNPanelsPerStation = 12,
    kNCrystals         = 674*2,
    kNCaloChannels     = kNCrystals*2,   // calorimeter
  };


  enum {
    kNEventHistSets      = 100,
    kNCrvHistSets        = 100,
    kNCrvdHistSets       = 100,
    kNCrvcHistSets       = 100,
    kNCrvpHistSets       = 100,
    kNRocHistSets        =  20,
    kNCalhHistSets       = 100,
    kNCalcHistSets       = 100,
    kNTrkHistSets        = 100,
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
  
  struct CrvIndex_t {
    int sel   {-1};
    int sbid  {-1};
    int sipm  {-1};
    int och   {-1};                           // offline channel - 4*sbid+sipm
    int roc   {-1};                           // 1-18 ??? 
    int feb   {-1};                           // 1-24 in offline domain
    int ch    {-1};                           // channel within the FEB (0-63)
    int sector{-1};
    int module{-1};
    int layer {-1};
  };
  
  struct CaloIndex_t {
    int sel   {-1};
    int crate {-1};
    int disk  {-1};
    int cid   {-1};                   //
  };

  struct crvb_param_t {
    int   intime;
    int   is;                           // sector number (for now, 0..10)
    int   sector_type;
    // float dt10;
    float dt32;
    float dt20;                         // resolution
    float dt31;                         // resolution

    float t;                            // estimated pulse time
    float trk_t;
    float trk_x;
    float trk_z;
    float dt;
    float dx;
    float dz;
    
    TEvdCrvCounter* counter;
  };
  
  struct trk_param_t {
    int                       intime;
    int                       intime_calc;
    int                       intime_crvc;
    float                     dtmin_tc;           // from the closest TC
    float                     x_disk[2];
    float                     y_disk[2];
    TStnCluster*              calc[2];            // closest, for each disk
    float                     dtmin_calc[2];      // from the closest CALC
    float                     dx_calc[2];
    float                     dy_calc[2];
    TStnTimeCluster*          tc;
    float                     dtmin_crvc;         // from the closest CRVC
    TCrvCoincidenceCluster*   crvc;
    float                     dxdz;
    float                     dydz;
    float                     xcrv;
    float                     zcrv;
    float                     dx_crvc;
    float                     dz_crvc;
  };

  struct crystal_t {
    int fCid;
    std::vector<TCaloHit*> fHits;
    
    int       Cid     () { return fCid; }
    int       NHits   () { return (int) fHits.size(); }
    TCaloHit* Hit(int i) { return fHits[i]; }
  };
                                        // a boilerplate for hit bookkeeping

  std::vector<crystal_t>  fCrystals;
  std::vector<crystal_t*> fHitCrystals;
  int                     fMaxHitsPerCrystal;

//-----------------------------------------------------------------------------
// CRV
//-----------------------------------------------------------------------------
  
  // CRV counter has 4 sipms
  struct CrvBar_t {
    std::vector<TCrvRecoPulse*>  list_of_pulses[4];
  };

  // a 
  struct CrvLayer_t {
    std::vector<CrvBar_t> counters;
  };
  
  struct CrvModule_t {
    CrvLayer_t   layers[4];
  };
  
  struct CrvSector_t {
    int n_modules;
    std::vector<CrvModule_t> modules;
  };
  
  struct CrvData_t {
    CrvSector_t sectors[21];

    CrvSector_t*  sector(int i) { return &sectors[i]; }
  } crv_data;
//-----------------------------------------------------------------------------
// histogram structures
//-----------------------------------------------------------------------------
  struct CrvsHist_t {
    TH1F* h_dz;
    TH1F* h_dx;
    TH1F* h_dt;
  };
    
  struct CrvHist_t {
    CrvsHist_t sector[11];
  };
  
  struct TrkHist_t {
    TH1F*         h_nhits;
    TH1F*         h_chi2d;
    TH1F*         h_t0;
    TH1F*         h_dxdz;
    TH1F*         h_dydz;

    TH1F*         h_dt_tc;
    
    TH1F*         h_xcrv;
    TH1F*         h_zcrv;
    TH1F*         h_dt_crvc;
    TH1F*         h_dx_crvc;
    TH1F*         h_dz_crvc;
    TH2F*         h_dx_crvc_vs_dxdy;
    TH2F*         h_dz_crvc_vs_dzdy;
  };

  struct Hist_t {
    TH1F*       h_ntrk;
    TH2F*       h_dt20_vs_sbid;         // CRV
    TH2F*       h_dt31_vs_sbid;         // CRV
    TrkHist_t*  trk[kNTrkHistSets];
    CrvHist_t*  crv[kNCrvHistSets];
  };

//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used
					// 0: TPR, 1: CPR
  TStnTimeClusterBlock*     fTcBlock;
  TCrvClusterBlock*         fCrvcBlock;
  TCrvPulseBlock*           fCrvpBlock;
  TComboHitBlock*           fChBlock;
  TStrTrackBlock*           fTrackBlock;
  TStnClusterBlock*         fCaloClusterBlock;
  TCaloHitBlock*            fCaloHitBlock;
  TCaloRecoDigiBlock*       fCaloRecoDigiBlock;

  TrkPanelMap*              fTpm;
  TCaloChannelMap*          fCaloChannelMap;
  TCrvChannelMap*           fCrvChannelMap;
  
  int                       fRunNumber;

  int                       fMaxEvent;   // for X-axis truncation
  int                       fNEvents;
 					
  Hist_t*                   fHist;       // histograms to be filled

  int                       fNTracks;
  int                       fNTc;
  int                       fNCrvc;
  int                       fNCrvp;
  int                       fNCrvd;

  std::vector<trk_param_t>  fListOfTrkParam;
  std::vector<crvb_param_t> fListOfCrvbParam;

  fit_result_t              fFr[kNCaloChannels];
  fit_result_t*             fFrRef;

  murat::TEvdCrv*           fEvdCrv;

  //  calorimeter_t             fCalo;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TCrvTimeAnaModule(const char* name="CrvTimeAna", const char* title="Stntuple CrvTimeAna");
  ~TCrvTimeAnaModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Hist_t*  GetHist        () { return fHist;        }
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  // void     SetPdgCode    (int Code ) { fPdgCode     = Code ; }
  // void     SetProcessCode(int Code ) { fProcessCode = Code ; }
  //  void     SetUseAllPulses(int Flag) { fUseAllPulses = Flag; }
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
  int              BookCrvsHistograms (CrvsHist_t*   Hist, CrvIndex_t*  Index, TFolder* Folder);
  int              BookTrkHistograms  (TrkHist_t*    Hist, TrkIndex_t*  Index, TFolder* Folder);
  int              BookHistograms     (Hist_t*       Hist, TFolder*     Folder);

  int              CalculateMissingParameters   ();
  int              CalculateMissingTrkParameters();
  
  int              FillCrvsHistograms (CrvsHist_t* Hist, TEvdCrvCounter* Crvb, crvb_param_t* Crvbp);
  int              FillTrkHistograms  (TrkHist_t*  Hist, TStrTrack*      Trk , trk_param_t*  Tp   );
  int              FillHistograms     ();

  int              FitCaloTimeOffsets (float TMin = 1., float TMax = -1.);
  //  int              PrintTimeCorrections();

  void             Debug();
  int              DisplayEvent();
 
  int              InitEvd();

  int              PrintTracks();
  int              PrintCrvBars();

  ClassDefOverride(murat::TCrvTimeAnaModule,0)
};
}
#endif
