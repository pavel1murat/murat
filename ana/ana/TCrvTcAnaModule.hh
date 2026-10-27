///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef __murat_ana_TCrvTcAnaModule_hh__
#define __murat_ana_TCrvTcAnaModule_hh__

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TCanvas.h"

#include "Stntuple/loop/TStnModule.hh"

#include "Stntuple/obj/TComboHitBlock.hh"

#include "Stntuple/obj/TCaloHitBlock.hh"
#include "Stntuple/obj/TCaloRecoDigiBlock.hh"
#include "Stntuple/obj/TStnClusterBlock.hh"

#include "Stntuple/obj/TCrvCoincidenceCluster.hh"
#include "Stntuple/obj/TCrvClusterBlock.hh"
#include "Stntuple/obj/TCrvPulseBlock.hh"
#include "Stntuple/obj/TCrvDigiBlock.hh"

#include "Stntuple/obj/TStnTimeClusterBlock.hh"
#include "Stntuple/obj/TStrTrackBlock.hh"

#include "Stntuple/geom/TCrvChannelMap.hh"
#include "Stntuple/geom/TCaloChannelMap.hh"
#include "Stntuple/geom/TrkPanelMap.hh"

#include "murat/gui/TEvdCrv.hh"

// subdetector-level alignment table
#include "murat/alg/DetPos_t.hh"

namespace murat {
class TCrvTcAnaModule: public TStnModule {
  
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
    kNCrvsHistSets       = 100,         // by CRV sector
    kNFebHistSets        = 100,
    kNRocHistSets        = 100,
  };
//-----------------------------------------------------------------------------
// indices
//-----------------------------------------------------------------------------
  struct EvtIndex_t {
    int sel;
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
//-----------------------------------------------------------------------------
// histogram structures
//-----------------------------------------------------------------------------
  struct ChannelHist_t {
    TH1F* h_ch;
    TH1F* h_dt;
  };

  struct CrvcHist_t {
    TH1F* h_crvc_tc_dt;
  };

                                        // by sector
  struct CrvsHist_t {
    TH1F* h_crvc_tc_dt;
    TH1F* h_crvp_tc_dt;
  };

  struct CrvpHist_t {
    TH1F* h_ph;
    TH1F* h_npes;
    TH1F* h_time;
    TH1F* h_feb;
    TH1F* h_ch;
    TH1F* h_crvp_tc_dt;
    TH2F* h_crvp_tc_dt_vs_feb;
  };

  struct FebHist_t {
    ChannelHist_t ch[64];
    CrvpHist_t*   crvp;               // for the whole FEB
    TH1F*         h_sbid;
    TH2F*         h_ch_vs_crvp_tc_dt;
  };
  
  struct RocHist_t {
    CrvpHist_t*   crvp;               // for the whole ROC
    FebHist_t*    feb [30];
    TH1F*         h_sbid;
    TH2F*         h_feb_vs_crvp_tc_dt;
  };

  struct EventHist_t {
    TH1F* h_evn;
  };
  
  struct Hist_t {
    EventHist_t evt[kNEventHistSets];
    CrvpHist_t* crvp[10];               // for all
    CrvcHist_t* crvc[10];               // for all
    CrvsHist_t* crvs[kNCrvsHistSets];   // by sector
    RocHist_t*  roc [10];               // by ROC
    TH1F*       h_chid;                  // occupancy offline channel
    TH1F*       h_sbid;
    TH2F*       h_feb_vs_chid;
    TH2F*       h_crvp_tc_dt_vs_chid;
    TH2F*       h_feb_vs_sbid[2];       // one per ROC
  };
//-----------------------------------------------------------------------------
//  additional parameters
//-----------------------------------------------------------------------------
  struct crvc_param_t {
    int    sector;
    float  dt_min;                          // T(CRVP)-T(TC)
  };

  struct crvp_param_t {
    int    sector;
    float  dt_min;                          // T(CRVP)-T(TC)
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
  TCrvDigiBlock*            fCrvdBlock;

  TrkPanelMap*              fTpm;
  TCaloChannelMap*          fCaloChannelMap;
  TCrvChannelMap*           fCrvChannelMap;
  
  int                       fRunNumber;

  int                       fMaxEvent;   // for X-axis truncation
  int                       fNEvents;
 					
  Hist_t*                   fHist;       // histograms to be filled

  int                       fNTc;
  int                       fNTc20;
  int                       fNCrvc;
  int                       fNCrvp;
  int                       fNCrvd;

  fit_result_t              fFr[kNCaloChannels];
  fit_result_t*             fFrRef;

  murat::TEvdCrv*           fEvdCrv;

  DetPos_t                  fDetPos;

  TCanvas*                  fCanvasWf;

  std::vector<crvp_param_t> fListOfCrvpParam;
  std::vector<crvc_param_t> fListOfCrvcParam;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TCrvTcAnaModule(const char* name="CrvTcAna", const char* title="Stntuple CrvTcAna");
  ~TCrvTcAnaModule();
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
  int              BookCrvcHistograms (CrvcHist_t*   Hist, CrvIndex_t*  Index, TFolder* Folder);
  int              BookCrvpHistograms (CrvpHist_t*   Hist, CrvIndex_t*  Index, TFolder* Folder);
  int              BookCrvsHistograms (CrvsHist_t*   Hist, CrvIndex_t*  Index, TFolder* Folder);
  int              BookEventHistograms(EventHist_t*  Hist, EvtIndex_t*  Index, TFolder* Folder);
  int              BookFebHistograms  (FebHist_t*    Hist, CrvIndex_t*  Index, TFolder* Folder);
  int              BookRocHistograms  (RocHist_t*    Hist, CrvIndex_t*  Index, TFolder* Folder);
  int              BookHistograms     (Hist_t*       Hist, TFolder*     Folder);

  int              FillCrvcHistograms (CrvcHist_t*  Hist, TCrvCoincidenceCluster* Crvp, crvc_param_t* Crvpp);
  int              FillCrvpHistograms (CrvpHist_t*  Hist, TCrvRecoPulse*          Crvp, crvp_param_t* Crvpp);
  int              FillCrvsHistograms ();
  int              FillEventHistograms(EventHist_t* Hist);
  int              FillFebHistograms  (FebHist_t*   Hist, TCrvRecoPulse* Crvp, crvp_param_t* Crvpp);
  int              FillRocHistograms  (RocHist_t*   Hist, TCrvRecoPulse* Crvp, crvp_param_t* Crvpp);
  int              FillHistograms     ();

  void             Debug();

  int              DisplayEvent();
  int              InitEvd();

  int              PrintTracks();
  int              PrintCrvBars();
  
  // by defaut open new canvas
  // if there are several waveforms, plot all of them
  void             PlotWaveform(int ChannelID, bool NewCanvas = true);

  ClassDefOverride(murat::TCrvTcAnaModule,0)
};
}
#endif
