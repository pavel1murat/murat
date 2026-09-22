///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef __murat_ana_TCaloTimeAnaModule_hh__
#define __murat_ana_TCaloTimeAnaModule_hh__

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

namespace murat {
class TCaloTimeAnaModule: public TStnModule {
  
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
    kNCalhHistSets       = 100,
    kNCalcHistSets       = 100,
  };
//-----------------------------------------------------------------------------
// indices
//-----------------------------------------------------------------------------
  struct CaloIndex_t {
    int sel   {-1};
    int crate {-1};
    int disk  {-1};
    int cid   {-1};                   //
  };
  
  struct calc_param_t {
    float                     dtmin_tc;           // from the closest TC
    TStnTimeCluster*          tc;
    float                     dtmin_crvc;     // from the closest CRVC
    float                     dtmin_trk;                // from the closest track
  };
                                                  // for now, a placeholder
  struct crystal_t {
    int fCid;
    std::vector<TCaloHit*>      fHits;
    std::vector<TCaloRecoDigi*> fCrd[2];
    
    int       Cid     () { return fCid; }
    int       NHits   () { return (int) fHits.size(); }
    TCaloHit* Hit(int i) { return fHits[i]; }
  };
                                        // a boilerplate for hit bookkeeping

  std::vector<crystal_t>  fCrystals;
  std::vector<crystal_t*> fHitCrystals;
  int                     fMaxHitsPerCrystal;
//-----------------------------------------------------------------------------
// histogram structures
//-----------------------------------------------------------------------------
  struct CalhHist_t {
    TH1F* h_dt;
  };

  struct CalcHist_t {
    TH1F*         h_edep;
    TH1F*         h_dt_tc;
    TH1F*         h_dt_crvc;
    TH2F*         h_dt_crvc_vs_dt_tc;
  };

  struct DiskHist_t {
    TH1F*         h_ch;
    TH2F*         h_board_vs_dt;
  };

  struct Hist_t {
    DiskHist_t* disk[2];                // 2 disks
    TH2F*       h_dt_vs_sipmid;
    TH2F*       h_dt10_vs_cid;
    TH1F*       h_dtpp[2];
    TH2F*       h_nsipms_vs_cid;
    TH1F*       h_sipmid;               // occupancy offline channel
    TH2F*       h_n2_vs_n1;
    TH1F*       h_ntrk[2];
    CalhHist_t* calh   [kNCalhHistSets];
    CalcHist_t* calc   [kNCalcHistSets];
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

  //  TrkPanelMap*              fTpm;
  TCaloChannelMap*          fCaloChannelMap;
  // TCrvChannelMap*           fCrvChannelMap;
  
  int                       fRunNumber;

  int                       fMaxEvent;   // for X-axis truncation
  int                       fNEvents;
 					
  Hist_t*                   fHist;       // histograms to be filled

  int                       fNCalh10[2]; // N(hits above 10 MeV)
  int                       fNCalh;
  int                       fNCalord;
  int                       fNCaloClusters;

  int                       fNCcDisk[2];

  std::vector<calc_param_t> fListOfCalcParam;

  fit_result_t              fFr[kNCaloChannels];
  fit_result_t*             fFrRef;

  //  calorimeter_t             fCalo;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TCaloTimeAnaModule(const char* name="CaloTimeAna", const char* title="murat::CaloTimeAna");
  ~TCaloTimeAnaModule();
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
  int              BookCalcHistograms (CalcHist_t*   Hist, CaloIndex_t* Index, TFolder* Folder);
  int              BookCalhHistograms (CalhHist_t*   Hist, CaloIndex_t* Index, TFolder* Folder);
  int              BookDiskHistograms (DiskHist_t*   Hist, CaloIndex_t* Index, TFolder* Folder);
  int              BookHistograms     (Hist_t*       Hist, TFolder*     Folder);

  int              CalculateMissingParameters   ();
  
  int              FillCalcHistograms (CalcHist_t* Hist, TStnCluster*   Calc, calc_param_t* Cp);
  int              FillDiskHistograms (DiskHist_t* Hist, TCaloRecoDigi* Calrd);
  int              FillHistograms     ();

  int              FitCaloTimeOffsets (float TMin = 1., float TMax = -1.);
  //  int              PrintTimeCorrections();

  void             Debug();
  
  int              PrintTracks();

  ClassDefOverride(murat::TCaloTimeAnaModule,0)
};
}
#endif
