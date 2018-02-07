///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef murat_ana_TTriggerAnaModule_hh
#define murat_ana_TTriggerAnaModule_hh

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

#include "Stntuple/loop/TStnModule.hh"

#include "Stntuple/obj/TStnTimePeakBlock.hh"
#include "Stntuple/obj/TStnHelixBlock.hh"
#include "Stntuple/obj/TStnTrackSeedBlock.hh"
#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TStnClusterBlock.hh"

#include "Stntuple/base/TStnArrayI.hh"

#include "murat/ana/HistBase_t.h"
#include "murat/ana/TrackPar_t.hh"

#include "murat/ana/AnaDefs.hh"

class TTriggerAnaModule: public TStnModule {
public:
//-----------------------------------------------------------------------------
//  histograms
//-----------------------------------------------------------------------------
  struct TimePeakHist_t : public HistBase_t {
    TH1F* fTime;
    TH1F* fNHits;
    TH1F* fEnergy;
    TH1F* fR;
  };

  struct HelixHist_t : public HistBase_t {
    TH1F* fNHits;
  };

  struct TrackSeedHist_t : public HistBase_t {
    TH1F* fP;
    TH1F* fNHits;
    TH1F* fChi2Dof;
    TH1F* fD0;
  };

  struct TrackHist_t : public HistBase_t {
    TH1F* fP;
    TH1F* fNActive;
    TH1F* fChi2Dof;
    TH1F* fT0;
    TH1F* fD0;
    TH1F* fZ0;
    TH1F* fTanDip;
    TH1F* fAlgMask;
  };

  struct EventHist_t : public HistBase_t {
    TH1F* fRunNumber;
    TH1F* fEventNumber;
    TH1F* fNTimeClusters;
    TH1F* fNHelices   ;
    TH1F* fNTrackSeeds[2];
    TH1F* fNGoodSeeds ;
    TH1F* fNTracks    ;
    TH1F* fPassed     ;
  };

//-----------------------------------------------------------------------------
  enum { kNTimePeakHistSets  = 100 };
  enum { kNHelixHistSets     = 100 };
  enum { kNTrackSeedHistSets = 100 };
  enum { kNTrackHistSets     = 200 };
  enum { kNEventHistSets     =  10 };

  struct Hist_t {
    TimePeakHist_t*  fTimePeak [kNTimePeakHistSets ];
    HelixHist_t*     fHelix    [kNHelixHistSets    ];
    TrackSeedHist_t* fTrackSeed[kNTrackSeedHistSets];
    TrackHist_t*     fTrack    [kNTrackHistSets    ];
    EventHist_t*     fEvent    [kNEventHistSets    ];
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used

  TStnTimePeakBlock*       fTimePeakBlock;
  TStnHelixBlock*          fHelixBlock;
  TStnTrackSeedBlock*      fTrackSeedBlock;
  TStnTrackBlock*          fTrackBlock;
  TStnClusterBlock*        fClusterBlock;

					// histograms filled
  Hist_t                   fHist;

  int                      fPassed;
  int                      fNTracks;
  int                      fNTimeClusters;
  int                      fNHelices;
  int                      fNTrackSeeds[10];
  int                      fNGoodSeeds;
  int                      fNGoodTracks;
  TrackPar_t               fTrackPar   [10];
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TTriggerAnaModule(const char* name="TriggerAna", const char* title="TriggerAna");
  ~TTriggerAnaModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Hist_t*            GetHist        () { return &fHist;        }
//-----------------------------------------------------------------------------
// setters
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// overloaded methods of TStnModule
//-----------------------------------------------------------------------------
  int     BeginJob();
  int     BeginRun();
  int     Event   (int ientry);
  int     EndJob  ();
//-----------------------------------------------------------------------------
// other methods
//-----------------------------------------------------------------------------
  void    BookTimePeakHistograms (HistBase_t* Hist, const char* Folder);
  void    BookHelixHistograms    (HistBase_t* Hist, const char* Folder);
  void    BookTrackSeedHistograms(HistBase_t* Hist, const char* Folder);
  void    BookTrackHistograms    (HistBase_t* Hist, const char* Folder);
  void    BookEventHistograms    (HistBase_t* Hist, const char* Folder);
  void    BookHistograms();

  void    FillTimePeakHistograms (HistBase_t* Hist, TStnTimePeak*  TPeak);
  void    FillHelixHistograms    (HistBase_t* Hist, TStnHelix*     Helix);
  void    FillTrackSeedHistograms(HistBase_t* Hist, TStnTrackSeed* Seed );
  void    FillTrackHistograms    (HistBase_t* Hist, TStnTrack*     Trk  );
  void    FillEventHistograms    (HistBase_t* Hist);

  void    FillHistograms();

  void    InitTrackPar();


  void    Debug();
//-----------------------------------------------------------------------------
// test
//-----------------------------------------------------------------------------
  void    Test001();

  ClassDef(TTriggerAnaModule,0)
};

#endif
