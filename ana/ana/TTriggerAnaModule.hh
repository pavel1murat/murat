///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef murat_ana_TTriggerAnaModule_hh
#define murat_ana_TTriggerAnaModule_hh

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

#include "Stntuple/obj/TGenpBlock.hh"
#include "Stntuple/obj/TSimpBlock.hh"
#include "Stntuple/obj/TStnTimeClusterBlock.hh"
#include "Stntuple/obj/TStnHelixBlock.hh"
#include "Stntuple/obj/TStnTrackSeedBlock.hh"
#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TStnClusterBlock.hh"
#include "Stntuple/obj/TStnTriggerBlock.hh"

#include "Stntuple/base/TStnArrayI.hh"

#include "murat/ana/HistBase_t.h"
#include "murat/ana/TrackPar_t.hh"
#include "murat/ana/HelixPar_t.hh"

#include "murat/ana/AnaDefs.hh"

#include "murat/ana/TAnaModule.hh"

namespace murat {

class TTriggerAnaModule: public TAnaModule {
public:
//-----------------------------------------------------------------------------
//  histograms
//-----------------------------------------------------------------------------
  struct TimeClusterHist_t : public HistBase_t {
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

  // struct TrackHist_t : public HistBase_t {
  //   TH1F* fP;
  //   TH1F* fNActive;
  //   TH1F* fChi2Dof;
  //   TH1F* fT0;
  //   TH1F* fD0;
  //   TH1F* fZ0;
  //   TH1F* fTanDip;
  //   TH1F* fAlgMask;
  // };

  struct TriggerHist_t : public HistBase_t {
    TH1F* fBits;
  };

  // struct EventHist_t : public HistBase_t {
  //   TH1F* fRunNumber;
  //   TH1F* fEventNumber;
  //   TH1F* fNTimeClusters;
  //   TH1F* fNHelices   ;
  //   TH1F* fNTrackSeeds[2];
  //   TH1F* fNGoodSeeds ;
  //   TH1F* fNTracks    ;
  //   TH1F* fPassed     ;
  //   TH1F* fMcMom      ;
  //   TH1F* fMcCosTh    ;
  // };

//-----------------------------------------------------------------------------
// histograms to be filled
//-----------------------------------------------------------------------------
  enum { kNTimeClusterHistSets = 100 };
  enum { kNHelixHistSets       = 100 };
  enum { kNTrackSeedHistSets   = 100 };
  enum { kNTrackHistSets       = 200 };
  enum { kNTriggerHistSets     = 100 };
  enum { kNEventHistSets       = 100 };

  struct Hist_t {
    TimeClusterHist_t*  fTimeCluster[kNTimeClusterHistSets];
    HelixHist_t*        fHelix      [kNHelixHistSets      ];
    TrackSeedHist_t*    fTrackSeed  [kNTrackSeedHistSets  ];
    TrackHist_t*        fTrack      [kNTrackHistSets      ];
    TriggerHist_t*      fTrigger    [kNTriggerHistSets    ];
    EventHist_t*        fEvent      [kNEventHistSets      ];
  } fHist;
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used
  TGenpBlock*              fGenpBlock;
  TSimpBlock*              fSimpBlock;
  TStnTimeClusterBlock*    fTimeClusterBlock;
  TStnHelixBlock*          fHelixBlock;
  TStnTrackSeedBlock*      fTrackSeedBlock;
  TStnTrackBlock*          fTrackBlock;
  TStnClusterBlock*        fClusterBlock;
  TStnTriggerBlock*        fTriggerBlock;

  SimPar_t                 fSimPar;     // additional parameters of the simulated MC particle

  int                      fPassed;
  int                      fNGenp;
  int                      fNTracks;
  int                      fNTimeClusters;
  int                      fNHelices;
  int                      fNTrackSeeds[10];
  int                      fNGoodSeeds;
  int                      fNGoodTracks;
  TrackPar_t               fTrackPar   [10];

  double                   fMinTrigMom;
  TGenParticle*            fParticle;
  double                   fMcMom;
  double                   fMcCosTh;
  double                   fWeight;

  TString                  fTrackBlockName;
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
  void               SetTrackBlockName(const char* Name) { fTrackBlockName = Name ; }
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
  void    BookTimeClusterHistograms(HistBase_t* Hist, const char* Folder);
  void    BookHelixHistograms      (HistBase_t* Hist, const char* Folder);
  void    BookTrackSeedHistograms  (HistBase_t* Hist, const char* Folder);
  //  void    BookTrackHistograms      (HistBase_t* Hist, const char* Folder);
  void    BookTriggerHistograms    (HistBase_t* Hist, const char* Folder);
  //  void    BookEventHistograms      (HistBase_t* Hist, const char* Folder);
  void    BookHistograms();

  void    FillTimeClusterHistograms(HistBase_t* Hist, TStnTimeCluster* TPeak);
  void    FillHelixHistograms      (HistBase_t* Hist, TStnHelix*       Helix, HelixPar_t* Help, double Weight = 1);
  void    FillTrackSeedHistograms  (HistBase_t* Hist, TStnTrackSeed*   Seed , double Weight = 1);
  //  void    FillTrackHistograms      (HistBase_t* Hist, TStnTrack*       Trk);
  void    FillTriggerHistograms    (HistBase_t* Hist);
  //  void    FillEventHistograms      (HistBase_t* Hist, double Weight = 1);

  void    FillHistograms();

  //  void    InitTrackPar();


  void    Debug();
//-----------------------------------------------------------------------------
// test
//-----------------------------------------------------------------------------
  void    Test001();

  ClassDef(TTriggerAnaModule,0)
};
}

#endif
