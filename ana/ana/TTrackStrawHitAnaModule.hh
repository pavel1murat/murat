///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef murat_ana_TTrackStrawHitAnaModule_hh
#define murat_ana_TTrackStrawHitAnaModule_hh

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

#include "Stntuple/loop/TStnModule.hh"

#include "Stntuple/obj/TTrackStrawHitBlock.hh"
#include "Stntuple/obj/TSimpBlock.hh"

#include "Stntuple/base/TStnArrayI.hh"

#include "Stntuple/alg/TStnTrackID.hh"

#include "murat/ana/HistBase_t.h"
#include "murat/ana/TrackParBase_t.hh"

#include "murat/ana/AnaDefs.hh"

class TTrackStrawHitAnaModule: public TStnModule {
public:
//-----------------------------------------------------------------------------
//  histograms
//-----------------------------------------------------------------------------
  struct TrackStrawHitHist_t : public HistBase_t {
    TH1F*    fMcDoca;		       //
    TH1F*    fEnergy;		       //
    TH1F*    fDriftRadius;		       //
    TH1F*    fMcMomentum;
    TH1F*    fTime;
    TH2F*    fDriftRadiusVsMcDoca;
    TH1F*    fDr; // drift-mcdoca
  };

  struct EventHist_t : public HistBase_t {
    TH1F* fRunNumber;
    TH1F* fEventNumber;
    TH1F* fNTracks;
  };

  struct TrackHist_t : public HistBase_t {
    TH1F* fNHits;
    TH1F* fNActive;
    TH1F* fN500;
  };

  struct TrackPar_t : public TrackParBase_t {
    int fNHits;
    int fNActive;
    int fN500;
  };

//-----------------------------------------------------------------------------
  enum { kNEventHistSets       = 10 };
  enum { kNTrackHistSets       = 10 };
  enum { kNTrkStrawHitHistSets = 10 };

  struct Hist_t {
    EventHist_t*          fEvent        [kNEventHistSets];
    TrackHist_t*          fTrack        [kNTrackHistSets];
    TrackStrawHitHist_t*  fTrackStrawHit[kNTrkStrawHitHistSets];
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used

  TString              fHitBlockName;
  TTrackStrawHitBlock* fTrackStrawHitDataBlock;

  int                   fNTracks;

  TrackPar_t            fTp[100];
					// histograms filled
  Hist_t                fHist;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TTrackStrawHitAnaModule(const char* name="TrackStrawHitAna", const char* title="TrackStrawHitAna");
  ~TTrackStrawHitAnaModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Hist_t*            GetHist        () { return &fHist;        }
//-----------------------------------------------------------------------------
// setters
//-----------------------------------------------------------------------------
  void    SetHitBlockName(const char* Name) { fHitBlockName = Name; }
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
  void    BookEventHistograms        (HistBase_t* Hist, const char* Folder);
  void    BookTrackStrawHitHistograms(HistBase_t* Hist, const char* Folder);
  void    BookTrackHistograms        (HistBase_t* Hist, const char* Folder);

  void    FillEventHistograms        (HistBase_t* Hist);
  void    FillTrackHistograms        (HistBase_t* Hist, TrackParBase_t*     Trk);
  void    FillTrackStrawHitHistograms(HistBase_t* Hist, TTrackStrawHitData* Hit);

  void    BookHistograms();
  void    FillHistograms();


  void    Debug();
//-----------------------------------------------------------------------------
// test
//-----------------------------------------------------------------------------
  void    Test001();

  ClassDef(TTrackStrawHitAnaModule,0)
};

#endif
