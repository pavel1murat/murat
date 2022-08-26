///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef murat_ana_TVstAnaModule_hh
#define murat_ana_TVstAnaModule_hh

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

#include "Stntuple/loop/TStnModule.hh"

#include "Stntuple/obj/TStrawHitBlock.hh"
#include "Stntuple/obj/TSimpBlock.hh"

#include "Stntuple/base/TStnArrayI.hh"

#include "Stntuple/alg/TStnTrackID.hh"

#include "murat/ana/HistBase_t.h"
#include "murat/ana/AnaDefs.hh"

class TVstAnaModule: public TStnModule {
public:
  struct  StrawHitPar_t {
    float    fSppTime;                  // mother sim particle production time
  };
//-----------------------------------------------------------------------------
//  histograms
//-----------------------------------------------------------------------------
  struct StrawHitHist_t : public HistBase_t {
    TH1F*    fEDep;
    TH1F*    fTime[2];
    TH1F*    fTOT[2];
    TH1F*    fDt;
    TH1F*    fStation;
    TH1F*    fFace;
    TH1F*    fPanel;
    TH1F*    fLayer;
    TH1F*    fStraw;
    TH1F*    fPreamp;
  };

  struct EventHist_t : public HistBase_t {
    TH1F* fRunNumber;
    TH1F* fEventNumber;
    TH1F* fInstLum;
    TH1F* fNStrawHits[2];

    TH1F* fNHitsPerPanel[6];
    TH1F* fNHitPanels;
    TH1F* fNHitsPerFace[2];
    TH1F* fNHitFaces;

    TH1F* fNGoodStrawHits;
    TH1F* fNGoodHitsPerPanel[6];
    TH1F* fNGoodHitPanels;
    TH1F* fNGoodHitsPerFace[2];
    TH1F* fNGoodHitFaces;
  };

  struct HitCounters_t {
    int fNHits;				// total
    int fNHitFaces;
    int fNHitPanels;
    int fNHitsPerPanel[6];
    int fNHitsPerFace [2];

    void Clear() {
      fNHits      = 0;
      fNHitFaces  = 0;
      fNHitPanels = 0;
      for (int i=0; i<6; i++) { fNHitsPerPanel[i] = 0; }
      for (int i=0; i<2; i++) { fNHitsPerFace [i] = 0; }
    }
  };
//-----------------------------------------------------------------------------
  enum { kNEventHistSets     =   10 };
  enum { kNStrawHitHistSets  = 1000 };

  enum { kMaxNStrawHits      = 20000 };

  struct Hist_t {
    EventHist_t*     fEvent   [kNEventHistSets];
    StrawHitHist_t*  fStrawHit[kNStrawHitHistSets];
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used
  TStrawHitBlock*       fStrawHitBlock;
  int                   fNStrawHits;
					// histograms filled
  Hist_t                fHist;
  StrawHitPar_t         fStrawHitPar[kMaxNStrawHits];

  HitCounters_t         fAllHits;
  HitCounters_t         fGoodHits;
//-----------------------------------------------------------------------------
// coming calibration constants
//-----------------------------------------------------------------------------
					// face for a given panel
  const int             kFace[6] = { 0, 1, 0, 1, 0, 1};

  float                 fMinEDep[6];
  float                 fMaxEDep[6];
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TVstAnaModule(const char* name="VstAna", const char* title="VstAna");
  ~TVstAnaModule();
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
  void    BookEventHistograms   (HistBase_t* Hist, const char* Folder);
  void    BookStrawHitHistograms(HistBase_t* Hist, const char* Folder);

  void    FillEventHistograms   (HistBase_t* Hist);
  void    FillStrawHitHistograms(HistBase_t* Hist, TStrawHit* Hit, StrawHitPar_t* Shp);

  void    BookHistograms();
  void    FillHistograms();

  void    ControlBar();                 // **MENU**


  void    Debug();
  void    PrintStrawHitBlock(); 
  void    PrintStrawHit(TStrawHit* Hit, const char* Option = "");
//-----------------------------------------------------------------------------
// test
//-----------------------------------------------------------------------------
//  void    Test001();

  ClassDef(TVstAnaModule,0)
};

#endif
