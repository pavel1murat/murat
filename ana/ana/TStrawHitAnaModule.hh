///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef murat_ana_TStrawHitAnaModule_hh
#define murat_ana_TStrawHitAnaModule_hh

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

#include "Stntuple/loop/TStnModule.hh"

#include "Stntuple/obj/TStrawDataBlock.hh"
#include "Stntuple/obj/TSimpBlock.hh"

#include "Stntuple/base/TStnArrayI.hh"

#include "Stntuple/alg/TStnTrackID.hh"

#include "murat/ana/HistBase_t.h"
#include "murat/ana/AnaDefs.hh"

class TStrawHitAnaModule: public TStnModule {
public:
  struct  StrawHitPar_t {
    float    fSppTime;                  // mother sim particle production time
  };
//-----------------------------------------------------------------------------
//  histograms
//-----------------------------------------------------------------------------
  struct StrawHitHist_t : public HistBase_t {
    TH1F*    fGeneratorCode;		       //
    TH1F*    fPdgCode;		       //
    TH1F*    fMotherPdgCode;		       //
    TH1F*    fMcMomentum;
    TH1F*    fEnergy;
    TH1F*    fTime;
    TH1F*    fSppTime[3];
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
    TH1F* fNStrawHits [4];
    TH1F* fNProtonHits[2];
    TH1F* fNStationsWithHits;
    TH1F* fDeltaSt;
    TH1F* fInstLum;
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
  TStrawDataBlock*      fStrawHitDataBlock;
  TSimpBlock*           fSimpBlock;

  int                   fNStrawHits;
  int                   fNsh200;          // n(straw hits) T>200
  int                   fNsh500;          // n(straw hits) T>500
  int                   fNProtonStrawHits;

					// histograms filled
  Hist_t                fHist;
  StrawHitPar_t         fStrawHitPar[kMaxNStrawHits];

  int                   fNStations;
  int                   fFirstStation;
  int                   fLastStation;
  int                   fNStationsWithHits;
  int                   fNHitsPerStation[50];
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TStrawHitAnaModule(const char* name="StrawHitAna", const char* title="StrawHitAna");
  ~TStrawHitAnaModule();
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
  void    FillStrawHitHistograms(HistBase_t* Hist, TStrawHitData* Hit, StrawHitPar_t* Shp);

  void    BookHistograms();
  void    FillHistograms();


  void    Debug();
  void    PrintStrawHit(TStrawHitData* Hit, StrawHitPar_t* Shp, const char* Option = "");
//-----------------------------------------------------------------------------
// test
//-----------------------------------------------------------------------------
  void    Test001();

  ClassDef(TStrawHitAnaModule,0)
};

#endif
