///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef murat_ana_TVdetAnaModule_hh
#define murat_ana_TVdetAnaModule_hh

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

#include "Stntuple/loop/TStnModule.hh"

#include "Stntuple/obj/TGenpBlock.hh"
#include "Stntuple/obj/TSimpBlock.hh"
#include "Stntuple/obj/TVdetDataBlock.hh"

#include "Stntuple/base/TStnArrayI.hh"

#include "Stntuple/alg/TStnTrackID.hh"
#include "Stntuple/alg/TEmuLogLH.hh"

class TVdetAnaModule: public TStnModule {
public:
//-----------------------------------------------------------------------------
//  histograms
//-----------------------------------------------------------------------------
  struct VdetHitHist_t {
    TH1F*    fIndex;
    TH1F*    fPdgCode;		       //
    TH1F*    fGenCode;		       // generator code
    TH1F*    fEnergy;
    TH1F*    fTime;
  };

  struct EventHist_t {
    TH1F*    fNVdetHits;		// N virtual detector hits in the event
    TH1F*    fNHitsTF;
    TH1F*    fNHitsTB;
    TH1F*    fMomTF;
    TH1F*    fMomTB;
    TH1F*    fMomLoss;
  };

//-----------------------------------------------------------------------------
  enum { kNEventHistSets    = 100 };
  enum { kNVdetHitHistSets  = 100 };

  struct Hist_t {
    EventHist_t*    fEvent   [kNEventHistSets];
    VdetHitHist_t*  fVdetHit [kNVdetHitHistSets];
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used
  TVdetDataBlock*       fVdetDataBlock;
					// histograms filled
  Hist_t                fHist;

  int                   fNVdetHits;
  int                   fNHitsTF;
  int                   fNHitsTB;

  float                 fMomTF;
  float                 fMomTB;

  int                   fPdgCode;
  int                   fGeneratorCode;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TVdetAnaModule(const char* name="VdetAna", const char* title="VdetAna");
  ~TVdetAnaModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Hist_t*            GetHist        () { return &fHist;        }
//-----------------------------------------------------------------------------
// setters
//-----------------------------------------------------------------------------
  void    SetPdgCode      (int Code) { fPdgCode       = Code; }
  void    SetGeneratorCode(int Code) { fGeneratorCode = Code; }
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
  void    BookVdetHitHistograms  (VdetHitHist_t*  Hist, const char* Folder);
  void    BookEventHistograms    (EventHist_t*    Hist, const char* Folder);

  void    FillVdetHitHistograms  (VdetHitHist_t* Hist, TVdetHitData*  Hit);
  void    FillEventHistograms    (EventHist_t*    Hist);

  void    BookHistograms();
  void    FillHistograms();


  void    Debug();
//-----------------------------------------------------------------------------
// test
//-----------------------------------------------------------------------------
  void    Test001();

  ClassDef(TVdetAnaModule,0)
};

#endif
