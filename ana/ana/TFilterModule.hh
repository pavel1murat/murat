///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef murat_ana_TFilterModule_hh
#define murat_ana_TFilterModule_hh

#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TStnTrackSeedBlock.hh"
#include "Stntuple/obj/TStnHelixBlock.hh"
#include "Stntuple/obj/TTrackStrawHitBlock.hh"
#include "Stntuple/obj/TStnClusterBlock.hh"
#include "Stntuple/obj/TCalDataBlock.hh"
#include "Stntuple/obj/TStrawHitBlock.hh"
#include "Stntuple/obj/TStepPointMCBlock.hh"
#include "Stntuple/obj/TGenpBlock.hh"
#include "Stntuple/obj/TSimpBlock.hh"

#include "Stntuple/base/TStnArrayI.hh"

//-----------------------------------------------------------------------------
// among other things, includes definitions of various histogram structures - 
// TrackHist_t etc
//-----------------------------------------------------------------------------
#include "murat/ana/TAnaModule.hh"

namespace murat {
class TFilterModule: public murat::TAnaModule {
public:
//-----------------------------------------------------------------------------
//  histograms
//-----------------------------------------------------------------------------
  enum { kNEventHistSets         = 100 };

  struct Hist_t {
    EventHist_t*           fEvent      [kNEventHistSets  ];
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
  TString               fTrackBlockName;      // 
  TString               fTrackStrawHitBlockName; // 
                                             // pointers to the data blocks used
  TStnTrackBlock*       fTrackBlock;
  TStnTrackSeedBlock*   fTrackSeedBlock;
  TStnHelixBlock*       fHelixBlock;
  TStnTimeClusterBlock* fTimeClusterBlock;
  // TTrackStrawHitBlock* fTrackStrawHitBlock;
  TStnClusterBlock*     fClusterBlock;
  TCalDataBlock*        fCalDataBlock;
  // TStrawHitBlock*      fStrawHitBlock;
  TStepPointMCBlock*    fVDetBlock;
  TGenpBlock*           fGenpBlock;
  TSimpBlock*           fSimpBlock;
                                        // additional parameters (assume nhelices/ntracks < 20)
                                        // histograms filled
  Hist_t               fHist;
                                        // fTrackNumber[i]: track number, 
                                        // corresponding to OBSP particle #i
                                        // or -1
  TStnTrack*           fTrack;

  TrackPar_t           fTrackPar[20];
  SimPar_t             fSimPar;

  int                  fMinNTracks;
  int                  fEventPassedSelections;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TFilterModule(const char* name="murat_Filter", const char* title="Filter");
  ~TFilterModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void            SetTrackBlockName(const char* Name) { fTrackBlockName = Name; }
  void            SetMinNTracks    (int         N   ) { fMinNTracks     = N;    }
//-----------------------------------------------------------------------------
// overloaded methods of TStnModule
//-----------------------------------------------------------------------------
  virtual int     BeginJob();
  virtual int     Event   (int ientry);
  virtual int     EndJob  ();
//-----------------------------------------------------------------------------
// other methods
//-----------------------------------------------------------------------------
  virtual void    BookEventHistograms(HistBase_t* Hist, const char* Folder);
  void    BookHistograms();

  virtual void    FillEventHistograms(HistBase_t* Hist, EventPar_t* Evp);
  void    FillHistograms();

  void    Debug();

  ClassDef(murat::TFilterModule,0)
};
}
#endif
