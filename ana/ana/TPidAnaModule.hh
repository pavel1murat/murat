///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef murat_ana_TPidAnaModule_hh
#define murat_ana_TPidAnaModule_hh

#include "Stntuple/obj/TGenpBlock.hh"
#include "Stntuple/obj/TSimpBlock.hh"
#include "Stntuple/obj/TStnPidBlock.hh"
#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TStnClusterBlock.hh"

#include "Stntuple/base/TStnArrayI.hh"

#include "Stntuple/alg/TStnTrackID.hh"
#include "Stntuple/alg/TEmuLogLH.hh"

#include "murat/ana/TAnaModule.hh"

namespace murat {
class TPidAnaModule: public TAnaModule {
public:

#include "murat/ana/SimPar_t.hh"
#include "murat/ana/TrackPar_t.hh"
#include "murat/ana/TrackHist_t.hh"

#include "murat/ana/TrackSeedHist_t.hh"
#include "murat/ana/EventHist_t.hh"
//-----------------------------------------------------------------------------
//  histograms
//-----------------------------------------------------------------------------
  struct PidHist_t {
    TH1F*    fTrkNumber;

    TH1F*    fNMatched;		       //
    TH1F*    fNMatchedAll;	       // generator code

    TH1F*    fNUsedOsH;
    TH1F*    fNUsedSsH;
    TH1F*    fNUsedOsD;

    TH1F*    fLogDedxProbEle;
    TH1F*    fLogDedxProbMuo;
    TH1F*    fLhrDedx;

    TH1F*    fDrdsVadim;
    TH1F*    fXdrdsVadim;

    TH1F*    fSumAvik;
    TH1F*    fSq2Avik;

    TH1F*    fDrdsOs;
    TH1F*    fXdrdsOs;

    TH1F*    fSumAvikOs;

    TH1F*    fDrdsSs;
    TH1F*    fXdrdsSs;
  };
//-----------------------------------------------------------------------------
// histograms
//-----------------------------------------------------------------------------
  enum { kNEventHistSets    = 100 };
  enum { kNPidHistSets      = 100 };

  struct Hist_t {
    EventHist_t*    fEvent   [kNEventHistSets];
    PidHist_t*      fPid     [kNPidHistSets  ];
  } fHist;
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used 0:DEM, 1:DMM
  TStnPidBlock*     fPidDataBlock  [2];
  TStnTrackBlock*   fTrackBlock    [2];
  TStnClusterBlock* fClusterBlock;

  TString           fTrackBlockName[2];

  SimPar_t          fSimPar;            // additional parameters of the simulated MC particle
					
  TrackPar_t        fTrackPar[2][20];   // additional track parameters (assume ntracks < 20)
  int               fNTracks [2];
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TPidAnaModule(const char* name="PidAna", const char* title="PidAna");
  ~TPidAnaModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Hist_t*            GetHist          () { return &fHist;        }
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  void    SetTrackBlockName(int I, const char* Name) { fTrackBlockName[I] = Name;}
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
  void    BookPidHistograms      (PidHist_t*  Hist, const char* Folder);
  //  void    BookEventHistograms    (EventHist_t*    Hist, const char* Folder);

  void    FillPidHistograms      (PidHist_t* Hist, TStnPid*  Hit);
  //  void    FillEventHistograms    (EventHist_t*    Hist);

  void    BookHistograms();
  void    FillHistograms();


  void    Debug();
//-----------------------------------------------------------------------------
// test
//-----------------------------------------------------------------------------
  void    Test001();

  ClassDef(murat::TPidAnaModule,0)
};
}
#endif
