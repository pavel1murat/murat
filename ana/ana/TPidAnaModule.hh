///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef murat_ana_TPidAnaModule_hh
#define murat_ana_TPidAnaModule_hh

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

#include "Stntuple/loop/TStnModule.hh"

#include "Stntuple/obj/TGenpBlock.hh"
#include "Stntuple/obj/TSimpBlock.hh"
#include "Stntuple/obj/TStnPidBlock.hh"
#include "Stntuple/obj/TStnTrackBlock.hh"

#include "Stntuple/base/TStnArrayI.hh"

#include "Stntuple/alg/TStnTrackID.hh"
#include "Stntuple/alg/TEmuLogLH.hh"

class TPidAnaModule: public TStnModule {
public:

#include "murat/ana/TrackPar_t.hh"
//-----------------------------------------------------------------------------
//  histograms
//-----------------------------------------------------------------------------
  struct PidHist_t {
    TH1F*    fTrkNumber;

    TH1F*    fNMatched;		       //
    TH1F*    fNMatchedAll;		       // generator code

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

  struct EventHist_t {
    TH1F*    fNTracks;
  };

//-----------------------------------------------------------------------------
  enum { kNEventHistSets    = 100 };
  enum { kNPidHistSets      = 100 };

  struct Hist_t {
    EventHist_t*    fEvent   [kNEventHistSets];
    PidHist_t*      fPid     [kNPidHistSets  ];
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used 0:DEM, 1:DMM

  TStnPidBlock*         fPidDataBlock[2];
  TStnTrackBlock*       fTrackBlock  [2];

					// additional track parameters (assume ntracks < 20)
  int                   fNID;
  TStnTrackID*          fTrackID [20];
  TrackPar_t            fTrackPar[20];
					// histograms filled
  Hist_t                fHist;
  int                   fNTracks[2];
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TPidAnaModule(const char* name="PidAna", const char* title="PidAna");
  ~TPidAnaModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Hist_t*            GetHist        () { return &fHist;        }
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
  void    BookEventHistograms    (EventHist_t*    Hist, const char* Folder);

  void    FillPidHistograms      (PidHist_t* Hist, TStnPid*  Hit);
  void    FillEventHistograms    (EventHist_t*    Hist);

  void    BookHistograms();
  void    FillHistograms();


  void    Debug();
//-----------------------------------------------------------------------------
// test
//-----------------------------------------------------------------------------
  void    Test001();

  ClassDef(TPidAnaModule,0)
};

#endif
