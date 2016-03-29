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
    TH1F*    fEleTrkNumber;
    TH1F*    fMuoTrkNumber;

    TH1F*    fNMatched;		       //
    TH1F*    fNMatchedAll;		       // generator code

    TH1F*    fNUsedOsEleH;
    TH1F*    fNUsedOsMuoH;
    TH1F*    fNUsedSsEleH;
    TH1F*    fNUsedSsMuoH;

    TH1F*    fNUsedOsEleD;
    TH1F*    fNUsedOsMuoD;

    TH1F*    fLogDedxProbEle;
    TH1F*    fLogDedxProbMuo;
    TH1F*    fLhrDedx;

    TH1F*    fDrdsVadimEle;
    TH1F*    fXdrdsVadimEle;
    TH1F*    fDrdsVadimMuo;
    TH1F*    fXdrdsVadimMuo;

    TH1F*    fSumAvikEle;
    TH1F*    fSumAvikMuo;
    TH1F*    fLhrAvikSum;

    TH1F*    fSq2AvikEle;
    TH1F*    fSq2AvikMuo;

    TH1F*    fDrdsOsEle;
    TH1F*    fXdrdsOsEle;
    TH1F*    fDrdsOsMuo;
    TH1F*    fXdrdsOsMuo;

    TH1F*    fSumAvikOsEle;
    TH1F*    fSumAvikOsMuo;

    TH1F*    fDrdsSsEle;
    TH1F*    fXdrdsSsEle;
    TH1F*    fDrdsSsMuo;
    TH1F*    fXdrdsSsMuo;
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
					// pointers to the data blocks used
  TStnPidBlock*         fPidDataBlock;
  TStnTrackBlock*       fTrackBlock;
					// additional track parameters (assume ntracks < 20)
  TStnTrackID*          fTrackID;
  TrackPar_t            fTrackPar[20];
					// histograms filled
  Hist_t                fHist;
  int                   fNTracks;
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
