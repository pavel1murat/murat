///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef murat_ana_TTrackPidAnaModule_hh
#define murat_ana_TTrackPidAnaModule_hh

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

#include "Stntuple/loop/TStnModule.hh"

#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TStnPidBlock.hh"
#include "Stntuple/obj/TStnPid.hh"

#include "Stntuple/base/TStnArrayI.hh"

#include "Stntuple/alg/TStnTrackID.hh"
#include "Stntuple/alg/TEmuLogLH.hh"

class TTrackPidAnaModule: public TStnModule {
public:

  struct TrackPidPar_t {
    int     fLHRDedx;
  };
//-----------------------------------------------------------------------------
//  histograms
//-----------------------------------------------------------------------------
  struct TrackPidHist_t {
    TH1F*    fLHEDedx;		       // 
    TH1F*    fLHMDedx;		       // 
    TH1F*    fLHRDedx;		       // 
    TH1F*    fDrdsVadimEle;
    TH1F*    fDrdsVadimMuo;
    TH1F*    fDxdsVadimEle;
    TH1F*    fDxdsVadimMuo;
    TH1F*    fSumAvikEle;
    TH1F*    fSumAvikMuo;

    TH1F*    fSq2AvikEle;
    TH1F*    fSq2AvikMuo;

    TH1F*    fDrdsOsEle;
    TH1F*    fDrdsOsMuo;
    TH1F*    fDxdsOsEle;
    TH1F*    fDxdsOsMuo;

    TH1F*    fDrdsSsEle;
    TH1F*    fDrdsSsMuo;
    TH1F*    fDxdsSsEle;
    TH1F*    fDxdsSsMuo;

    TH1F*    fSumAvikOsEle;
    TH1F*    fSumAvikOsMuo;
    TH1F*    fNUsedOsEle;
    TH1F*    fNUsedOsMuo;
    TH1F*    fLHROs;
  };

  struct EventHist_t {
    TH1F*    fNTracks[2];		// N(reconstructed tracks)
    TH2F*    fNDmmVsNDem;
    TH1F*    fNDem;
    TH1F*    fNDmm;
  };

  struct TrackHist_t {
    TH1F*    fP;			// total momentum
  };

//-----------------------------------------------------------------------------
//  fTrackPidHist[ 0]: all the tracks
//  fTrackPidHist[ 1]: Set C tracks
//-----------------------------------------------------------------------------
  enum { kNEventHistSets    = 100 };
  enum { kNTrackHistSets    = 100 };
  enum { kNTrackPidHistSets =  10 };

  struct Hist_t {
    EventHist_t*    fEvent   [kNEventHistSets   ];
    TrackHist_t*    fTrack   [kNTrackHistSets   ];
    TrackPidHist_t* fTrackPid[kNTrackPidHistSets];
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used
  TStnTrackBlock*   fTrackBlockDem;
  TStnTrackBlock*   fTrackBlockDmm;
  TStnPidBlock*     fTrackPidBlock;
					// additional track PID parameters (for simplicity, assume Ntracks < 20)
  TrackPidPar_t     fTrackPidPar[20];
					// histograms filled
  Hist_t            fHist;
					// cut values
  double            fMinPt;
  double            fMinT0;

  int               fNTracks[10];
  int               fNTracksDem;
  int               fNTracksDmm;

  int               fNPid;

  TStnTrack*        fTrack;

  TStnTrackID*      fTrackID;

  TEmuLogLH*        fLogLH;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TTrackPidAnaModule(const char* name="TrackPidAna", const char* title="TrackPidAna");
  ~TTrackPidAnaModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Hist_t*            GetHist         () { return &fHist;         }
  TStnTrackBlock*    GetTrackBlockDem() { return fTrackBlockDem; }
  TStnTrackBlock*    GetTrackBlockDmm() { return fTrackBlockDmm; }
  TStnPidBlock*      GetTrackPidBlock() { return fTrackPidBlock; }

  TStnTrackID*       GetTrackID      () { return fTrackID; }
  TEmuLogLH*         GetLogLH        () { return fLogLH; }
//-----------------------------------------------------------------------------
// accessors
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
  void    BookEventHistograms   (EventHist_t*    Hist, const char* Folder);
  void    BookTrackHistograms   (TrackHist_t*    Hist, const char* Folder);
  void    BookTrackPidHistograms(TrackPidHist_t* Hist, const char* Folder);

  void    FillEventHistograms    (EventHist_t*    Hist);
  void    FillTrackHistograms    (TrackHist_t*    Hist, TStnTrack*  Trk    );
  void    FillTrackPidHistograms (TrackPidHist_t* Hist, TStnPid*    Cluster);

  void    BookHistograms();
  void    FillHistograms();


  void    Debug();
//-----------------------------------------------------------------------------
// tests : not used so far
//-----------------------------------------------------------------------------
  void    Test001();

  ClassDef(TTrackPidAnaModule,0)
};

#endif
