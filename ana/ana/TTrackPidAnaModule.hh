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

#include "murat/ana/TAnaPart.hh"
#include "murat/ana/AnaDefs.hh"

class TTrackPidAnaModule: public TStnModule {
public:

  struct TrackPidPar_t {
    int     fLHRDedx;
  };
//-----------------------------------------------------------------------------
//  histograms
//-----------------------------------------------------------------------------
  struct PidHist_t {
    TH1F*    fLHEDedx;		       // 
    TH1F*    fLHMDedx;		       // 
    TH1F*    fLHRDedx;

    TH1F*    fNUsedOsH;
    TH1F*    fNUsedSsH;
    TH1F*    fNUsedOsD;

    TH1F*    fDrdsVadim;
    TH1F*    fDxdsVadim;
    TH1F*    fSumAvik;
    TH1F*    fMeanAvik;

    TH1F*    fSq2Avik;
    TH1F*    fMq2Avik;

    TH1F*    fDrdsOs;
    TH1F*    fDxdsOs;

    TH1F*    fDrdsSs;
    TH1F*    fDxdsSs;

    TH1F*    fSumAvikOs;
  };

  struct EventHist_t {
    TH1F*    fNTracks[8];		// N(reconstructed tracks)
  };

  struct TrackHist_t {
    TH1F*    fP;			// total momentum
  };

//-----------------------------------------------------------------------------
//  fTrackPidHist[ 0]: all the tracks
//  fTrackPidHist[ 1]: Set C tracks
//-----------------------------------------------------------------------------
  enum { kNEventHistSets    = 800 };
  enum { kNTrackHistSets    = 800 };
  enum { kNPidHistSets      = 800 };

  struct Hist_t {
    EventHist_t* fEvent[kNEventHistSets];
    TrackHist_t* fTrack[kNTrackHistSets];
    PidHist_t*   fPid  [kNPidHistSets  ];
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used

  TStnTrackBlock*      fTrackBlock   [8];
  TTrackStrawHitBlock* fTrackHitBlock[8];
  TStnPidBlock*        fPidBlock     [8];

  double               fLogProb      [8];  // 

					   // assume less than 20 particles
  TAnaPart             fPart[20];
					   // additional track PID parameters (for simplicity, assume Ntracks < 20)
  TrackPidPar_t        fTrackPidPar[20];
					   // histograms filled
  Hist_t               fHist;
					   // cut values
  double               fMinPt;
  int                  fNTracks[8];
  int                  fNID;
  TStnTrackID*         fTrackID[20];
  int                  fBestID;

  TEmuLogLH*           fLogLH;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TTrackPidAnaModule(const char* name="TrackPidAna", const char* title="TrackPidAna");
  ~TTrackPidAnaModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Hist_t*            GetHist         () { return &fHist;   }
  TEmuLogLH*         GetLogLH        () { return fLogLH;   }
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
  void    BookEventHistograms(EventHist_t*  Hist, const char* Folder);
  void    BookTrackHistograms(TrackHist_t*  Hist, const char* Folder);
  void    BookPidHistograms  (PidHist_t*    Hist, const char* Folder);

  void    FillEventHistograms(EventHist_t* Hist);
  void    FillTrackHistograms(TrackHist_t* Hist, TStnTrack*  Trk    );
  void    FillPidHistograms  (PidHist_t*   Hist, TStnPid*    Cluster);

  void    BookHistograms();
  void    FillHistograms();

  int     FindTrack(TAnaPart* Part, int Index);
  int     MakeParticles();


  void    Debug();
//-----------------------------------------------------------------------------
// tests : not used so far
//-----------------------------------------------------------------------------
  void    Test001();

  ClassDef(TTrackPidAnaModule,0)
};

#endif
