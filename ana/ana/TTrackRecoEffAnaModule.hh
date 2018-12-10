///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef murat_ana_TTrackRecoEffAnaModule_hh
#define murat_ana_TTrackRecoEffAnaModule_hh

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

#include "Stntuple/loop/TStnModule.hh"

#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TStnClusterBlock.hh"
#include "Stntuple/obj/TCalDataBlock.hh"
#include "Stntuple/obj/TStrawDataBlock.hh"
#include "Stntuple/obj/TGenpBlock.hh"
#include "Stntuple/obj/TSimpBlock.hh"
#include "Stntuple/obj/TVDetDataBlock.hh"

#include "Stntuple/base/TStnArrayI.hh"

#include "Stntuple/alg/TStnTrackID.hh"
#include "Stntuple/alg/TEmuLogLH.hh"

class TTrackRecoEffAnaModule: public TStnModule {
public:
//-----------------------------------------------------------------------------
//  histograms
//-----------------------------------------------------------------------------
  struct StrawHitHist_t {
    TH1F*    fPdgCode;		       //
    TH1F*    fGenCode;		       // generator code
    TH1F*    fEnergy;
    TH1F*    fTime;
    TH1F*    fDt;
    TH1F*    fMcMomentum;
  };

  struct EventHist_t {
    TH1F*    fNStrawHits[2];		// N straw hits in the event, different scales
    TH1F*    fNHitsSignal;
    TH1F*    fMomTF;			// 
    TH1F*    fNTracks;
    TH1F*    fFitCons;
    TH1F*    fT0Err;
    TH1F*    fNActive;
    TH1F*    fTanDipMC;
    TH1F*    fT0;
    TH1F*    fD0;
    TH1F*    fRMax;
    TH1F*    fTanDip;
    TH1F*    fP;
    TH1F*    fAlgMask;
    TH1F*    fClusterE;
  };

  struct RecoEffHist_t {
    TH1F*    fNTracks  [5];
    TH1F*    fFitCons  [5];
    TH1F*    fNActive  [5];
    TH1F*    fFitMomErr[5];
    TH1F*    fT0Err    [5];
    TH1F*    fT0       [5];
    TH1F*    fTanDip   [5];
    TH1F*    fD0       [5];
    TH1F*    fRMax     [5];
    TH1F*    fP        [5];
    TH1F*    fFailedBits;
    TH1F*    fPassed;
  };

  struct GenpHist_t {
    TH1F*    fP;
    TH1F*    fPdgCode[2];		// same distribution in different scale
    TH1F*    fGenID;
    TH1F*    fZ0;
    TH1F*    fT0;
    TH1F*    fR0;
    TH1F*    fCosTh;
  };

  struct TrackHist_t {
    TH1F*    fP;			// total momentum, 3 hists with different binning
    TH1F*    fP0;
    TH1F*    fP2;
  };
//-----------------------------------------------------------------------------
  enum { kNEventHistSets    = 100 };
  enum { kNStrawHitHistSets = 100 };
  enum { kNGenpHistSets     = 100 };
  enum { kNTrackHistSets    = 100 };
  enum { kNRecoEffHistSets  = 10  };

  struct Hist_t {
    EventHist_t*    fEvent    [kNEventHistSets];
    StrawHitHist_t* fStrawHit [kNStrawHitHistSets];
    GenpHist_t*     fGenp     [kNGenpHistSets];
    TrackHist_t*    fTrack    [kNTrackHistSets];
    RecoEffHist_t*  fRecoEff  [kNRecoEffHistSets];
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used
  TStnTrackBlock*       fTrackBlock;
  TStrawDataBlock*      fStrawDataBlock;
  TGenpBlock*           fGenpBlock;
  TVDetDataBlock*       fVDetDataBlock;
					// histograms filled
  Hist_t                fHist;

  TGenParticle*         fGenp;
  TSimParticle*         fSimp;

  int                   fNGenp;		// N(generated particles)
  int                   fNVDetHits;
  int                   fNStrawHits;
  int                   fNHitsSignal;
  int                   fPdgCode;
  int                   fGeneratorCode;

  int                   fNHitsTF;
  int                   fNHitsTB;
  float                 fMomTF;
  float                 fMomTB;
  float                 fPitchTF;

  int                   fNTracks;
  TStnTrack*            fTrack;
  float                 fNActive;
  float                 fFitCons;
  float                 fT0;
  float                 fT0Err;
  float                 fFitMomErr;
  float                 fD0;
  float                 fRMax;
  float                 fP;
  float                 fTanDipMC;
  float                 fTanDip;
  int                   fAlgMask;
  float                 fClusterE;
  int                   fRecoAlgFlag;
  int                   fTrkPatRecOnly;

  int                   fMinNMCHits;
  float                 fMinMCMomentum;
  float			fMinMCPitch;
  float			fMaxMCPitch;

  enum {
    kNTracksBit   = 0x1 <<  0,
    kFitConsBit   = 0x1 <<  1,
    kNActiveBit   = 0x1 <<  2,
    kFitMomErrBit = 0x1 <<  3,
    kT0ErrBit     = 0x1 <<  4,
    kT0Bit        = 0x1 <<  5,
    kTanDipBit    = 0x1 <<  6,
    kD0Bit        = 0x1 <<  7,
    kRMaxBit      = 0x1 <<  8,
    kPBit         = 0x1 <<  9
  };

  int                   fIDWord;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TTrackRecoEffAnaModule(const char* name="TrackRecoEffAna", const char* title="TrackRecoEffAna");
  ~TTrackRecoEffAnaModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Hist_t*            GetHist        () { return &fHist;        }
//-----------------------------------------------------------------------------
// setters
//-----------------------------------------------------------------------------
  void    SetPdgCode      (int Code) { fPdgCode       = Code; }
  void    SetGeneratorCode(int Code) { fGeneratorCode = Code; }
  void    SetTrkPatRecOnly(int Code) { fTrkPatRecOnly = Code; }

  void    SetMinNMCHits   (int   N ) { fMinNMCHits    = N;    }
  void    SetMinMCMomentum(float P ) { fMinMCMomentum = P;    }
  void    SetMinMCPitch   (float X ) { fMinMCPitch    = X;    }
  void    SetMaxMCPitch   (float X ) { fMaxMCPitch    = X;    }
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
  void    BookStrawHitHistograms (StrawHitHist_t* Hist, const char* Folder);
  void    BookEventHistograms    (EventHist_t*    Hist, const char* Folder);
  void    BookGenpHistograms     (GenpHist_t*     Hist, const char* Folder);
  void    BookTrackHistograms    (TrackHist_t*    Hist, const char* Folder);
  void    BookRecoEffHistograms  (RecoEffHist_t*  Hist, const char* Folder);

  void    FillStrawHitHistograms (StrawHitHist_t* Hist, TStrawHitData*  Hit);
  void    FillEventHistograms    (EventHist_t*    Hist);
  void    FillGenpHistograms     (GenpHist_t*     Hist, TGenParticle* Genp );
  void    FillTrackHistograms    (TrackHist_t*    Hist, TStnTrack*    Track);
  void    FillRecoEffHistograms  (RecoEffHist_t*  Hist);

  void    BookHistograms();
  void    FillHistograms();


  void    Debug();
//-----------------------------------------------------------------------------
// test
//-----------------------------------------------------------------------------
  void    Test001();

  ClassDef(TTrackRecoEffAnaModule,0)
};

#endif
