///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef murat_ana_TTrackCompModule_hh
#define murat_ana_TTrackCompModule_hh

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

#include "Stntuple/base/TStnArrayI.hh"

#include "Stntuple/obj/TDiskCalorimeter.hh"

#include "Stntuple/alg/TStnTrackID.hh"
#include "Stntuple/alg/TEmuLogLH.hh"

class TTrackCompModule: public TStnModule {
public:

  struct TrackPar_t {
    int     fNHPl;
    int     fNEPl;
    int     fNDPl;

    int     fIDWord_2025;		// ID word for 20 <= N(active) < 25 tracks
    int     fIDWord_30;	                // ID word for 30 <= N(active)

    float   fDpF ;                      // tracker-only resolution
    float   fDp0 ;
    float   fDp2 ;
    float   fDpFSt;
    double  fDioWt;

    double  fEcl;
    double  fEp;
    double  fDx;
    double  fDy;
    double  fDz;
    double  fDt;
    double  fDu;			// rotated residuals
    double  fDv;
    double  fChi2Match;
    double  fChi2XY;
    double  fChi2T;
    double  fPath;

    double  fSinTC;			// angle between the cluster and the track
  };
//-----------------------------------------------------------------------------
//  histograms
//-----------------------------------------------------------------------------
  struct EventHist_t {
    TH1F*    fRv;			// MC truth information
    TH1F*    fZv;
    TH1F*    fEleMom;
    TH1F*    fEleCosTh;
    TH1F*    fNTracks[2];
    TH1F*    fNStrawHits[2];
    TH1F*    fNGoodSH;
    TH1F*    fDtClT;
    TH1F*    fDtClS;
    TH1F*    fSHTime;
    TH1F*    fNHyp;
    TH1F*    fBestHyp[2];		// [0]: by chi2, [1]: by fit consistency
    TH1F*    fNGenp;                    // N(particles in GENP block)
    TH1F*    fNClusters;
    TH1F*    fEClMax;			// energy of the first (highest) reconstructed cluster
    TH1F*    fTClMax;			// time   of the first (highest) reconstructed cluster
  };

  struct TrackHist_t {
    TH1F*    fP[3];			// total momentum, 3 hists with different binning
    TH1F*    fP0;
    TH1F*    fP2;
    TH1D*    fPDio;
    TH1F*    fPt;
    TH1F*    fFitMomErr;
    TH1F*    fPFront;
    TH1F*    fDpFront;
    TH1F*    fXDpF;                     // DpF/MomErr
    TH1F*    fDpFDio;
    TH1F*    fDpFront0;
    TH1F*    fDpFront2;
    TH2F*    fDpFVsZ1;
    TH1F*    fPStOut;
    TH1F*    fDpFSt;			// P(TT_Hollow) - P(ST_Out)
    TH1F*    fCosTh;
    TH1F*    fChi2;
    TH1F*    fNDof;
    TH1F*    fChi2Dof;
    TH1F*    fChi2DofC;
    TH1F*    fNActive;
    TH1F*    fNWrong;
    TH1F*    fT0;
    TH1F*    fT0Err;
    TH1F*    fQ;
    TH1F*    fFitCons[2];		// fit consistency (0 to 1)
    TH1F*    fD0;
    TH1F*    fZ0;
    TH1F*    fTanDip;
    TH1F*    fResid;
    TH1F*    fAlgMask;
    TH2F*    fFConsVsNActive;
  };

  struct SimpHist_t {
    TH1F*    fPdgCode;
    TH1F*    fMomTargetEnd;
    TH1F*    fMomTrackerFront;
    TH1F*    fNStrawHits;
  };
//-----------------------------------------------------------------------------
//  fTrackHist[  0]: all tracks
//  fTrackHist[100]: Set C tracks
//-----------------------------------------------------------------------------
  enum { kNEventHistSets   = 100 };
  enum { kNTrackHistSets   = 400 };
  enum { kNSimpHistSets    = 100 };

  struct Hist_t {
    EventHist_t*   fEvent  [kNEventHistSets];
    TrackHist_t*   fTrack  [kNTrackHistSets];
    SimpHist_t*    fSimp   [kNSimpHistSets];
  };


  struct Cut_t {
    double fXMin;
    double fXMax;
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used
  TStnTrackBlock*   fTrackBlock[2];	// [0]: TrkPatRec, [1]:CalPatRec
  TStnClusterBlock* fClusterBlock;
  TGenpBlock*       fGenpBlock;
  TSimpBlock*       fSimpBlock;
					// additional track parameters (assume ntracks < 10)
  TrackPar_t        fTrackPar[2][10];
					// histograms filled
  Hist_t            fHist;
					// cut values
  double            fPtMin;

  Cut_t             fDebugCut[100];

  TGenParticle*     fParticle;		// electron or muon
  int               fPdgCode;		// determines which one
  int               fGeneratorCode;      

  TSimParticle*     fSimp;
  double            fEleE;		// electron energy

  int               fNTracks[2];        // 0:TrkPatRec 1:CalPatRec
  int               fNGoodTracks[2];
  int               fNGenp;		// N(generated particles)

					// fTrackNumber[i]: track number, 
					// corresponding to OBSP particle #i
					// or -1
  TStnArrayI        fTrackNumber;

  TStnTrack*        fTrack;
  TStnTrackID*      fTrackID;
  TStnTrackID*      fTrackID_2025;
  TStnTrackID*      fTrackID_30;

  double            fMinT0;

  TStnCluster*      fCluster;
  int               fNClusters;
  double            fEClMax;
  double            fTClMax;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TTrackCompModule(const char* name="TrackComp", const char* title="TrackComp");
  ~TTrackCompModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Hist_t*            GetHist      () { return &fHist;      }
  TStnTrackID*       GetTrackID   () { return fTrackID; }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void               SetPdgCode      (int Code ) { fPdgCode       = Code ; }
  void               SetGeneratorCode(int Code ) { fGeneratorCode = Code ; }
  void               SetDebugCut(int I, double XMin, double XMax) {
    fDebugCut[I].fXMin = XMin;
    fDebugCut[I].fXMax = XMax;
  }
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
  void    BookEventHistograms   (EventHist_t*   Hist, const char* Folder);
  void    BookSimpHistograms    (SimpHist_t*    Hist, const char* Folder);
  void    BookTrackHistograms   (TrackHist_t*   Hist, const char* Folder);

  void    FillEventHistograms    (EventHist_t*  Hist);
  void    FillSimpHistograms     (SimpHist_t*   Hist, TSimParticle* Simp);
  void    FillTrackHistograms    (TrackHist_t*  Hist, TStnTrack*    Trk , TrackPar_t* Tp);

  void    BookHistograms();
  void    FillHistograms();


  void    Debug();
//-----------------------------------------------------------------------------
// test
//-----------------------------------------------------------------------------
  void    Test001();

  ClassDef(TTrackCompModule,0)
};

#endif
