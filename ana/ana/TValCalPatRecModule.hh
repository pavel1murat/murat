///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef murat_ana_TValCalPatRecModule_hh
#define murat_ana_TValCalPatRecModule_hh

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

#include "Stntuple/loop/TStnModule.hh"

#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TStnClusterBlock.hh"

#include "Stntuple/base/TStnArrayI.hh"

#include "Stntuple/obj/TDiskCalorimeter.hh"

#include "Stntuple/alg/TStnTrackID.hh"
#include "Stntuple/alg/TEmuLogLH.hh"

class TValCalPatRecModule: public TStnModule {
//-----------------------------------------------------------------------------
//  histograms
//-----------------------------------------------------------------------------
public:

  struct TrackPar_t {
    int     fNHPl;
    int     fNEPl;
    int     fNDPl;
    float   fDpF ;    // tracker-only resolution
    float   fDp0 ;
    float   fDp2 ;
    float   fDpFSt;
    double  fDioWt;
  };

  struct EventHist_t {
    TH1F*    fRv;			// MC truth information
    TH1F*    fZv;
    TH1F*    fEleCosTh;
    TH1F*    fNClusters;
    TH2F*    fNT2VsNT1;
    TH2F*    fNGT2VsNGT1;
  };

  struct TrackHist_t {
    TH1F*    fP;			// total momentum, 3 hists with different binning
    TH1F*    fPt;
    TH1F*    fPFront;
    TH1F*    fDpFront;
    TH1F*    fPStOut;
    TH1F*    fDpFSt;			// P(TT_Hollow) - P(ST_Out)
    TH1F*    fCosTh;
    TH1F*    fChi2;
    TH1F*    fNDof;
    TH1F*    fChi2Dof;
    TH1F*    fNActive;
    TH1F*    fT0;
    TH1F*    fQ;
    TH1F*    fFitCons[2];		// fit consistency (0 to 1)
    TH1F*    fD0;
    TH1F*    fZ0;
    TH1F*    fTanDip;
    TH1F*    fResid;
					// matching histograms
    TH1F*    fNClusters;
    TH1F*    fVaneID;
    TH1F*    fXCal;
    TH1F*    fYCal;
    TH1F*    fZCal;
    TH1F*    fXTrk;
    TH1F*    fYTrk;
    TH1F*    fZTrk;
    TH1F*    fRTrk;
    TH1F*    fDt;			// track-cluster residuals
    TH1F*    fDt_eMinus;
    TH1F*    fDt_ePlus;
    TH1F*    fDt_muMinus;
    TH1F*    fDt_muPlus;
    TH1F*    fDx;
    TH1F*    fDy;
    TH1F*    fDz;
    TH1F*    fEp;
    TH1F*    fEp_eMinus;
    TH1F*    fEp_ePlus;
    TH1F*    fEp_muMinus;
    TH1F*    fEp_muPlus;
    TH2F*    fNHVsStation;
    TH2F*    fNHVsNSt;
					// MC truth
    TH1F*    fPdgCode;	// PDG code of the particle produced most hits
    TH1F*    fFrGH;			// fraction of hits produced by the particle

    TH2F*    fNEPlVsNHPl;
    TH2F*    fNDPlVsNHPl;
    TH2F*    fChi2dVsNDPl;
    TH2F*    fDpFVsNDPl;

  };

  enum { kNEventHistSets   = 100 };
  enum { kNTrackHistSets   = 400 };

  struct Hist_t {
    EventHist_t*   fEvent  [kNEventHistSets];
    TrackHist_t*   fTrack1 [kNTrackHistSets];
    TrackHist_t*   fTrack2 [kNTrackHistSets];
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used
  TStnTrackBlock*   fTrackBlock[2];
  TStnClusterBlock* fClusterBlock;
					// additional track parameters (assume ntracks < 20)
  TrackPar_t        fTrackPar[2][20];
					// histograms filled
  Hist_t            fHist;

  int               fNClusters;
  int               fNTracks[10];
  int               fNGoodTracks[10];
  int               fNMatchedTracks[10];

  TStnTrack*        fTrack;
  TStnCluster*      fCluster;

  TStnTrackID*      fTrackID;

  double            fMinT0;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TValCalPatRecModule(const char* name="ValCalPatRec", const char* title="ValCalPatRec");
  ~TValCalPatRecModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Hist_t*            GetHist        () { return &fHist;        }
  TStnClusterBlock*  GetClusterBlock() { return fClusterBlock; }

  TStnTrackID*       GetTrackID     () { return fTrackID; }
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
  void    BookEventHistograms   (EventHist_t*  Hist, const char* Folder);
  void    BookTrackHistograms   (TrackHist_t*  Hist, const char* Folder);

  void    FillEventHistograms    (EventHist_t* Hist);
  void    FillTrackHistograms    (TrackHist_t* Hist, TStnTrack*  Trk, TrackPar_t* Tp);

  void    BookHistograms();
  void    FillHistograms();
  void    Debug();
//-----------------------------------------------------------------------------
// test
//-----------------------------------------------------------------------------
  void    Test001();

  ClassDef(TValCalPatRecModule,0)
};

#endif
