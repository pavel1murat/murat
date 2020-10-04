///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef Stntuple_ana_TRMCAnaModule_hh
#define Stntuple_ana_TRMCAnaModule_hh

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

#include "Stntuple/loop/TStnModule.hh"

#include "Stntuple/obj/TStnTimeClusterBlock.hh"
#include "Stntuple/obj/TStnHelixBlock.hh"
#include "Stntuple/obj/TStnTrackSeedBlock.hh"
#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TStnClusterBlock.hh"
#include "Stntuple/obj/TCalDataBlock.hh"
#include "Stntuple/obj/TStrawDataBlock.hh"
#include "Stntuple/obj/TGenpBlock.hh"
#include "Stntuple/obj/TSimpBlock.hh"

#include "Stntuple/base/TStnArrayI.hh"

#include "Stntuple/geom/TStnCrystal.hh"
#include "Stntuple/alg/TStnTrackID.hh"
#include "Stntuple/alg/TEmuLogLH.hh"

class TRMCAnaModule: public TStnModule {
public:

#include "murat/ana/TrackPar_t.hh"

  struct EventHist_t {
    TH1F*    fRv;			// MC truth information
    TH1F*    fZv;
    TH1F*    fEleMom;
    TH1D*    fDioMom;
    TH1F*    fEleCosTh;
    TH1F*    fNTracks;
    TH1F*    fNStrawHits[2];
    TH1F*    fNGoodSH;
    TH1F*    fDtClT;
    TH1F*    fEMax;			// energy of the first reco cluster
    TH1F*    fDtClS;
    TH1F*    fSHTime;
    TH1F*    fNHyp;
    TH1F*    fBestHyp[2];		// [0]: by chi2, [1]: by fit consistency
    TH1F*    fNGenp;                    // N(particles in GENP block)
  };

  struct TrackHist_t {

    TH1F*    fP[3];			// total momentum, 3 hists with different binning
    TH1F*    fP0;
    TH1F*    fP2;
    TH1F*    fPt;
    TH1F*    fFitMomErr;
    TH1F*    fPFront;
    TH1F*    fDpFront;
    TH1F*    fXDpF;                     // DpF/MomErr
    TH1F*    fDpFront0;
    TH1F*    fDpFront2;
    TH2F*    fDpFVsZ1;
    TH1F*    fPStOut;
    TH1F*    fDpFSt;			// P(TT_Hollow) - P(ST_Out)
    TH1F*    fCosTh;
    TH1F*    fChi2;
    TH1F*    fChi2Dof;

    TH1F*    fNActive;
    TH1F*    fNaFract;
    TH1F*    fDNa;
    TH1F*    fNWrong;			// MC-only histogram: N(hits) with wrong drift signs
    TH1F*    fNDoublets;
    TH1F*    fNadOverNd;		// fraction of doublets with all hits active
    TH1F*    fNSSD;
    TH1F*    fNOSD;
    TH1F*    fNdOverNa;

    TH1F*    fNssdOverNa;
    TH1F*    fNosdOverNa;
    TH1F*    fNZeroAmb;
    TH1F*    fNzaOverNa;
    TH1F*    fNMatActive;
    TH1F*    fNmaOverNa;
    TH1F*    fNBend;

    TH1F*    fT0;
    TH1F*    fT0Err;
    TH1F*    fQ;
    TH1F*    fFitCons[2];		// fit consistency (0 to 1)
    TH1F*    fD0;
    TH1F*    fZ0;
    TH1F*    fTanDip;
    TH1F*    fRMax;
    TH1F*    fDtZ0;			// MC truth: T0-T(MC TMid)
    TH1F*    fXtZ0;                     // pull(dt) at Z=0

    TH1F*    fResid;
    TH1F*    fAlgMask;
					// matching
    TH1F*    fChi2Tcm;
    TH1F*    fChi2XY;
    TH1F*    fChi2T;

    TH1F*    fDt;			// track-cluster residuals
    TH1F*    fDx;
    TH1F*    fDy;
    TH1F*    fDz;
    TH1F*    fDu;
    TH1F*    fDv;
    TH1F*    fPath;

    TH1F*    fECl;
    TH1F*    fEClEKin;
    TH1F*    fEp;
    TH1F*    fDrDzCal;
    TH1F*    fDtClZ0;                   // T(cluster back at Z0)-T_true(Z0)
    TH2F*    fDtClZ0VsECl;              // 
    TH2F*    fDtClZ0VsP;              // 

    TH2F*    fFConsVsNActive;
    TH1F*    fDaveTrkQual;
    TH1F*    fMVAOut;			// output of our MVA
    TH1F*    fDeltaMVA;			// DaveTrkQual-MVAOut[0]
  };

  struct GenpHist_t {
    TH1F*    fPdgCode[2];		// same distribution in different scale
    TH1F*    fGenID;			// 
    TH1F*    fZ0;			// 
    TH1F*    fT0;			// 
    TH1F*    fR0;			// 
    TH1F*    fP;			// 
    TH1F*    fCosTh;			// 
  };
					// histograms for the simulated CE
  struct SimpHist_t {
    TH1F*    fPdgCode;
    TH1F*    fMomTargetEnd;
    TH1F*    fMomTrackerFront;
    TH1F*    fNStrawHits;
  };

  enum { kNEventHistSets   = 100 };
  enum { kNTrackHistSets   = 400 };
  enum { kNGenpHistSets    = 100 };
  enum { kNSimpHistSets    = 100 };

  struct Hist_t {
    EventHist_t*     fEvent     [kNEventHistSets];
    TrackHist_t*     fTrack     [kNTrackHistSets];
    GenpHist_t*      fGenp      [kNGenpHistSets];
    SimpHist_t*      fSimp      [kNSimpHistSets];
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used

  TGenpBlock*       fGenpBlock;
  TSimpBlock*       fSimpBlock;
					// PAR and DAR
  TStnTrackBlock*   fTrackBlock[2];
  TrackPar_t        fTrackPar[2][10];	// additional track parameters (assume ntracks < 10)

					// histograms filled
  Hist_t            fHist;
					// cut values
  double            fPtMin;

  TGenParticle*     fParticle;		// electron or muon
  int               fPdgCode;		// determines which one
  int               fGeneratorCode;      

  TSimParticle*     fSimp;
  double            fEleE;		// electron energy

  int               fNTrackSeeds[5];
  int               fNTracks[10];
  int               fNGoodTracks;
  int               fNMatchedTracks;
  int               fNGenp;		// N(generated particles)

  int               fNHyp;
  int               fBestHyp[10];
  int               fFillDioHist;
					// fTrackNumber[i]: track number, 
					// corresponding to OBSP particle #i
					// or -1
  TStnArrayI        fTrackNumber;

  TStnTrack*        fTrack;

  TStnTrackID*      fTrackID;
  TStnTrackID*      fBestTrackID[2];
  TEmuLogLH*        fLogLH;

  double            fMinT0;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TRMCAnaModule(const char* name="StnTrackAna", const char* title="StnTrackAna");
  ~TRMCAnaModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Hist_t*            GetHist        () { return &fHist;        }
  TStnTrackBlock*    GetTrackBlock  () { return fTrackBlock;   }

  TStnTrackID*       GetTrackID     () { return fTrackID; }
  TEmuLogLH*         GetLogLH       () { return fLogLH; }
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  void               SetFillDioHist  (int YesNo) { fFillDioHist   = YesNo; }
  void               SetPdgCode      (int Code ) { fPdgCode       = Code ; }
  void               SetGeneratorCode(int Code ) { fGeneratorCode = Code ; }
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
  void    BookGenpHistograms    (GenpBase_t*  Hist, const char* Folder);
  void    BookEventHistograms   (EventBase_t* Hist, const char* Folder);
  void    BookSimpHistograms    (SimpBase_t*  Hist, const char* Folder);
  void    BookTrackHistograms   (TrackBase_t* Hist, const char* Folder);

  void    FillEventHistograms   (GenpBase_t*  Hist);
  void    FillGenpHistograms    (EventBase_t* Hist, TGenParticle* Genp);
  void    FillSimpHistograms    (SimpBase_t*  Hist, TSimParticle* Simp);
  void    FillTrackHistograms   (TrackBase_t* Hist, TStnTrack*    Trk , TrackPar_t* Tp, double Weight);

  void    BookHistograms();
  void    FillHistograms();


  void    Debug();
//-----------------------------------------------------------------------------
// test
//-----------------------------------------------------------------------------
  void    Test001();

  ClassDef(TRMCAnaModule,0)
};

#endif
