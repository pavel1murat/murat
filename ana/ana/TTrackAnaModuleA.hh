///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef murat_ana_TTrackAnaModuleA_hh
#define murat_ana_TTrackAnaModuleA_hh

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

#include "murat/ana/TTrackAnaModuleBase.hh"

#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TStnClusterBlock.hh"
#include "Stntuple/obj/TCalDataBlock.hh"
#include "Stntuple/obj/TStrawDataBlock.hh"
#include "Stntuple/obj/TVdetDataBlock.hh"
#include "Stntuple/obj/TGenpBlock.hh"
#include "Stntuple/obj/TSimpBlock.hh"

#include "Stntuple/base/TStnArrayI.hh"

#include "Stntuple/obj/TDiskCalorimeter.hh"

#include "Stntuple/alg/TStnTrackID.hh"
#include "Stntuple/alg/TEmuLogLH.hh"

class TTrackAnaModuleA: public TTrackAnaModuleBase {
public:
//-----------------------------------------------------------------------------
// track and sim particle additional parameters
//-----------------------------------------------------------------------------
#include "murat/ana/TrackPar_t.hh"
#include "murat/ana/SimPar_t.hh"

  enum { kNDisks = 2 };
//-----------------------------------------------------------------------------
//  histograms
//-----------------------------------------------------------------------------
  struct CaloHist_t: public HistData_t {
    TH1F*    fDiskID;		       // per crystal hit
    TH1F*    fEnergy  [kNDisks];
    TH1F*    fTime    [kNDisks];
    TH1F*    fNHits   [kNDisks];
    TH1F*    fRadius  [kNDisks];
    TH1F*    fRadiusWE[kNDisks];
    TH1F*    fE700    [kNDisks];
    TH1F*    fT700    [kNDisks];
    TH1F*    fN700    [kNDisks];
    TH1F*    fR700    [kNDisks];
    TH1F*    fRWE700  [kNDisks];
  };

  struct ClusterHist_t: public HistData_t {
    TH1F*    fDiskID;
    TH1F*    fEnergy;
    TH1F*    fT0;
    TH1F*    fRow;
    TH1F*    fCol;
    TH1F*    fX;
    TH1F*    fY;
    TH1F*    fZ;
    TH1F*    fR;
    TH1F*    fNCr0;			// all clustered
    TH1F*    fNCr1;			// above 1MeV
    TH1F*    fYMean;
    TH1F*    fZMean;
    TH1F*    fSigY;
    TH1F*    fSigZ;
    TH1F*    fSigR;
    TH1F*    fFrE1;
    TH1F*    fFrE2;
    TH1F*    fSigE1;
    TH1F*    fSigE2;
  };

  struct EventHist_t: public HistData_t  {
    TH1F*    fRv;			// MC truth information
    TH1F*    fZv;
    TH1F*    fEleMom;
    TH1D*    fDioMom;
    TH1F*    fEleCosTh;
    TH1F*    fNClusters;
    TH1F*    fNTracks;
    TH1F*    fNStrawHits[2];
    TH1F*    fNGoodSH;
    TH1F*    fMomTF;                    // signal particle momentuum @ Tracker Front
    TH1F*    fPitchTF;		        // SIM tan(pitch) @ Tracker Front
    TH1F*    fDtClT;
    TH1F*    fEMax;			// energy of the first reco cluster
    TH1F*    fDtClS;
    TH1F*    fSHTime;
    TH1F*    fNHyp;
    TH1F*    fBestHyp[2];		// [0]: by chi2, [1]: by fit consistency
    TH1F*    fNGenp;                    // N(particles in GENP block)

    TH1F*    fNCaloCrystalHits[kNDisks];
    TH2F*    fNCaloHitsVsDisk [kNDisks];
    TH2F*    fNCaloHitsVsRow  [kNDisks];
    TH2F*    fNCaloHitsVsCol  [kNDisks];
    // calorimeter hit histograms

    TH1F*    fETot        [kNDisks];            // total energy/event 
    TH2F*    fECrVsR      [kNDisks];            // total energy_per_crystal/event vs radius
    TH2F*    fNCrVsR      [kNDisks];            // total energy_per_crystal/event vs radius

    TH2F*    fNCrystalHitsVsR[kNDisks];            //
    TH2F*    fNHitCrystalsVsR[kNDisks];            //

    TH1F*    fNHitCrystalsTot;
    TH1F*    fECal;
    TH1F*    fECalOverEKin;
    TH1F*    fInstLumi;
  };

  struct TrackHist_t: public HistData_t  {
    TH1F*    fP[3];			// total momentum, 3 hists with different binning
    TH1F*    fP0;
    TH1F*    fP2;
    TH1F*    fPt;
    TH1D*    fPDio;                     // momentum dist weighted with the DIO weight
    TH1F*    fFitMomErr;
    TH1F*    fPFront;
    TH1F*    fDpFront;
    TH1F*    fXDpF;
    TH1F*    fDpFDio;
    TH1F*    fXDpFDio;
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
					// matching histograms
    TH1F*    fNClusters;
    TH1F*    fDiskID;
    TH1F*    fXCal;
    TH1F*    fYCal;
    TH1F*    fZCal;
    TH1F*    fXTrk;
    TH1F*    fYTrk;
    TH1F*    fZTrk;
    TH1F*    fRTrk;
    TH1F*    fDt;			// track-cluster residuals
    TH1F*    fChi2Match;
    TH1F*    fChi2XY;
    TH1F*    fChi2T;
    TH1F*    fDt_eMinus;
    TH1F*    fDt_ePlus;
    TH1F*    fDt_muMinus;
    TH1F*    fDt_muPlus;
    TH1F*    fDx;
    TH1F*    fDy;
    TH1F*    fDz;
    TH1F*    fDu;
    TH1F*    fDv;
    TH2F*    fDvVsDu;
    TH1F*    fPath;
    TH2F*    fDuVsPath;
    TH2F*    fDvVsPath;
    TH2F*    fDtVsPath;
    TH2F*    fDuVsTDip;
    TH2F*    fDvVsTDip;
    TH1F*    fZ1;
    TH1F*    fECl;
    TH1F*    fEClEKin;
    TH1F*    fEp;
    TH2F*    fEpVsPath;
    TH1F*    fEp_eMinus;
    TH1F*    fEp_ePlus;
    TH1F*    fEp_muMinus;
    TH1F*    fEp_muPlus;
    TH2F*    fNHVsStation;
    TH2F*    fNHVsNSt;

    TH1F*    fRSlope;
    TH1F*    fXSlope;
					// likelihoods
    TH2F*    fEpVsDt;
    TH1F*    fEleLogLHCal;
    TH1F*    fMuoLogLHCal;
    TH1F*    fLogLHRCal;
    TH1F*    fLogLHRDeDx;
    TH1F*    fLogLHRXs;
    TH1F*    fLogLHRTrk;
    TH1F*    fLogLHR;
					// MC truth
    TH1F*    fPdgCode;	                // PDG code of the particle produced most hits
    TH1F*    fFrGH;			// fraction of hits produced by the particle

    TH2F*    fNEPlVsNHPl;
    TH2F*    fNDPlVsNHPl;
    TH2F*    fChi2dVsNDPl;
    TH2F*    fDpFVsNDPl;

    TH1F*    fFrE1;
    TH1F*    fFrE2;

    TH1F*    fSinTC;			// sin(track-cluster angle)
    TH1F*    fDrTC;                     // deltaR(cluster-track)
    TH1F*    fSInt;                     // calculated interaction length
    TH1F*    fDaveTrkQual;		// 
    TH1F*    fNMcStrawHits;             // N(straw hits) produced in the tracker by the MC particle
  };

  struct GenpHist_t: public HistData_t  {
    TH1F*    fPdgCode[2];		// same distribution in different scale
    TH1F*    fGenID;			// 
    TH1F*    fZ0;			// 
    TH1F*    fT0;			// 
    TH1F*    fR0;			// 
    TH1F*    fP;			// 
    TH1F*    fCosTh;			// 
  };
					// histograms for the simulated CE
  struct SimpHist_t: public HistData_t  {
    TH1F*    fPdgCode;
    TH1F*    fMomTargetEnd;
    TH1F*    fMomTrackerFront;
    TH1F*    fNStrawHits;
  };

  struct TrackEffHist_t: public HistData_t  {
    TH1F*    fPtMc;			// denominator
    TH1F*    fPtReco;			// numerator
  };
//-----------------------------------------------------------------------------
  enum { kNEventHistSets   = 100 };
  enum { kNTrackHistSets   = 400 };
  enum { kNClusterHistSets = 100 };
  enum { kNCaloHistSets    = 100 };
  enum { kNGenpHistSets    = 100 };
  enum { kNSimpHistSets    = 100 };

  struct Hist_t {
    TH1F*          fCrystalR[2];	          // crystal radius
    EventHist_t*   fEvent  [kNEventHistSets  ];
    TrackHist_t*   fTrack  [kNTrackHistSets  ];
    ClusterHist_t* fCluster[kNClusterHistSets];
    CaloHist_t*    fCalo   [kNCaloHistSets   ];
    GenpHist_t*    fGenp   [kNGenpHistSets   ];
    SimpHist_t*    fSimp   [kNSimpHistSets   ];
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used
  TStnTrackBlock*   fTrackBlock;
  TStnClusterBlock* fClusterBlock;
  TCalDataBlock*    fCalDataBlock;
  TStrawDataBlock*  fStrawDataBlock;
  TVdetDataBlock*   fVdetBlock;
  TGenpBlock*       fGenpBlock;
  TSimpBlock*       fSimpBlock;
					// histograms filled
  Hist_t            fHist;
					// cut values
  double            fPtMin;

  TGenParticle*     fParticle;		// electron or muon
  int               fPdgCode;		// determines which one
  int               fGeneratorCode;      

  TSimParticle*     fSimp;
  double            fEleE;		// electron energy

  int               fNClusters;
  int               fNTracks[10];
  int               fNGoodTracks;
  int               fNCalPatRec;
  int               fNMatchedTracks;
  int               fNStrawHits;
  int               fNCalHits;
  int               fNGenp;		// N(generated particles)

  int               fNHyp;
  int               fBestHyp[10];
  int               fFillDioHist;

  double            fMinT0;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TTrackAnaModuleA(const char* name="TrackAnaA", const char* title="TrackAnaA");
  ~TTrackAnaModuleA();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Hist_t*            GetHist        () { return &fHist;        }
  TStnTrackBlock*    GetTrackBlock  () { return fTrackBlock;   }
  TStnClusterBlock*  GetClusterBlock() { return fClusterBlock; }

  TStnTrackID*       GetTrackID(int I) { return fTrackID[I];   }
  TEmuLogLH*         GetLogLH       () { return fLogLH;        }
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
  void    BookCaloHistograms    (HistData_t* Hist, const char* Folder);
  void    BookClusterHistograms (HistData_t* Hist, const char* Folder);
  void    BookGenpHistograms    (HistData_t* Hist, const char* Folder);
  void    BookEventHistograms   (HistData_t* Hist, const char* Folder);
  void    BookSimpHistograms    (HistData_t* Hist, const char* Folder);
  void    BookTrackHistograms   (HistData_t* Hist, const char* Folder);

  void    FillEventHistograms    (HistData_t* Hist);
  void    FillCaloHistograms     (HistData_t* Hist, TStnCrystal*  Crystal);
  void    FillClusterHistograms  (HistData_t* Hist, TStnCluster*  Cluster);
  void    FillGenpHistograms     (HistData_t* Hist, TGenParticle* Genp   );
  void    FillSimpHistograms     (HistData_t* Hist, TSimParticle* Simp   );
  void    FillTrackHistograms    (HistData_t* Hist, TStnTrack*    Trk    );

  void    BookHistograms();
  void    FillHistograms();

  int     InitTrackPar(TStnTrackBlock*   TrackBlock  , 
		       TStnClusterBlock* ClusterBlock, 
		       TrackPar_t*       TrackPar    );

  void    Debug();
//-----------------------------------------------------------------------------
// test
//-----------------------------------------------------------------------------
  void    Test001();

  ClassDef(TTrackAnaModuleA,0)
};

#endif