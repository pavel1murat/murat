///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef murat_ana_TTrackCompModule_hh
#define murat_ana_TTrackCompModule_hh

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

#include "Stntuple/loop/TStnModule.hh"

#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TStnClusterBlock.hh"
#include "Stntuple/obj/TCalDataBlock.hh"
#include "Stntuple/obj/TStrawDataBlock.hh"
#include "Stntuple/obj/TGenpBlock.hh"
#include "Stntuple/obj/TSimpBlock.hh"
#include "Stntuple/obj/TVdetDataBlock.hh"

#include "Stntuple/base/TStnArrayI.hh"

#include "Stntuple/obj/TDiskCalorimeter.hh"

#include "Stntuple/alg/TStnTrackID.hh"
#include "Stntuple/alg/TEmuLogLH.hh"

#include "murat/ana/prob_dist.hh"

namespace mu2e { 
  class MVATools;
};

class TTrackCompModule: public TStnModule {
public:
#include "murat/ana/TrackPar_t.hh"
#include "murat/ana/SimPar_t.hh"
#include "murat/ana/HistBase_t.h"
//-----------------------------------------------------------------------------
//  histograms
//-----------------------------------------------------------------------------
  struct EventHist_t : public HistBase_t {

    TH1D*    fLumWt;		        // luminosity related MC weight
    TH1F*    fRv;			// MC truth information
    TH1F*    fZv;

    TH1F*    fPdgCode;
    TH1F*    fMomTargetEnd;
    TH1F*    fMomTrackerFront;
    TH1F*    fNshCE;

    TH1F*    fEleMom;
    TH1F*    fEleCosTh;
    TH1F*    fNTracks[2];
    TH1F*    fNshTot [2];
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
    TH1F*    fDp;                       // P(TrkPatRec)-P(CalPatRec)
    TH1F*    fInstLumi;                 // lumi
    TH1F*    fWeight;			// weight, need with statistics
  };

  struct TrackHist_t : public HistBase_t {

    TH1F*    fP[3];			// total momentum, 3 hists with different binning
    TH1F*    fP0;
    TH1F*    fP2;
    TH1D*    fPDio;
    TH1D*    fPlw;			// lumi-weighted momentum
    TH1D*    fPDiolw;			// lumi- and LO DIO-weighted momentum
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
    TH1F*    fChi2Dof;

    TH1F*    fNActive;
    TH1F*    fNaFract;
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

    TH2F*    fFConsVsNActive;
    TH1F*    fDaveTrkQual;
    TH1F*    fMVAOut;			// output of our MVA
    TH1F*    fDeltaMVA;			// DaveTrkQual-MVAOut[0]
  };
//-----------------------------------------------------------------------------
//  fTrackHist[  0]: all tracks
//  fTrackHist[100]: Set C tracks
//-----------------------------------------------------------------------------
  enum { kNEventHistSets   =  100 };
  enum { kNTrackHistSets   =  500 };
  enum { kNSimpHistSets    =  100 };

  struct Hist_t {
    EventHist_t*   fEvent  [kNEventHistSets];
    TrackHist_t*   fTrack  [kNTrackHistSets];
    //    SimpHist_t*    fSimp   [kNSimpHistSets];
  };


  struct Cut_t {
    double fXMin;
    double fXMax;
  };

  struct TmvaTrainingData_t {
    float    fP;
    float    fPMC;
    float    fTanDip;
    float    fNActive;
    float    fNaFract;
    float    fChi2Dof;
    float    fFitCons;
    float    fMomErr;
    float    fT0Err;
    float    fD0;
    float    fRMax;
    float    fNdaOverNa;
    float    fNzaOverNa;
    float    fNmaOverNa;
    float    fZ1;			// Z-coordinate of the first hit
    float    fWeight;
  };

  struct TmvaTrainingBranches_t {
    TBranch*  fP;
    TBranch*  fPMC;
    TBranch*  fTanDip;
    TBranch*  fNActive;
    TBranch*  fNaFract;
    TBranch*  fChi2Dof;
    TBranch*  fFitCons;
    TBranch*  fMomErr;
    TBranch*  fT0Err;
    TBranch*  fD0;
    TBranch*  fRMax;
    TBranch*  fNdaOverNa;
    TBranch*  fNzaOverNa;
    TBranch*  fNmaOverNa;
    TBranch*  fZ1;			// Z-coordinate of the first hit
    TBranch*  fWeight;			// for background only
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used
  TStnTrackBlock*   fTrackBlock[2];	// [0]: TrkPatRec tracks, [1]:all CalPatRec
  TStnClusterBlock* fClusterBlock;
  TGenpBlock*       fGenpBlock;
  TSimpBlock*       fSimpBlock;
  TVdetDataBlock*   fVdetBlock;
					
  TrackPar_t        fTrackPar[2][10];	// additional track parameters (assume ntracks < 10)
  SimPar_t          fSimPar;		// additional parameters of the simulated MC particle
  Hist_t            fHist;		// histograms filled

					// cut values
  double            fPtMin;

  Cut_t             fDebugCut[100];

  TGenParticle*     fParticle;		// electron or muon
  int               fPdgCode;		// determines which one
  int               fGeneratorCode;      

  TSimParticle*     fSimp;
  double            fEleE;		// electron energy

  int               fNTracks    [2];    // 0:TrkPatRec 1:CalPatRec
  int               fNGoodTracks[2];
  int               fNGenp;		// N(generated particles)

  TStnTrack*        fTrack;
  int               fFillDioHist;
					// [0]: SetC, [1-6]: TrkQual 0.1 ... 0.6
  int               fNID;
  TStnTrackID*      fTrackID[20];
  TStnTrackID*      fBestTrackID[2];
  int               fBestID[2];

  TEmuLogLH*        fLogLH;

  double            fMinT0;
  double            fLumWt;

  TStnCluster*      fCluster;
  int               fNClusters;
  double            fEClMax;
  double            fTClMax;
//-----------------------------------------------------------------------------
// TMVA training ntuples
//-----------------------------------------------------------------------------
  int                     fWriteTmvaTree;
  int                     fTmvaAlgorithm;   // write TMVS training tree for 0:TrkPatRec, 1:CalPatRec 

  TFile*                  fTmvaFile;

  TTree*                  fSigTree;
  //  TTree*                  fBgrTree;

  TmvaTrainingData_t      fTmvaData;
  TmvaTrainingBranches_t  fSigBranch;
  //  TmvaTrainingBranches_t  fBgrBranch;

  int                     fUseMVA;
  int                     fNMVA;	// number of MVA classifiers used for tracks of the same type
//-----------------------------------------------------------------------------
// MVA-based classifiers for TrkPatRec and CalPatRec tracks separately
//-----------------------------------------------------------------------------
  TString                 fTprWeightsFile;
  TString                 fCprWeightsFile;

  mu2e::MVATools*         fTprQualMva;
  mu2e::MVATools*         fCprQualMva;

  TString                 fTrkQualFile;
  prob_dist*              fTrackProb[2];// 0:TrkPatRec, 1:CalPatRec
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
  //  TStnTrackID*       GetTrackID   () { return fTrackID; }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void    SetPdgCode      (int Code ) { fPdgCode       = Code ; }
  void    SetGeneratorCode(int Code ) { fGeneratorCode = Code ; }
  void    SetDebugCut(int I, double XMin, double XMax) {
    fDebugCut[I].fXMin = XMin;
    fDebugCut[I].fXMax = XMax;
  }
  
  void    SetUseMVA        (int Flag) { fUseMVA        = Flag; }

  void    SetWriteTmvaTree (int Algo) { 
    fWriteTmvaTree = 1; 
    fTmvaAlgorithm = Algo; 
  }

  void    SetTprWeightsFile(const char* Fn) { if (Fn[0] != 0) fTprWeightsFile = Fn; }
  void    SetCprWeightsFile(const char* Fn) { if (Fn[0] != 0) fCprWeightsFile = Fn; }
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
  void    BookEventHistograms   (HistBase_t*   Hist, const char* Folder);
  void    BookTrackHistograms   (HistBase_t*   Hist, const char* Folder);

  void    FillEventHistograms    (HistBase_t*  Hist);
  void    FillTrackHistograms    (HistBase_t*  Hist, TStnTrack*    Trk , TrackPar_t* Tp);

  void    FillEfficiencyHistograms(TStnTrackBlock* TrackBlock, 
				   TStnTrackID*    TrackID   , 
				   TrackPar_t*     TPar      , 
				   int             HistSet   );

  int     FillTmvaTree();

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

  ClassDef(TTrackCompModule,0)
};

#endif
