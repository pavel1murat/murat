//
// Read the tracks added to the event by the track finding and fitting code.
//
// $Id: TCalm002_module.hh,v 1.19 2015/01/06 18:36:46 murat Exp $
// $Author: murat $
// $Date: 2015/01/06 18:36:46 $
//
// Contact person Pavel Murat
//
#ifndef murat_inc_TCalm002_module_hh
#define murat_inc_TCalm002_module_hh

// Mu2e includes.
// #include "CLHEP/Geometry/HepPoint.h"
#include "BTrk/BbrGeom/HepPoint.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Vector.h"

#include "BTrk/TrkBase/TrkHelixUtils.hh"
#include "BTrk/TrkBase/HelixParams.hh"
#include "BTrk/KalmanTrack/KalHit.hh"

#include "RecoDataProducts/inc/KalRepPtrCollection.hh"
// storable objects (data products)
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"

#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"


// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

// C++ includes.
#include <iostream>

#include "TString.h"
#include "TFolder.h"
#include "TFile.h"

#include "Stntuple/obj/TStnClusterBlock.hh"
#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TCalDataBlock.hh"
#include "Stntuple/obj/TStrawDataBlock.hh"

#include "Stntuple/geom/TDiskCalorimeter.hh"

#include "Stntuple/mod/THistModule.hh"

#include "Stntuple/alg/TEmuLogLH.hh"

class TStnTrackID;
class TStnCrystal;

namespace mu2e {

  class TCalm002 : public THistModule {
  public:

    struct CaloHist_t {
      TH1F*    fVaneID;		       // per crystal hit
      TH1F*    fEnergy  [4];
      TH1F*    fTime    [4];
      TH1F*    fNHits   [4];
      TH1F*    fRadius  [4];
      TH1F*    fRadiusWE[4];
      TH1F*    fE700    [4];
      TH1F*    fT700    [4];
      TH1F*    fN700    [4];
      TH1F*    fR700    [4];
      TH1F*    fRWE700  [4];
    };

    struct ClusterHist_t {
      TH1F*    fVaneID;
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

    struct TrackHist_t {
      TH1F*    fP[3];			// total momentum, 3 hists with different binning
      TH1F*    fPt;
      TH1F*    fPFront;
      TH1F*    fDpFront;
      TH1F*    fPStOut;
      TH1F*    fDpFSt;			// P(TT_Hollow) - P(ST_Out)
      TH1F*    fCosTh;
      TH1F*    fChi2;
      TH1F*    fNDof;
      TH1F*    fChi2Dof;
      TH1F*    fChi2DofC;
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

      TH1F*    fRSlope;
      TH1F*    fXSlope;
					// likelihoods
      TH1F*    fEleLogLHCal;
      TH1F*    fMuoLogLHCal;
      TH1F*    fLogLHRCal;
      TH1F*    fLogLHRDeDx;
      TH1F*    fLogLHRXs;
      TH1F*    fLogLHRTrk;
					// MC truth
      TH1F*    fPdgCode;		// PDG code of the particle produced most hits
      TH1F*    fFrGH;			// fraction of hits produced by the particle
    };

    struct EventHist_t {
      TH1F*    fRv;			// MC truth information
      TH1F*    fZv;
      TH1F*    fEleCosTh;
      TH1F*    fNClusters;
      TH1F*    fNTracks;
      TH1F*    fNStrawHits[2];
      TH1F*    fNGoodSH;
      TH1F*    fDtClT;
      TH1F*    fEMax;			// energy of the first reco cluster
      TH1F*    fDtClS;
      TH1F*    fSHTime;
      TH1F*    fNHyp;
      TH1F*    fBestHyp[2];		// [0]: by chi2, [1]: by fit consistency

      TH1F*    fNCaloCrystalHits[2];
      TH2F*    fNCaloHitsVsVane[2];
      TH2F*    fNCaloHitsVsRow[2];
      TH2F*    fNCaloHitsVsCol[2];
				  	// calorimeter hit histograms

      TH1F*    fETot        [4];            // total energy/event 
      TH2F*    fECrVsR      [4];            // total energy_per_crystal/event vs radius
      TH2F*    fNCrVsR      [4];            // total energy_per_crystal/event vs radius

      TH2F*    fNCrystalHitsVsR[4];            //
      TH2F*    fNHitCrystalsVsR[4];            //

      TH1F*    fNHitCrystalsTot;
      
    };

    enum { kNEventHistSets   = 100 };
    enum { kNTrackHistSets   = 400 };
    enum { kNClusterHistSets = 100 };
    enum { kNCaloHistSets    = 100 };

    struct Hist_t {
      TH1F*          fCrystalR[2];	          // crystal radius
      EventHist_t*   fEvent  [kNEventHistSets];
      TrackHist_t*   fTrack  [kNTrackHistSets];
      ClusterHist_t* fCluster[kNClusterHistSets];
      CaloHist_t*    fCalo   [kNCaloHistSets];
    };

    Hist_t       fHist;
//-----------------------------------------------------------------------------
// parameters: module names
//-----------------------------------------------------------------------------
    std::string                   fHistFileName;
    std::string                   fG4ModuleLabel;
    std::string                   fStrawHitMaker;
    std::string                   fTrkPatRecDem;
    std::string                   fTrkPatRecUem;
    std::string                   fCrystalHitMaker;
    std::string                   fCaloClusterMaker;
    std::string                   fTrkExtrapol;
    std::string                   fTrkCalMatch;
    std::string                   fPidDem;

    double                        fMinTActive  ;  // start of the active window
    double                        fMinECrystal ;  // 
    double                        fMinCrystalFr;  // min fraction of the crystal area

    KalRepPtrCollection*                fListOfTracks[4]; // hypotheses: de-, ue+, dmu-, umu+
    GenParticleCollection*              fListOfGenParticles;
    StrawHitCollection*                 fListOfStrawHits;
    CaloCrystalHitCollection*           fListOfCaloCrystalHits;
    PtrStepPointMCVectorCollection*     fListOfMcStrawHits;

    TStnTrack*                    fTrack;
    TStnCluster*                  fCluster;
    GenParticle*                  fEle;

    int                           fNClusters;   // number of reconstructed clusters
    int                           fNTracks[4];
    int                           fNStrawHits;  // Nhits in the straw tracker
    int                           fNCaloCrystalHits;
    int                           fNHyp;        // number of hyp's with successfull fits
    int                           fBestHyp[2];  // hypothesis with the best chi2/ndof
    int                           fNGoodTracks; // passed sel C
    int                           fNMatchedTracks; // passed sel C and |dt|<2.5ns

    int                           fCalorimeterType;  // 1:vanes, 2:disks

					// analysis objects

    TStnTrackBlock*               fTrackBlock;
    TCalDataBlock*                fCalDataBlock;
    TStrawDataBlock*              fStrawDataBlock;
    TStnClusterBlock*             fClusterBlock;
    TStnTrackID*                  fTrackID;
    TEmuLogLH*                    fLogLH;

    TDiskCalorimeter*             fDiskCalorimeter;
//-----------------------------------------------------------------------------
// methods
//-----------------------------------------------------------------------------
    explicit TCalm002(fhicl::ParameterSet const& pset);
    ~TCalm002();

    void BookCaloHistograms   (CaloHist_t*    Hist, const char* Folder);
    void BookClusterHistograms(ClusterHist_t* Hist, const char* Folder);
    void BookEventHistograms  (EventHist_t*   Hist, const char* Folder);
    void BookTrackHistograms  (TrackHist_t*   Hist, const char* Folder);
    void BookHistograms();

    void FillCaloHistograms   (CaloHist_t*    Hist, TStnCrystal* Cr );
    void FillClusterHistograms(ClusterHist_t* Hist, TStnCluster* Cls);
    void FillEventHistograms  (EventHist_t*   Hist);
    void FillTrackHistograms  (TrackHist_t*   Hist, TStnTrack* Trk);
    void FillHistograms();

    void DefineBestHypothesis();

    void Debug(const art::Event& Evt);

    void Init(art::Event& Evt);

    //    int  TrackPassedSelectionC(TStnTrack* Track);
//-----------------------------------------------------------------------------
// overloaded methods of TAnaModule
//-----------------------------------------------------------------------------
    void getData(const art::Event& Evt);
//-----------------------------------------------------------------------------
// overloaded framework module methods
//-----------------------------------------------------------------------------
    virtual bool filter  (art::Event& e);
    virtual bool beginRun(art::Run&   r);
    virtual void beginJob();
    virtual void endJob  ();

  };

}

#endif
