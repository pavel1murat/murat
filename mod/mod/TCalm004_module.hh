///////////////////////////////////////////////////////////////////////////////
// rejection of teh cosmics-associated background
///////////////////////////////////////////////////////////////////////////////
//
// $Id: TCalm004_module.hh,v 1.14 2015/01/06 19:41:31 murat Exp $
// $Author: murat $
// $Date: 2015/01/06 19:41:31 $
//
// Contact person Pavel Murat
//
#ifndef murat_mod_mod_TCalm004_module_hh
#define murat_mod_mod_TCalm004_module_hh

// Mu2e includes.
// #include "CLHEP/Geometry/HepPoint.h"
#include "BTrk/BbrGeom/HepPoint.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Vector.h"

#include "BTrk/TrkBase/TrkHelixUtils.hh"
#include "BTrk/TrkBase/HelixParams.hh"
#include "BTrk/KalmanTrack/KalHit.hh"

#include "KalmanTests/inc/KalRepPtrCollection.hh"
// storable objects (data products)
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"

#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"


// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "TrackCaloMatching/inc/TrkToCaloExtrapolCollection.hh"

// C++ includes.
#include <iostream>

// #include "TString.h"
#include "TFolder.h"
#include "TFile.h"

#include "Stntuple/obj/TStnClusterBlock.hh"
#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/mod/THistModule.hh"

#include "Stntuple/alg/TEmuLogLH.hh"

class TStnTrackID;

namespace mu2e {

  class TCalm004 : public THistModule {

  public:

    enum {
      kClusterBit     = 0x1 << 0,
      kClusterDxBit   = 0x1 << 1,
      kTrackD0Bit     = 0x1 << 2,
      kLHRCalBit      = 0x1 << 3,
      kLHRTrkBit      = 0x1 << 4
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
      TH1F*    fP;			// total momentum
      TH1F*    fPt;
      TH1F*    fCosTh;
      TH1F*    fChi2;
      TH1F*    fNDof;
      TH1F*    fChi2Dof;
      TH1F*    fNActive;
      TH1F*    fT0;
      TH1F*    fQ;
      TH1F*    fFitCons;		// fit consistency (0 to 1)
      TH1F*    fD0;
      TH1F*    fZ0;
      TH1F*    fTanDip;
					// matching histograms
      TH1F*    fNClusters;
      TH1F*    fVaneID;
      TH1F*    fXCal;
      TH1F*    fYCal;
      TH1F*    fZCal;
      TH1F*    fYTrk;
      TH1F*    fZTrk;
      TH1F*    fDt;			// track-cluster residuals
      TH1F*    fDy;
      TH1F*    fDz;
      TH1F*    fEp;
      TH2F*    fNHVsStation;
      TH2F*    fNHVsNSt;
      
      TH1F*    fRSlope;
      TH1F*    fXSlope;

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
      TH1F*    fEleCosTh;
      TH1F*    fRv;
      TH1F*    fZv;
      TH1F*    fNClusters;
      TH1F*    fNTracksDem;
      TH1F*    fNTracksUem;
      TH1F*    fNStrawHits[2];
      TH1F*    fNGoodSH;
      TH1F*    fDtClT;
      TH1F*    fEMax;			// energy of the first reco cluster
      TH1F*    fDtClS;
      TH1F*    fNHyp;
      TH1F*    fBestHyp[2];		// [0]: by chi2, [1]: by fit consistency
    };

    enum { kNEventHistSets   = 100 };
    enum { kNTrackHistSets   = 400 };
    enum { kNClusterHistSets = 100 };

    struct Hist_t {
      EventHist_t*   fEvent  [kNEventHistSets];
      TrackHist_t*   fTrack  [kNTrackHistSets];
      ClusterHist_t* fCluster[kNClusterHistSets];
    };

    Hist_t       fHist;
//-----------------------------------------------------------------------------
// parameters
//-----------------------------------------------------------------------------
    std::string                     fHistFileName;
    std::string                     fStrawHitMaker;
    std::string                     fClusterMaker;

    std::string                     fTrkPatRecDem;
    std::string                     fTrkPatRecUem; // the 2nd leg of the electron, if any
    std::string                     fTrkExtrapolDem;
    std::string                     fTrkCalMatchDem;
    std::string                     fTrkExtrapolUem;
    std::string                     fTrkCalMatchUem;

    std::string                     fPidDem; // Vadim's particle ID module label

    double                          fMinTActive  ;  // start of the active window
    double                          fMinECrystal ;  // 

    CaloClusterCollection*          fListOfClusters;
    KalRepPtrCollection*            fListOfTracks[4]; // hypotheses: de-, ue+, dmu-, umu+
    GenParticleCollection*          fListOfGenParticles;
    TrkToCaloExtrapolCollection*    fListOfExtrapolatedTracks;
    StrawHitCollection*             fListOfStrawHits;
    PtrStepPointMCVectorCollection* fListOfMcStrawHits;

    TEmuLogLH*                      fLogLH;

    TStnTrack*                      fTrack;
    TStnCluster*                    fCluster;

    GenParticle*                    fEle;

    int                             fNClusters;   // number of reconstructed clusters
    int                             fNTracks[2];  // 
    int                             fNStrawHits;  // Nhits in the straw tracker
    int                             fNHyp;        // number of hyp's with successfull fits
    int                             fBestHyp[2];  // hypothesis with the best chi2/ndof
    int                             fNGoodTracks; // passed sel C
    int                             fNMatchedTracks; // passed sel C and |dt|<2.5ns

    int                             fCosmicFlag;
    int                             fPassed;

					           // analysis objects
//-----------------------------------------------------------------------------
// given there are multiple lists of objects , out of which one has to assemble the 
// analysis objects, it seems easier to maintain 2 track blocks, one - for the downstream 
// negative tracks and the second one - for the upstream tracks.
// the cluster block, however, is just one...
//-----------------------------------------------------------------------------
    TStnTrackBlock*                 fTrackBlock[2]; // [0]: DEM, [1]: UEM
    TStnClusterBlock*               fClusterBlock;
    TStnTrackID*                    fTrackIDSetB;
    TStnTrackID*                    fTrackIDSetC;
    int                             fTrackIDWordSetB[100]; // assume N(tracks) < 100
    int                             fElectronIDWord [100];
//-----------------------------------------------------------------------------
// methods
//-----------------------------------------------------------------------------
    explicit TCalm004(fhicl::ParameterSet const& pset);
    ~TCalm004();

    void BookClusterHistograms(ClusterHist_t* Hist, const char* Folder);
    void BookEventHistograms  (EventHist_t*   Hist, const char* Folder);
    void BookTrackHistograms  (TrackHist_t*   Hist, const char* Folder);
    void BookHistograms();

    void FillClusterHistograms(ClusterHist_t* Hist, TStnCluster* Cls);
    void FillEventHistograms  (EventHist_t*   Hist);
    void FillTrackHistograms  (TrackHist_t*   Hist, TStnTrack* Trk);
    void FillHistograms();

    void Debug(const art::Event& Evt);

    void Init(art::Event& Evt);

    //    int  TrackPassedSelectionC(TAnaTrack* Track);
//-----------------------------------------------------------------------------
// overloaded methods of TAnaModule
//-----------------------------------------------------------------------------
    void getData(const art::Event& Evt);
//-----------------------------------------------------------------------------
// overloaded framework module methods
//-----------------------------------------------------------------------------
    virtual bool filter  (art::Event& e);
    virtual void beginJob();
    virtual void endJob  ();

  };

}

#endif
