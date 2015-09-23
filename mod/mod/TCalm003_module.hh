///////////////////////////////////////////////////////////////////////////////
// Optimize Z-position of the disks
//
// $Id: TCalm003_module.hh,v 1.7 2015/01/06 19:41:31 murat Exp $
// $Author: murat $
// $Date: 2015/01/06 19:41:31 $
//
// Contact person Pavel Murat
///////////////////////////////////////////////////////////////////////////////
#ifndef murat_inc_TCalm003_module_hh
#define murat_inc_TCalm003_module_hh

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
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"


// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

// #include "TrackCaloMatching/inc/TrackClusterLink.hh"
#include "TrackCaloMatching/inc/TrkToCaloExtrapolCollection.hh"

// C++ includes.
#include <iostream>

#include "TString.h"
#include "TFolder.h"
#include "TFile.h"

#include "Stntuple/obj/TStnClusterBlock.hh"
#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/alg/TStnTrackID.hh"

#include "Stntuple/mod/THistModule.hh"

class TAnaTrackID;


namespace mu2e {
  class TTracker;
  class VaneCalorimeter;
}

namespace mu2e {

  class TCalm003 : public THistModule {
  public:

    enum { kNPlanes = 60 };

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
					// MC truth
      TH1F*    fPdgCode;		// PDG code of the particle produced most hits
      TH1F*    fFrGH;			// fraction of hits produced by the particle
    };

    struct EventHist_t {
      TH1F*    fEleCosTh;
      TH1F*    fNTracks;
      TH1F*    fRTrack [kNPlanes];
      TH1F*    fRTrack0[kNPlanes][kNPlanes]; // miss in i1, dist for i2 [i1][i2]
      TH1F*    fRTrack1[kNPlanes][kNPlanes]; // hit  in i1, dist for i2 [i1][i2]

      TH1F*    fRTrack330[kNPlanes][kNPlanes]; // miss in i1@R=330, dist for i2 [i1][i2]
      TH1F*    fRTrack331[kNPlanes][kNPlanes]; // hit  in i1@R=330, dist for i2 [i1][i2]
    };

    enum { kNEventHistSets   = 100 };
    enum { kNTrackHistSets   = 400 };

    struct Hist_t {
      EventHist_t*   fEvent  [kNEventHistSets];
      TrackHist_t*   fTrack  [kNTrackHistSets];
    };

    Hist_t       fHist;
//-----------------------------------------------------------------------------
// parameters
//-----------------------------------------------------------------------------
    std::string                   fHistFileName;
    std::string                   fStrawHitMaker;
    std::string                   fTrkExtrapol;
    std::string                   fTrkCalMatch;
    double                        fMinTActive  ;  // start of the active window

    CaloClusterCollection*              fListOfClusters;
    KalRepPtrCollection*                fListOfTracks[4]; // hypotheses: de-, ue+, dmu-, umu+
    GenParticleCollection*              fListOfGenParticles;
    TrkToCaloExtrapolCollection*        fListOfExtrapolatedTracks;
    StrawHitCollection*                 fListOfStrawHits;
    PtrStepPointMCVectorCollection*     fListOfMcStrawHits;

    TTracker*                     fTracker;

    KalRep*                       fTrack;
    GenParticle*                  fEle;

    int                           fNClusters;   // number of reconstructed clusters
    int                           fNTracks[4];
    int                           fNStrawHits;  // Nhits in the straw tracker
    int                           fNHyp;        // number of hyp's with successfull fits
    int                           fBestHyp[2];  // hypothesis with the best chi2/ndof
    int                           fNGoodTracks; // passed sel C
    int                           fNMatchedTracks; // passed sel C and |dt|<2.5ns

					// analysis objects

    TStnTrackBlock*               fTrackBlock;
    TStnClusterBlock*             fClusterBlock;
    TStnTrackID*                  fTrackID;

    double                        fZPlane[kNPlanes];
    double                        fZDisk1;

					// assume less than 100 tracks
    HepPoint                      fPos[100][kNPlanes];
    HepPoint                      fPosDisk1[100];
    double                        fRDisk1  [100];
//-----------------------------------------------------------------------------
// methods
//-----------------------------------------------------------------------------
    explicit TCalm003(fhicl::ParameterSet const& pset);
    ~TCalm003();

    void BookEventHistograms  (EventHist_t*   Hist, const char* Folder);
    void BookTrackHistograms  (TrackHist_t*   Hist, const char* Folder);
    void BookHistograms();

    void FillEventHistograms  (EventHist_t*   Hist);
    void FillTrackHistograms  (TrackHist_t*   Hist, TStnTrack* Trk);
    void FillHistograms();

    void Debug(art::Event* Evt);
    void Init (art::Event* Evt);
//-----------------------------------------------------------------------------
// overloaded methods of TAnaModule
//-----------------------------------------------------------------------------
    void getData(art::Event* Evt);
//-----------------------------------------------------------------------------
// overloaded framework module methods
//-----------------------------------------------------------------------------
    virtual bool filter  (art::Event& e);
    virtual void beginJob();
    virtual void endJob  ();

  };

}

#endif
