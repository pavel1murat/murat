///////////////////////////////////////////////////////////////////////////////
// rejection of teh cosmics-associated background
///////////////////////////////////////////////////////////////////////////////
//
// $Id: PetAna001_module.hh,v 1.9 2014/01/19 18:50:48 murat Exp $
// $Author: murat $
// $Date: 2014/01/19 18:50:48 $
//
// Contact person Pavel Murat
//
#ifndef murat_mod_mod_PetAna001_module_hh
#define murat_mod_mod_PetAna001_module_hh

// Mu2e includes.
#include "CLHEP/Geometry/HepPoint.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Vector.h"

// storable objects (data products)
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"

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

// #include "TString.h"
#include "TFolder.h"
#include "TFile.h"

#include "Stntuple/mod/THistModule.hh"
#include "murat/pet/BrainImager.hh"

namespace mu2e {

  class PetAna001 : public THistModule {

  public:

    enum {
      kClusterBit     = 0x1 << 0,
      kClusterDxBit   = 0x1 << 1,
      kTrackD0Bit     = 0x1 << 2,
      kLHRBit         = 0x1 << 3
    };

    struct PhotonHist_t {
      TH1F*    fTime;
      TH1F*    fCosTh;
      TH1F*    fPhi;
      TH1F*    fZ;			// Z
      TH1F*    fR;			// radius
      TH2F*    fRVsZ;			// production point z:radius
    };


    struct EventHist_t {
      TH1F*    fNPhotons;
      TH1F*    fNHitCrystals;
      TH1F*    fNHitModules;
      TH1F*    fEMax;
      TH1F*    fETot;
      TH2F*    fE2VsE1Cr;
      TH2F*    fE2VsE1Wed;
      TH2F*    fW2VsW1;
      TH1F*    fHitTime;
      TH1F*    fDtMin;
      TH1F*    fDist;
      TH1F*    fCType;
    };

    enum { 
      kNEventHistSets   = 100,
      kNPhotonHistSets  = 100
    };

    struct Hist_t {
      PhotonHist_t*  fPhoton [kNPhotonHistSets];
      EventHist_t*   fEvent  [kNEventHistSets ];
    };

    Hist_t       fHist;
//-----------------------------------------------------------------------------
// parameters
//-----------------------------------------------------------------------------
    std::string                     fHistFileName;
    std::string                     fGenParticleMaker;
    std::string                     fCaloCrystalHitMaker;
    double                          fMinPhotopeakE;
    double                          fMaxPhotopeakE;
    double                          fTriggerDt;
    int                             fPrintFrequency;
//-----------------------------------------------------------------------------
// other variables
//-----------------------------------------------------------------------------
    GenParticleCollection*          fListOfGenParticles;
    GenParticle*                    fPho;

    const CaloCrystalHitCollection* fListOfCaloCrystalHits;
    const StepPointMCCollection*    fStepPointMCCollection;

    

    BrainImager*                    fBrainImager;
    int                             fNWedges;

    int                             fNHitCrystals;   // number of crystals with hits
    int                             fNHitModules;
    int                             fNHitsPerModule[50];
    double                          fETotModule    [50];
    double                          fEMaxModule    [50];
    double                          fTMaxModule    [50];
    int                             fModIndex      [50];
    CLHEP::Hep3Vector               fHitPos        [50];
    const mu2e::CaloCrystalHit*     fHit           [50];

    int                             fNPhotons;

    int                             fTrigIndex[2]; // 2 "trigger" modules
    double                          fDtMin;
    double                          fDistance;     // vertex to straight line
    int                             fCoincidenceType;
//-----------------------------------------------------------------------------
// methods
//-----------------------------------------------------------------------------
    explicit PetAna001(fhicl::ParameterSet const& pset);
    ~PetAna001();

    void BookEventHistograms  (EventHist_t*   Hist, const char* Folder);
    void BookPhotonHistograms (PhotonHist_t*  Hist, const char* Folder);
    void BookHistograms();

    void FillPhotonHistograms (PhotonHist_t*  Hist, mu2e::GenParticle* Photon);
    void FillEventHistograms  (EventHist_t*   Hist);
    void FillHistograms();

    void Debug(const art::Event& Evt);

    void Init(art::Event& Evt);
					// find annihilation photon hitting 
					// given crystal

    const mu2e::SimParticle* FindPhoton(int CrystalID);
//-----------------------------------------------------------------------------
// overloaded methods of TAnaModule - to come...
//-----------------------------------------------------------------------------
    int getData(const art::Event& Evt);
//-----------------------------------------------------------------------------
// overloaded framework module methods
//-----------------------------------------------------------------------------
    virtual bool filter  (art::Event& Evt);
    virtual void beginJob();
    virtual bool beginRun(art::Run& Run);
    virtual void endJob  ();

  };

}

#endif
