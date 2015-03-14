///////////////////////////////////////////////////////////////////////////////
// study origin of the particles hitting the Mu2e calorimeter
///////////////////////////////////////////////////////////////////////////////
//
// $Id: TCalm006_module.hh,v 1.3 2015/01/06 19:41:31 murat Exp $
// $Author: murat $
// $Date: 2015/01/06 19:41:31 $
//
// Contact person Pavel Murat
//
#ifndef murat_mod_mod_TCalm006_module_hh
#define murat_mod_mod_TCalm006_module_hh

// Mu2e includes.
#include "CLHEP/Geometry/HepPoint.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Vector.h"

// storable objects (data products)
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"

//-----------------------------------------------------------------------------
// keep deleted from CalorimeterGeometry class locally
//-----------------------------------------------------------------------------
#include "murat/obj/HexMap.hh"

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

// C++ includes.
#include <iostream>

// #include "TString.h"
#include "TFolder.h"
#include "TFile.h"

#include "Stntuple/obj/TStnClusterBlock.hh"
#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/mod/THistModule.hh"

namespace mu2e {

  class TCalm006 : public THistModule {

  public:

    struct EventHist_t {
      TH1F*    fEventNumber;
      TH1F*    fRunNumber;
      TH1F*    fEEle;
      TH1F*    fEnergy[100];
      TH1F*    fESum  [100];
    };

    enum { kNEventHistSets   = 100 };

    struct Hist_t {
      EventHist_t*       fEvent      [kNEventHistSets];
    };

    Hist_t       fHist;
//-----------------------------------------------------------------------------
// parameters
//-----------------------------------------------------------------------------
    std::string                     fG4RunModuleLabel; // default: "g4run"
    std::string                     fHistFileName;
    std::string                     fProcessName;      // default : ""
    std::string                     fProductName;      // "calorimeter"
    double                          fMinTActive  ;     // start of the active window
    std::string                     fCrystalHitMaker;

    int                             fNCaloCrystalHits;
    CaloCrystalHitCollection*       fListOfCaloCrystalHits;
    HexMap*                         fHexMap;

    double                          fEnergy[10];
    double                          fESum  [10];
//-----------------------------------------------------------------------------
// methods
//-----------------------------------------------------------------------------
    explicit TCalm006(fhicl::ParameterSet const& pset);
    ~TCalm006();

    void BookEventHistograms      (EventHist_t*       Hist, const char* Folder);
    void BookHistograms();

    void FillEventHistograms  (EventHist_t*   Hist);
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
    virtual bool beginRun(art::Run&   r);
    virtual void beginJob();
    virtual void endJob  ();

  };

}

#endif
