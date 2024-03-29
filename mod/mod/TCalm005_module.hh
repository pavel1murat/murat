///////////////////////////////////////////////////////////////////////////////
// study origin of the particles hitting the Mu2e calorimeter
///////////////////////////////////////////////////////////////////////////////
//
// $Id: TCalm005_module.hh,v 1.5 2015/01/06 19:41:31 murat Exp $
// $Author: murat $
// $Date: 2015/01/06 19:41:31 $
//
// Contact person Pavel Murat
//
#ifndef murat_mod_mod_TCalm005_module_hh
#define murat_mod_mod_TCalm005_module_hh

// Mu2e includes.
// #include "CLHEP/Geometry/HepPoint.h"
#include "BTrk/BbrGeom/HepPoint.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Vector.h"

#include "BTrk/TrkBase/TrkHelixUtils.hh"
#include "BTrk/TrkBase/HelixParams.hh"
#include "BTrk/KalmanTrack/KalHit.hh"

// storable objects (data products)

#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"

#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"

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

  class TCalm005 : public THistModule {

  public:

    struct StepPointMCHist_t {
      TH2F*  fRVsZ;
      TH2F*  fR0VsZ0;
      TH1F*  fPdgCode[2];		// same histogram, different ranges
      TH1F*  fTime;
      TH1F*  fLength;
      TH1F*  fEDep;
    };


    struct EventHist_t {
      TH1F*    fEventNumber;
      TH1F*    fRunNumber;
      TH1F*    fNSteps;
    };

    enum { kNEventHistSets   = 100 };
    enum { kNStepPointMCHistSets = 100 };

    struct Hist_t {
      EventHist_t*       fEvent      [kNEventHistSets];
      StepPointMCHist_t* fStepPointMC[kNStepPointMCHistSets];
    };

    Hist_t       fHist;
//-----------------------------------------------------------------------------
// parameters
//-----------------------------------------------------------------------------
    std::string                     fG4RunModuleLabel; // default: "g4run"
    std::string                     fHistFileName;
    std::string                     fProcessName;    // default : ""
    std::string                     fProductName;   // "calorimeter"
    double                          fMinTActive  ;  // start of the active window

    const StepPointMCCollection*    fStepPointMCColl;

    int                             fNSteps; // number of step points
//-----------------------------------------------------------------------------
// methods
//-----------------------------------------------------------------------------
    explicit TCalm005(fhicl::ParameterSet const& pset);
    ~TCalm005();

    void BookStepPointMCHistograms(StepPointMCHist_t* Hist, const char* Folder);
    void BookEventHistograms      (EventHist_t*       Hist, const char* Folder);
    void BookHistograms();

    void FillStepPointMCHistograms(StepPointMCHist_t* Hist, StepPointMC* Step);
    void FillEventHistograms  (EventHist_t*   Hist);
    void FillHistograms();

    void Debug(const art::Event* Evt);
    void Init (const art::Event* Evt);
//-----------------------------------------------------------------------------
// overloaded methods of TAnaModule
//-----------------------------------------------------------------------------
    void getData(const art::Event* Evt);
//-----------------------------------------------------------------------------
// overloaded framework module methods
//-----------------------------------------------------------------------------
    virtual void analyze (const art::Event& e);
    virtual void beginRun(const art::Run&   r);
    virtual void beginJob();
    virtual void endJob  ();

  };

}

#endif
