///////////////////////////////////////////////////////////////////////////////
// A half-interactive 2D event display. 
//
// $Id: GaasSandwichAna_module.cc,v 1.11 2014/10/02 17:15:09 murat Exp $
// $Author: murat $
// $Date: 2014/10/02 17:15:09 $
//
// Contact person:  Pavel Murat
//
// Debug_003: look at the systematics between the StepPointMCs and StrawHitPosition's
// Debug_004: look at various hit-level MC distributions
//
// .fcl file to use: murat/test/trackerMCCheck.fcl
///////////////////////////////////////////////////////////////////////////////

// C++ includes.
#include <iostream>
#include <string>

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Principal/Selector.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes.
#include "GeometryService/inc/VirtualDetector.hh"

#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"

#include "Stntuple/mod/StntupleModule.hh"

// ROOT includes
// #include "TApplication.h"
// #include "TArc.h"
// #include "TArrow.h"
// #include "TCanvas.h"
// #include "TDirectory.h"
// #include "TGraph.h"
#include "TH1F.h"
// #include "TLine.h"
// #include "TBox.h"
// #include "TMarker.h"
// #include "TEllipse.h"
// #include "TText.h"
// #include "TNtuple.h"

// Other includes
// #include "CLHEP/Units/SystemOfUnits.h"

using namespace std;
using CLHEP::Hep3Vector;

namespace mu2e {

  class GaasSandwichAna : public StntupleModule {
  private:
//-----------------------------------------------------------------------------
// Module labels 
//-----------------------------------------------------------------------------
    std::string        fModuleLabel;	             // this module label
    std::string        processName_;
    std::string        fG4ModuleLabel;
    
    std::string        producerName_;
    std::string        fStrawHitMaker;
    std::string        fStrawHitPosMaker;
    std::string        fFlagBgrHitsModuleLabel;
    
    int                fPdgCode;
    int                fGeneratorCode;

    const mu2e::StepPointMCCollection*           fSteps;            //

    struct Hist_t {
      TH1F*   fEDepActive;                  // 
      TH1F*   fEDepPassive;                  // 
      TH1F*   fEDepTotal;                  // 
      TH1F*   fNActive;                    // 
      TH1F*   fFirst;                    // 
    } ;

    Hist_t fHist;


  public:
    explicit GaasSandwichAna(fhicl::ParameterSet const& pset);
    virtual ~GaasSandwichAna();

    void     getData(const art::Event* Evt);
    void     Init   (art::Event* Evt);
    void     Debug_003();   // handles fDr
//-----------------------------------------------------------------------------
// overloaded virtual methods of the base class
//-----------------------------------------------------------------------------
    virtual void     beginJob();
    virtual void     endJob  ();
    virtual bool     filter (art::Event& Evt);
  };


//-----------------------------------------------------------------------------
  GaasSandwichAna::GaasSandwichAna(fhicl::ParameterSet const& pset): 
    StntupleModule            (pset,"GaasSandwichAna"),
    fModuleLabel              (pset.get<std::string>("module_label"                 )),
    processName_              (pset.get<std::string>("processName"          ,""     )),
    fG4ModuleLabel            (pset.get<std::string>("g4ModuleLabel"        ,"g4run"))
  {

  }

//-----------------------------------------------------------------------------
  GaasSandwichAna::~GaasSandwichAna() { 
  }


//-----------------------------------------------------------------------------
  void GaasSandwichAna::endJob() {
    art::ServiceHandle<art::TFileService> tfs;
  }
//-----------------------------------------------------------------------------
  void GaasSandwichAna::beginJob() {

    art::ServiceHandle<art::TFileService> tfs;

//     art::TFileDirectory tfdir = tfs->mkdir( "CosmicDYB" );
//     _cosmicMultiplicityH = tfdir.make<TH1D>( "MultiplicityH", "Cosmic Multiplicity", 20, -0.5, 19.5);

    fHist.fEDepActive  = tfs->make<TH1F>("edepa" ,"Energy Deposition in GaAs ", 200,0, 1);
    fHist.fEDepPassive = tfs->make<TH1F>("edepp" ,"Energy Deposition in Lead ", 200,0, 1);
    fHist.fEDepTotal   = tfs->make<TH1F>("edept" ,"Energy Deposition Total   ", 200,0, 1);
    fHist.fFirst       = tfs->make<TH1F>("first" ,"First active layer"        , 500,0, 500);
    fHist.fNActive     = tfs->make<TH1F>("nactv" ,"N(active layers w/signals)", 500,0, 500);
//-----------------------------------------------------------------------------
// define collection names to be used for initialization
//-----------------------------------------------------------------------------
  }

//-----------------------------------------------------------------------------
// get data from the event record
//-----------------------------------------------------------------------------
  void GaasSandwichAna::getData(const art::Event* Evt) {
    //    int   rc (0);
    //    const char* oname = "GaasSandwichAna::getData";

    art::Handle<StepPointMCCollection> stepsHandle;
    art::Selector sel(art::ProductInstanceNameSelector("stepper") &&
		      art::ProcessNameSelector(processName_) &&
		      art::ModuleLabelSelector(fG4ModuleLabel)  );
    Evt->get(sel, stepsHandle);
    if (stepsHandle.isValid()) fSteps =  (const mu2e::StepPointMCCollection*) &(*stepsHandle);
    else                       fSteps = NULL;
  }

//-----------------------------------------------------------------------------
  void GaasSandwichAna::Init(art::Event* Evt) {
  }


  //-----------------------------------------------------------------------------
  bool GaasSandwichAna::filter(art::Event& Evt) {
    const char* oname = "GaasSandwichAna::filter";

    printf("[%s] RUN: %10i EVENT: %10i\n",oname,Evt.run(),Evt.event());

//-----------------------------------------------------------------------------
// get event data and initialize data blocks
//-----------------------------------------------------------------------------
    getData(&Evt);

    Debug_003();
    //    if  (DebugBit(3)) Debug_003();

    return true;
  }

//-----------------------------------------------------------------------------
// plot distribution in radial distance between the StepPointMC and the 
// corresponding hit - fHist.fDr
// so far assume 400 layers...
//-----------------------------------------------------------------------------
  void GaasSandwichAna::Debug_003() {
    
    int nsteps;

    const mu2e::StepPointMC        *step;

    nsteps = fSteps->size();

    int hit_layer[400];

    for (int i=0; i<400; i++) hit_layer[i] = 0;
    
    double edep_active = 0;
    double edep_passive = 0;

    for (int i=0; i<nsteps; i++) {
      step =  &fSteps->at(i);
      int vid  = step->volumeId();
      double edep = step->totalEDep();
      if ((vid >= 100) && (vid < 500)) {
	edep_active += edep;
	if (edep > 0) hit_layer[vid-100] = 1;
      }
      else if ( vid >= 500) {
	edep_passive += edep;
      }
    }

    int nactv = 0;
    int first(-1);
    for (int i=0; i<400; i++) {
      if (hit_layer[i] != 0) {
	if (nactv == 0) first = i;
	nactv += 1;
      }
    }

    fHist.fEDepActive->Fill(edep_active);
    fHist.fEDepPassive->Fill(edep_passive);
    fHist.fEDepTotal->Fill(edep_active+edep_passive);
    fHist.fNActive->Fill(nactv);
    fHist.fFirst->Fill(first);
  }
    
}

using mu2e::GaasSandwichAna;
DEFINE_ART_MODULE(GaasSandwichAna);
