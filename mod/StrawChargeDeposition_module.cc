///////////////////////////////////////////////////////////////////////////////
//
// .fcl file to use: murat/test/strawchargeDeposition.fcl
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Selector.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/VirtualDetector.hh"

#include "TrackerGeom/inc/Tracker.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "Mu2eUtilities/inc/SimParticlesWithHits.hh"
#include "Mu2eUtilities/inc/SortedStepPoints.hh"
#include "Mu2eUtilities/inc/TrackTool.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"

#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"

#include "BTrk/TrkBase/HelixParams.hh"
#include "BTrk/KalmanTrack/KalHit.hh"

#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"

#include "Stntuple/mod/StntupleModule.hh"

#include "TH1F.h"

using namespace std;
using CLHEP::Hep3Vector;

namespace mu2e {

  class StrawChargeDeposition : public art::EDFilter {
  private:
//-----------------------------------------------------------------------------
// Module labels 
//-----------------------------------------------------------------------------
    std::string       _spmcCollTag;
    const Tracker*    _tracker;     // straw tracker

					// hit flag bits which should be ON and OFF

    const mu2e::StepPointMCCollection*  _spmcc;            //

    struct Hist_t {
      TH1F*   _energy;                      // step point energy deposition
      TH1F*   _energyVsPlane;               // 
      TH2F*   _strawVsRho[40];
    } ;

    Hist_t _hist;

  public:
    explicit StrawChargeDeposition(fhicl::ParameterSet const& pset);
    virtual ~StrawChargeDeposition();

    void     bookHistograms();
    void     getData(const art::Event* Evt);
    void     Init   (art::Event* Evt);
//-----------------------------------------------------------------------------
// overloaded virtual methods of the base class
//-----------------------------------------------------------------------------
    virtual void     beginJob();
    virtual bool     beginRun(art::Run& );
    virtual void     endJob  ();
    virtual bool     filter (art::Event& Evt);
  };


//-----------------------------------------------------------------------------
  StrawChargeDeposition::StrawChargeDeposition(fhicl::ParameterSet const& pset): 
    EDFilter            (pset),
    _spmcCollTag              (pset.get<std::string>("spmcCollTag"            ))
  {

  }

//-----------------------------------------------------------------------------
  StrawChargeDeposition::~StrawChargeDeposition() { 
    //    delete fHist.fDr;
  }

//-----------------------------------------------------------------------------
  void StrawChargeDeposition::endJob() {
    //    art::ServiceHandle<art::TFileService> tfs;
  }


//-----------------------------------------------------------------------------
  void StrawChargeDeposition::bookHistograms() {

    art::ServiceHandle<art::TFileService> tfs;

    _hist._energy        = tfs->make<TH1F>("energy"    ,"Step point energy deposition", 1500,0,  0.15);
    _hist._energyVsPlane = tfs->make<TH1F>("e_vs_plane","E vs plane", 50,0,  50);

    for (int i=0; i<36; i++) {
      _hist._strawVsRho[i]    = tfs->make<TH2F>(Form("straw_vs_rho_%02i",i), Form("straw vs rhoplane %02i",i),
						100,-1000, 1000,100,0,100);
    }
  }
//-----------------------------------------------------------------------------
  void StrawChargeDeposition::beginJob() {
    bookHistograms();
  }

//-----------------------------------------------------------------------------
  bool StrawChargeDeposition::beginRun(art::Run& ) {
    mu2e::GeomHandle<mu2e::Tracker> th;
    _tracker = th.get();
    return true;
  }


//-----------------------------------------------------------------------------
// get data from the event record
//-----------------------------------------------------------------------------
  void StrawChargeDeposition::getData(const art::Event* Evt) {

    art::Handle<mu2e::StepPointMCCollection> spmcch;

    _spmcc = nullptr;

    Evt->getByLabel(_spmcCollTag,spmcch);
    
    if (spmcch.isValid()) _spmcc = spmcch.product();
  }

//-----------------------------------------------------------------------------
  void StrawChargeDeposition::Init(art::Event* Evt) {
  }


  //-----------------------------------------------------------------------------
  bool StrawChargeDeposition::filter(art::Event& Evt) {
    // const char* oname = "StrawChargeDeposition::filter";

    //    printf("[%s] RUN: %10i EVENT: %10i\n",oname,Evt.run(),Evt.event());

    int plane, panel, layer, ist;

    getData(&Evt);

    int ns = _spmcc->size();

    for (int i=0; i<ns; i++) {
      const StepPointMC* step = &_spmcc->at(i);
      
      const Straw& straw = _tracker->getStraw(step->strawId());

      plane = straw.id().getPlane();
      panel = straw.id().getPanel();
      layer = straw.id().getLayer();
      ist   = straw.id().getStraw();

      if (i < 0) printf("%5i %5i %5i %5i \n",plane, panel,layer,ist);

      const CLHEP::Hep3Vector& smp = straw.getMidPoint();

      double dx = step->position().x()-smp.x();
      double dy = step->position().y()-smp.y();

      double rho = dx*smp.unit().y()-dy*smp.unit().x();

      //      double rho = sqrt(dx*dx+dy*dy);

      _hist._energy->Fill(step->ionizingEdep());
      _hist._energyVsPlane->Fill(plane,step->ionizingEdep());
      _hist._strawVsRho[plane]->Fill(rho,ist,step->ionizingEdep());
    }

    return true;
  }
    
}

using mu2e::StrawChargeDeposition;
DEFINE_ART_MODULE(StrawChargeDeposition);
