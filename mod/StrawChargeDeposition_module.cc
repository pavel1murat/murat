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
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/VirtualDetector.hh"

#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/Mu2eUtilities/inc/SortedStepPoints.hh"
#include "Offline/Mu2eUtilities/inc/TrackTool.hh"
#include "Offline/Mu2eUtilities/inc/TwoLinePCA.hh"

#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/DataProducts/inc/VirtualDetectorId.hh"

#include "BTrk/TrkBase/HelixParams.hh"
#include "BTrk/KalmanTrack/KalHit.hh"

#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"

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
    int               _debugBit;
    const Tracker*    _tracker;         // straw tracker

					// hit flag bits which should be ON and OFF

    const mu2e::StepPointMCCollection*  _spmcc;            //

    struct Hist_t {
      TH1F*   _stepEnergy;                   // step point energy deposition
      TH1F*   _totEnergyVsPlane;               // 
      TH1F*   _stepEnergyPlane[40];        // 
      TH2F*   _strawVsRho[40];
      TH2F*   _y_vs_x[40];
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
    _spmcCollTag        (pset.get<std::string>("spmcCollTag")),
    _debugBit           (pset.get<int>        ("debugBit"   ))
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

    _hist._stepEnergy        = tfs->make<TH1F>("eStep"    ,"Step point energy deposition", 1500,0,0.15);
    _hist._totEnergyVsPlane = tfs->make<TH1F>("e_vs_plane","total E vs plane", 50,0,  50);

    for (int i=0; i<36; i++) {
      _hist._stepEnergyPlane[i] = tfs->make<TH1F>(Form("eStep_%02i",i),Form("E step in plane %02i",i),1500,0,0.15);

      _hist._strawVsRho[i]      = tfs->make<TH2F>(Form("straw_vs_rho_%02i",i), Form("straw vs rho plane %02i",i),
						100,-1000, 1000,100,0,100);

      _hist._y_vs_x[i]          = tfs->make<TH2F>(Form("y_vs_x_%02i",i), Form("Y vs X plane %02i",i),
						  100,-1000, 1000,100,-1000,1000);
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
    const char* oname = "StrawChargeDeposition::filter";

    //    printf("[%s] RUN: %10i EVENT: %10i\n",oname,Evt.run(),Evt.event());

    int plane, panel, layer, ist;

    getData(&Evt);

    int ns = _spmcc->size();

    for (int i=0; i<ns; i++) {
      const StepPointMC* step = &_spmcc->at(i);
      const SimParticle* simp = step->simParticle().get();
      
      const Straw& straw = _tracker->getStraw(step->strawId());

      plane = straw.id().getPlane();
      panel = straw.id().getPanel();
      layer = straw.id().getLayer();
      ist   = straw.id().getStraw();

      if (i < 0) printf("%5i %5i %5i %5i \n",plane, panel,layer,ist);

      const CLHEP::Hep3Vector& smp = straw.getMidPoint();
      const CLHEP::Hep3Vector& dir = straw.getDirection();

      double x  = step->position().x();
      double y  = step->position().y();
      double dx = x-smp.x();
      double dy = y-smp.y();

      double rho  = dx*smp.unit().y()-dy*smp.unit().x();
      double rho1 = dx*dir.x()+dy*dir.y();


      if ((_debugBit == 1) && (fabs(rho) > 450.)) {
	printf("[%s] RUN: %10i subrun: %10i EVENT: %10i PDG mom iplane panel ist rho: %10i %12.4f %5i %5i %5i %12.4f %12.4f\n",
	       oname,Evt.run(),Evt.subRun(),Evt.event(),
	       simp->pdgId(),simp->startMomentum().vect().mag(),plane,panel,ist,rho,rho1
	       );
      }

      _hist._stepEnergy->Fill(step->ionizingEdep());

      _hist._totEnergyVsPlane->Fill(plane,step->ionizingEdep());

      _hist._stepEnergyPlane[plane]->Fill(step->ionizingEdep());

      _hist._strawVsRho[plane]->Fill(rho1,ist,step->ionizingEdep());
      _hist._y_vs_x    [plane]->Fill(x,y,step->ionizingEdep());
    }

    return true;
  }
    
}

using mu2e::StrawChargeDeposition;
DEFINE_ART_MODULE(StrawChargeDeposition);
