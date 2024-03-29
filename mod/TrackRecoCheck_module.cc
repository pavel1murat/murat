///////////////////////////////////////////////////////////////////////////////
// $Id: TrackRecoCheck_module.cc,v 1.3 2014/11/17 21:29:51 murat Exp $
// $Author: murat $
// $Date: 2014/11/17 21:29:51 $
//
//
// 2014-08-25: check various aspects of the track and hit reconstruction
//
// bit  1: print dw, dx, dy
// bit 11: plot distribution in hit-to-track distance along the wire
//
// .fcl file to use: murat/test/trackRecoCheck.fcl
///////////////////////////////////////////////////////////////////////////////

// C++ includes.
#include <iostream>
#include <string>

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Selector.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes.
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/ConditionsService/inc/ConditionsHandle.hh"
// #include "ConditionsService/inc/TrackerCalibrations.hh"
#include "Offline/TrackerConditions/inc/Mu2eDetector.hh"

// #include "BTrkHelper/inc/BTrkHelper.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"

#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/Mu2eUtilities/inc/SortedStepPoints.hh"
#include "Offline/Mu2eUtilities/inc/TrackTool.hh"

#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"

#include "BTrk/TrkBase/HelixParams.hh"
#include "BTrk/KalmanTrack/KalHit.hh"

#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"

#include "Offline/BTrkData/inc/TrkStrawHit.hh"
#include "Offline/RecoDataProducts/inc/KalRepPtrCollection.hh"

#include "Stntuple/mod/StntupleModule.hh"

// #include "Mu2eBTrk/inc/Mu2eDetectorModel.hh"

#include "TH1F.h"

using namespace std;
using CLHEP::Hep3Vector;

namespace mu2e {

  class TrackRecoCheck : public StntupleModule {
  private:
//-----------------------------------------------------------------------------
// Module labels 
//-----------------------------------------------------------------------------
    std::string        fProcessName;
    std::string        fG4ModuleLabel;
    std::string        fTrkPatRecModuleLabel;
    
    struct Hist_t {
      TH1F*   fHitDw;                   // radial distance between the StepPointMC and the corresponding hit
      TH1F*   fHitRes;			// hit residual
      TH1F*   fPath;
      TH2F*   fEHitVsPath;
    } ;

    Hist_t fHist;

    ProditionsHandle<Mu2eDetector> _mu2eDetector_h;

  public:
    explicit TrackRecoCheck(fhicl::ParameterSet const& pset);
    virtual ~TrackRecoCheck();

    void     Debug_11(const art::Event* Evt);   // handles fDr
//-----------------------------------------------------------------------------
// overloaded virtual methods of the base class
//-----------------------------------------------------------------------------
    virtual void     beginJob();
    virtual void     endJob  ();
    virtual void     analyzer(const art::Event& Evt);
  };


//-----------------------------------------------------------------------------
  TrackRecoCheck::TrackRecoCheck(fhicl::ParameterSet const& pset): 
    StntupleModule(pset,"TrackRecoCheck"),
    fProcessName              (pset.get<std::string>("processName", "")),
    fG4ModuleLabel            (pset.get<std::string>("g4ModuleLabel"               )),
    fTrkPatRecModuleLabel     (pset.get<std::string>("trkPatRecModuleLabel"        ))
  {

  }

//-----------------------------------------------------------------------------
  TrackRecoCheck::~TrackRecoCheck() { 
    //    delete fHist.fDr;
  }


//-----------------------------------------------------------------------------
  void TrackRecoCheck::endJob() {
    art::ServiceHandle<art::TFileService> tfs;
  }
//-----------------------------------------------------------------------------
  void TrackRecoCheck::beginJob() {

    art::ServiceHandle<art::TFileService> tfs;

    fHist.fHitDw        = tfs->make<TH1F>("hit_dw","Hit Dw",200,-1000,1000);
    fHist.fHitRes       = tfs->make<TH1F>("hit_res","Hit Residual",1000,-5,5);
    fHist.fPath         = tfs->make<TH1F>("path","Track path in a straw",200,0,20);
    fHist.fEHitVsPath   = tfs->make<TH2F>("ehit_vs_path","EHit Vs Path",200,0,20,1000,0.,0.005);
//-----------------------------------------------------------------------------
// define collection names to be used for initialization
//-----------------------------------------------------------------------------
  }


  //-----------------------------------------------------------------------------
  void TrackRecoCheck::analyzer(const art::Event& Evt) {
    const char* oname = "TrackRecoCheck::filter";

    if (DebugBit(11)) printf("[%s] RUN: %10i EVENT: %10i\n",oname,Evt.run(),Evt.event());
//-----------------------------------------------------------------------------
// do the rest
//-----------------------------------------------------------------------------
    if (DebugBit(11)) Debug_11(&Evt);
  }

//-----------------------------------------------------------------------------
// fill E(hit) vs path distributions
//-----------------------------------------------------------------------------
  void TrackRecoCheck::Debug_11(const art::Event* Evt) {
    const char* oname = "TrackRecoCheck::Debug_11";

    //     GeomHandle<Mu2eDetectorModel> detmodel;
    //    Mu2eDetectorModel const& detmodel { art::ServiceHandle<BTrkHelper>()->detectorModel() };

    auto detmodel = _mu2eDetector_h.getPtr(Evt->id());

    art::Handle<mu2e::KalRepPtrCollection> krepsHandle;

    art::Selector  selector(art::ProcessNameSelector(fProcessName)         && 
			    art::ModuleLabelSelector(fTrkPatRecModuleLabel)   );
    Evt->get(selector,krepsHandle);
//-----------------------------------------------------------------------------
// make sure collection exists
//-----------------------------------------------------------------------------
    if (! krepsHandle.isValid()) {
      printf("ERROR in %s: no KalRepPtrCollection for module %s, BAIL OUT\n",
	     oname, fTrkPatRecModuleLabel.data());
      return;
    }

    int ntrk = krepsHandle->size();

    if (DebugBit(1)) {
      printf("%s: ntrk = %i\n",oname, ntrk);
    }

    const KalRep             *trk;
    HepPoint                 tpos;
    Hep3Vector               hpos;
    const TrkHitVector*      hot_list;
    const mu2e::TrkStrawHit* hit;
    const mu2e::ComboHit*    sh;

    double   res, resErr, len, dx, dy, dw, path, ehit; // , sigw;

    //    int banner_printed = 0;
    for (int i=0; i<ntrk; i++) {
      trk = krepsHandle->at(i).get();

//       double mom   = trk->momentum().mag();
//       double pt    = trk->momentum().perp();
//       double costh = trk->momentum().cosTheta();
//       double chi2  = trk->chisq();
//       int    ndof  = trk->nDof ();
//       int    nact  = trk->nActive();
//       double t0    = trk->t0().t0();
//       double fit_consistency = trk->chisqConsistency().consistency();
//       int q        = trk->charge();

      hot_list = &trk->hitVector();

      for(auto it=hot_list->begin(); it!=hot_list->end(); it++) {
      // TrkStrawHit inherits from TrkHitOnTrk

	hit  = (const mu2e::TrkStrawHit*) (*it);
	sh   = &hit->comboHit();

	const mu2e::Straw* straw     = &hit->straw();
	const CLHEP::Hep3Vector* dir = &straw->direction();

					// gasPath() returns the half-path
	const DetStrawElem* strawelem = detmodel->strawElem(*straw);

	path = 2.*strawelem->gasPath(hit->driftRadius(),hit->trkTraj()->direction( hit->fltLen()));
	ehit = sh->energyDep();

	//	const mu2e::StrawHit* sh = &hit->strawHit();
	//	mu2e::Straw*   straw = (mu2e::Straw*) &hit->straw();

	hit->resid(res,resErr,true);
	hit->hitPosition(hpos);

	len  = hit->fltLen();
	tpos = trk->position(len);

	dx = hpos.x()-tpos.x();
	dy = hpos.y()-tpos.y();

	dw = sqrt(dx*dx+dy*dy);

	double dww = dx*dir->x() +dy*dir->y();

	if (DebugBit(1) == 1) {
	  printf("dx,dy,dw,len,dww : %10.3f %10.3f %10.3f %10.3f %10.3f\n",dx,dy,dw,len,dww);
	}

	//	sigw  = hitpos->posRes(StrawHitPosition::phi); 

	fHist.fHitDw->Fill(dww);
	fHist.fHitRes->Fill(res);
	fHist.fPath->Fill(path);
	fHist.fEHitVsPath->Fill(path,ehit);
      }
    }
  }

}

using mu2e::TrackRecoCheck;
DEFINE_ART_MODULE(TrackRecoCheck)
