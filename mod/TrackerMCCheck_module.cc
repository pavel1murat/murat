///////////////////////////////////////////////////////////////////////////////
// A half-interactive 2D event display. 
//
// $Id: TrackerMCCheck_module.cc,v 1.11 2014/10/02 17:15:09 murat Exp $
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
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"

#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/VirtualDetector.hh"

#include "TTrackerGeom/inc/TTracker.hh"
#include "CalorimeterGeom/inc/VaneCalorimeter.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CaloCluster/inc/CaloClusterUtilities.hh"
#include "Mu2eUtilities/inc/SimParticlesWithHits.hh"
#include "Mu2eUtilities/inc/SortedStepPoints.hh"
#include "Mu2eUtilities/inc/TrackTool.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"

#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/VirtualDetectorId.hh"

#include "TrkBase/HelixParams.hh"
#include "TrkBase/TrkHotList.hh"
#include "KalmanTrack/KalHit.hh"

#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"

#include "KalmanTests/inc/TrkStrawHit.hh"

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

  class TrackerMCCheck : public StntupleModule {
  private:
//-----------------------------------------------------------------------------
// Module labels 
//-----------------------------------------------------------------------------
    std::string        fModuleLabel;	             // this module label
    std::string        processName_;
    std::string        fGeneratorModuleLabel; 
    std::string        fG4ModuleLabel;
    
    std::string        producerName_;
    std::string        fStrawHitMaker;
    std::string        fStrawHitPosMaker;
    std::string        fFlagBgrHitsModuleLabel;
    
    int                fPdgCode;
    int                fGeneratorCode;

					// hit flag bits which should be ON and OFF

    mu2e::StrawHitFlag         fGoodHitMask;
    mu2e::StrawHitFlag         fBadHitMask; 

    const mu2e::StrawHitCollection*              fStrawHitColl;     // 
    const mu2e::StrawHitPositionCollection*      fStrawHitPosColl;  //
    const mu2e::StrawHitFlagCollection*          fStrawHitFlagColl; // 
    const mu2e::StepPointMCCollection*           fSteps;            //
    const mu2e::PtrStepPointMCVectorCollection*  fStepPointMCVectorCollection;

    struct Hist_t {
      TH1F*   fDr;                      // radial distance between the StepPointMC and the corresponding hit
      TH1F*   fNStepsPerHit;		//
      TH1F*   fNStrawHits[2];		// same distribution, different ranges 
      TH1F*   fEHitCE;			// CE straw hit energy 
      TH1F*   fEHitMu;			// CE straw hit energy 
      TH2F*   fEHitCEVsPath;		// CE straw hit energy 
      TH2F*   fEHitMuVsPath;		// muon straw hit energy 
      TH1F*   fEHitProt;		// proton straw hit energy
      TH1F*   fEHitDelta;		// delta hit energy
      TH1F*   fMomCE;			// CE     momentum 
      TH1F*   fMomMu;			// muon   momentum 
      TH1F*   fMomProt;			// proton momentum 
      TH1F*   fMomDelta;		// delta  momentum 
      TH1F*   fNStrawHitsCE[4];		//
      TH1F*   fDt;			// time division difference
      TH1F*   fWPos;
      TH2F*   fDtVsWPos;
      TH1F*   fErrPos;			// hit pos err along the wire

      TH1F*   fPStOut;
      TH1F*   fPFront [2];
      TH1F*   fCEPitch[4];
    } ;

    Hist_t fHist;


    double   fPStOut;	      // MC particle momentum out of the stopping traget
    double   fPFront;	      // MC particle momentum on tracker entry
    double   fCePitch;

					// identify the particle of choice
  public:
    explicit TrackerMCCheck(fhicl::ParameterSet const& pset);
    virtual ~TrackerMCCheck();

    void     getData(const art::Event* Evt);
    void     Init   (art::Event* Evt);
    void     Debug_003();   // handles fDr
    void     Debug_004();
//-----------------------------------------------------------------------------
// overloaded virtual methods of the base class
//-----------------------------------------------------------------------------
    virtual void     beginJob();
    virtual void     endJob  ();
    virtual bool     filter (art::Event& Evt);
  };


//-----------------------------------------------------------------------------
  TrackerMCCheck::TrackerMCCheck(fhicl::ParameterSet const& pset): 
    StntupleModule            (pset,"TrackerMCCheck"),
    fModuleLabel              (pset.get<std::string>("module_label"                )),
    processName_              (pset.get<std::string>("processName"          ,""    )),
    fGeneratorModuleLabel     (pset.get<std::string>("generatorModuleLabel"        )),
    fG4ModuleLabel            (pset.get<std::string>("g4ModuleLabel"               )),

    fStrawHitMaker            (pset.get<std::string>("strawHitMakerModuleLabel"    )),
    fStrawHitPosMaker         (pset.get<std::string>("strawHitPosMakerModuleLabel" )),
    fFlagBgrHitsModuleLabel   (pset.get<std::string>("flagBgrHitsModuleLabel"      )),
    
    fPdgCode                  (pset.get<int>        ("pdgCode"                     )),
    fGeneratorCode            (pset.get<int>        ("generatorCode"               ))
  {

  }

//-----------------------------------------------------------------------------
  TrackerMCCheck::~TrackerMCCheck() { 
    //    delete fHist.fDr;
  }


//-----------------------------------------------------------------------------
  void TrackerMCCheck::endJob() {
    art::ServiceHandle<art::TFileService> tfs;
  }
//-----------------------------------------------------------------------------
  void TrackerMCCheck::beginJob() {

    art::ServiceHandle<art::TFileService> tfs;

//     art::TFileDirectory tfdir = tfs->mkdir( "CosmicDYB" );
//     _cosmicMultiplicityH = tfdir.make<TH1D>( "MultiplicityH", "Cosmic Multiplicity", 20, -0.5, 19.5);

    fHist.fNStrawHits[0] = tfs->make<TH1F>("nsh_0" ,"Total N straw hits [0]", 200,0,  200);
    fHist.fNStrawHits[1] = tfs->make<TH1F>("nsh_1" ,"Total N straw hits [1]",1000,0,10000);

    fHist.fNStrawHitsCE[0]  = tfs->make<TH1F>("nshce_0" ,"N(straw hits) CE [0]"       , 200,0,  200);
    fHist.fNStrawHitsCE[1]  = tfs->make<TH1F>("nshce_1" ,"N(straw hits) CE [1]"       , 200,0,  200);
    fHist.fNStrawHitsCE[2]  = tfs->make<TH1F>("nshce_2" ,"N(straw hits) CE [2]"       , 200,0,  200);
    fHist.fNStrawHitsCE[3]  = tfs->make<TH1F>("nshce_3" ,"N(straw hits) CE [3]"       , 200,0,  200);

    fHist.fPStOut        = tfs->make<TH1F>("pstout","P(STOUT)"               ,1000,0,  200);
    fHist.fPFront[0]     = tfs->make<TH1F>("pfront_0","P(front)[0]"          ,1000,0,  200);
    fHist.fPFront[1]     = tfs->make<TH1F>("pfront_1","P(front)[1]"          ,1000,0,  200);

    fHist.fCEPitch[0]    = tfs->make<TH1F>("ce_pitch_0","CE Pitch [0]"       , 200,0,  2  );
    fHist.fCEPitch[1]    = tfs->make<TH1F>("ce_pitch_1","CE Pitch [1]"       , 200,0,  2  );
    fHist.fCEPitch[2]    = tfs->make<TH1F>("ce_pitch_2","CE Pitch [2]"       , 200,0,  2  );
    fHist.fCEPitch[3]    = tfs->make<TH1F>("ce_pitch_3","CE Pitch [3]"       , 200,0,  2  );

    fHist.fNStepsPerHit = tfs->make<TH1F>("nsph","N steps per hit"          , 100,0,  100);
    fHist.fDr           = tfs->make<TH1F>("dr"  ,"Hit DR"            , 100,-100,100);
    fHist.fEHitCE       = tfs->make<TH1F>("ehce","E(hit) CE"         , 300,0,0.03);
    fHist.fEHitMu       = tfs->make<TH1F>("ehmu","E(hit) Muon"       , 300,0,0.03);
    fHist.fEHitProt     = tfs->make<TH1F>("ehpr","E(hit) Proton"     , 300,0,0.03);
    fHist.fEHitDelta    = tfs->make<TH1F>("ehdl","E(hit) Delta"      , 300,0,0.03);
    fHist.fMomCE        = tfs->make<TH1F>("pce" ,"Momentum CE"       , 220,0,110);
    fHist.fMomProt      = tfs->make<TH1F>("ppr" ,"Momentum Proton"   , 500,0,500);
    fHist.fMomDelta     = tfs->make<TH1F>("pdl" ,"Momentum Delta"    , 200,0, 10);
    fHist.fDt           = tfs->make<TH1F>("dt"    ,"Delta(T) left-right", 1000,-10,10);
    fHist.fWPos         = tfs->make<TH1F>("wpos"  ,"Hit Position along the Wire", 200,-1000,1000);
    fHist.fDtVsWPos     = tfs->make<TH2F>("dt_vs_wpos" ,"Hit Dt vs WPos", 200,-1000,1000, 200, -10,10);
    fHist.fErrPos       = tfs->make<TH1F>("errpos","Hit Pos Err along the Wire", 100,0,200);

    fHist.fEHitCEVsPath = tfs->make<TH2F>("ehce_vs_path","E(hit) CE   Vs Path" , 100,0,10,300,0,0.03);
    fHist.fEHitMuVsPath = tfs->make<TH2F>("ehmu_vs_path","E(hit) muon Vs Path" , 100,0,10,300,0,0.03);

//-----------------------------------------------------------------------------
// define collection names to be used for initialization
//-----------------------------------------------------------------------------
  }

//-----------------------------------------------------------------------------
// get data from the event record
//-----------------------------------------------------------------------------
  void TrackerMCCheck::getData(const art::Event* Evt) {
    //    int   rc (0);
    const char* oname = "TrackerMCCheck::getData";

    art::Handle<mu2e::StrawHitFlagCollection> shflagH;

    art::Handle<StepPointMCCollection> stepsHandle;
    art::Selector getTrackerSteps(art::ProductInstanceNameSelector("tracker") &&
				  art::ProcessNameSelector(processName_) &&
				  art::ModuleLabelSelector(fG4ModuleLabel)  );
    Evt->get(getTrackerSteps, stepsHandle);
    if (stepsHandle.isValid()) fSteps =  (const mu2e::StepPointMCCollection*) &(*stepsHandle);
    else                       fSteps = NULL;
//-----------------------------------------------------------------------------
//  straw hit information
//-----------------------------------------------------------------------------
    art::Handle<StrawHitCollection> shH;
    Evt->getByLabel(fStrawHitMaker,shH);
    if (shH.isValid()) fStrawHitColl = shH.product();
    else               fStrawHitColl = NULL;
    
    art::Handle<mu2e::StrawHitPositionCollection> shpH;
    Evt->getByLabel(fStrawHitPosMaker,shpH);
    if (shpH.isValid()) fStrawHitPosColl = shpH.product();
    else                fStrawHitPosColl = NULL;
//-----------------------------------------------------------------------------
// get straw hit flags (half-hack)
//-----------------------------------------------------------------------------
    Evt->getByLabel(fFlagBgrHitsModuleLabel.data(),shflagH);
    if (shflagH.isValid()) fStrawHitFlagColl = shflagH.product();
    else {
      printf(">>> ERROR in %s: Straw Hit Flag Collection by \"%s\" doesn't exist. Bail Out.\n",
	     oname,fFlagBgrHitsModuleLabel.data());
      fStrawHitFlagColl = NULL;
    }

    // 12 - 11 -2013 giani added some MC info of the straws
    art::Handle<mu2e::PtrStepPointMCVectorCollection> mcptrHandleStraw;
    Evt->getByLabel(fStrawHitMaker,"StrawHitMCPtr",mcptrHandleStraw);

    if (mcptrHandleStraw.isValid()) fStepPointMCVectorCollection = mcptrHandleStraw.product();
    else                            fStepPointMCVectorCollection = NULL;
 
//     art::Handle<mu2e::StrawHitFlagCollection> shfH;
//     Evt->getByLabel(fStrawHitFlagMaker,shfH);
//     fStrawHitFlagColl = shfH.product();
  }

//-----------------------------------------------------------------------------
  void TrackerMCCheck::Init(art::Event* Evt) {
//    TStnCluster*    cluster;
//    int             id_word, ntrk;

  }


  //-----------------------------------------------------------------------------
  bool TrackerMCCheck::filter(art::Event& Evt) {
    const char* oname = "TrackerMCCheck::filter";

    printf("[%s] RUN: %10i EVENT: %10i\n",oname,Evt.run(),Evt.event());

//-----------------------------------------------------------------------------
// get event data and initialize data blocks
//-----------------------------------------------------------------------------
    getData(&Evt);



//-----------------------------------------------------------------------------
// particle parameters at virtual detectors
//-----------------------------------------------------------------------------
    mu2e::GeomHandle<mu2e::VirtualDetector>   vdg;
    art::Handle<mu2e::StepPointMCCollection>  vdhits;
    const mu2e::StepPointMC                   *hit; 
    const mu2e::SimParticle                   *sim;
    int   pdg_id, gen_id; 

    fPFront  = -1.;
    fPStOut  = -1.;
    fCePitch = -1.;

    if (vdg->nDet() > 0) {
      Evt.getByLabel(fG4ModuleLabel,"virtualdetector",vdhits);
      if (!vdhits.isValid()) {
	printf("ERROR %s : No valid VD step points\n",oname);
      }
      else {
	int nvdhits = vdhits->size();
//-----------------------------------------------------------------------------
// loop over the virtual detector hits
//-----------------------------------------------------------------------------
	for (int i=0; i<nvdhits; i++) {
	  hit = &vdhits->at(i);
	  mu2e::VirtualDetectorId vdid(hit->volumeId());
	  sim = hit->simParticle().get();
	  if (sim == NULL) {
	    printf(">>> ERROR: %s sim == NULL\n",oname);
	                                                    goto NEXT_VD_HIT;
	  }
	  
	  pdg_id        = sim->pdgId();

	  if (sim->fromGenerator()) gen_id = sim->genParticle()->generatorId().id();
	  else                      gen_id = -1;
	  
	  if (vdid.id() == mu2e::VirtualDetectorId::ST_Out) {
	    if ((pdg_id == fPdgCode) && (gen_id == fGeneratorCode)) {
	      fPStOut = hit->momentum().mag();
	    }
	  }
	  else if (vdid.isTrackerFront()) {
	    if ((pdg_id == fPdgCode) && (gen_id == fGeneratorCode)) {
	      fPFront  = hit->momentum().mag();
	      fCePitch = hit->momentum().perp()/hit->momentum().z();
	    }
	  }
	NEXT_VD_HIT:;
	}
      }
    }



    if (DebugBit(3)) {
      Debug_003();
    }
    else if (DebugBit(4)) {
      Debug_004();
    }

    return true;
  }

//-----------------------------------------------------------------------------
  void TrackerMCCheck::Debug_004() {

//     mu2e::GeomHandle<mu2e::TTracker> ttHandle;
//     const mu2e::TTracker* tracker = ttHandle.get();

    float                  dt, p, ehit, wpos, errpos;

    int                    nhits, nhits_ce, pdg_id, mother_pdg_id, nsteps_per_hit;
    int                    gen_code; //, sim_id;

    const mu2e::StrawHit           *hit;
    const mu2e::StrawHitPosition   *shp;
    const mu2e::StepPointMC        *step; 

    static const double MIN_PITCH = 1;
    static const double MAX_PITCH = sqrt(3);

    nhits = fStrawHitColl->size();

    fHist.fNStrawHits[0]->Fill(nhits);
    fHist.fNStrawHits[1]->Fill(nhits);

    nhits_ce = 0;

    for (int i=0;  i<nhits;  i++) {
      mu2e::PtrStepPointMCVector const& mcptr(fStepPointMCVectorCollection->at(i));
      nsteps_per_hit = mcptr.size();
      step = mcptr[0].get();
    
      hit   = &fStrawHitColl->at(i);
      shp   = &fStrawHitPosColl->at(i);

      art::Ptr<mu2e::SimParticle> const& simptr = step->simParticle(); 
      art::Ptr<mu2e::SimParticle> mother = simptr;
      while(mother->hasParent())  mother = mother->parent();
      const mu2e::SimParticle *   sim    = mother.operator ->();
      
      p             = step->momentum().mag();
      ehit          = hit->energyDep();
      pdg_id        = simptr->pdgId();
      mother_pdg_id = sim->pdgId();
      dt            = hit->dt();
      //      sim_id        = simptr->id().asInt();

      if (simptr->fromGenerator()) gen_code = simptr->genParticle()->generatorId().id();
      else                         gen_code = -1;

      if ((pdg_id == fPdgCode) && (gen_code == fGeneratorCode)) {
	nhits_ce += 1;
      }

      wpos          = shp->wireDist();
      errpos        = shp->posRes(StrawHitPosition::phi);   

      fHist.fDt->Fill(dt);
      fHist.fNStepsPerHit->Fill(nsteps_per_hit);
      fHist.fWPos->Fill(wpos);
      fHist.fDt->Fill(dt);
      fHist.fDtVsWPos->Fill(wpos,dt);
      fHist.fErrPos->Fill(errpos);

      if (pdg_id == 11) {
	if (mother_pdg_id == 11) {
					// electrons, not from photon conversions
	  fHist.fEHitCE->Fill(ehit);
	  //	  fHist.fEHitCEVsPath->Fill(path,ehit);
	  fHist.fMomCE->Fill(p);
	}
	else if (mother_pdg_id == 22) {
					// delta-electrons(from photon conversions)
	  fHist.fEHitDelta->Fill(ehit);
	  fHist.fMomDelta->Fill(p);
	}
      }
      else if (pdg_id == 2212) {
					// protons
	fHist.fEHitProt->Fill(ehit);
	fHist.fMomProt->Fill(p);
      }
      else if (pdg_id == 13) {
					// protons
	fHist.fEHitMu->Fill(ehit);
	//	fHist.fEHitMuVsPath->Fill(path,ehit);
	fHist.fMomMu->Fill(p);
      }
    }
//-----------------------------------------------------------------------------
// Track efficiency cut chain
// 2. cut on N(CE hits)
//-----------------------------------------------------------------------------
    fHist.fPStOut->Fill(fPStOut);
    fHist.fNStrawHitsCE[0]->Fill(nhits_ce);
    fHist.fPFront      [0]->Fill(fPFront);
    fHist.fCEPitch     [0]->Fill(fCePitch);
    if (nhits_ce >= 20) {
      fHist.fNStrawHitsCE[1]->Fill(nhits_ce);
      fHist.fPFront      [1]->Fill(fPFront);
      fHist.fCEPitch     [1]->Fill(fCePitch);
//-----------------------------------------------------------------------------
// 2. cut on P(front)
//-----------------------------------------------------------------------------
      if (fPFront > 100.) {
	fHist.fNStrawHitsCE[2]->Fill(nhits_ce);
	fHist.fCEPitch     [2]->Fill(fCePitch);
//-----------------------------------------------------------------------------
// 3. cut on the CE pitch  30 < lambda < 45 deg
//-----------------------------------------------------------------------------
	if ((fCePitch > MIN_PITCH) && (fCePitch < MAX_PITCH)) {
	  fHist.fNStrawHitsCE[3]->Fill(nhits_ce);
	  fHist.fCEPitch     [3]->Fill(fCePitch);
	  // continue
	}
      }
    }
  }
    
//-----------------------------------------------------------------------------
// plot distribution in radial distance between the StepPointMC and the 
// corresponding hit - fHist.fDr
//-----------------------------------------------------------------------------
  void TrackerMCCheck::Debug_003() {
    
    int nsteps, nhits, jclosest;

    //    static TH1F*  hist(0);

    //    static TCanvas* c;

    const mu2e::StepPointMC        *step;
    const mu2e::StrawHitPosition   *hitp, *closest_hitp;
    //    const mu2e::StrawHit           *sh;
    //    const mu2e::Straw              *straw;

    double   zstep, zhit, dz, dz_min, /* dx, dy, drho, */ rh, rs, dr;

    //    mu2e::GeomHandle<mu2e::TTracker> ttHandle;

    //    const mu2e::TTracker* tracker = ttHandle.get();

    nsteps = fSteps->size();
    nhits  = fStrawHitPosColl->size();
    
    for (int i=0; i<nsteps; i++) {
      step =  &fSteps->at(i);

      art::Ptr<mu2e::SimParticle> const& simptr = step->simParticle();
      const mu2e::SimParticle* sim  = simptr.operator ->();
      int gen_id = (int) sim->generatorIndex();

      if (gen_id != 0) goto NEXT_STEP;

      zstep = step->position().z();

      // for each steppoint find closest StrawHit

      jclosest = -1;
      dz_min   = 1e10;

      for (int j=0; j<nhits; j++) {
	hitp   = &fStrawHitPosColl->at(j);
	zhit   = hitp->pos().z();

	dz = zstep-zhit;
	if (fabs(dz) < dz_min) {
	  dz_min   = fabs(dz);
	  jclosest = j;
	}
      }

      if (jclosest >= 0) {
	closest_hitp = &fStrawHitPosColl->at(jclosest);

	//	dx = step->position().x()-closest_hitp->pos().x();
	//	dy = step->position().y()-closest_hitp->pos().y();

	//	sh = &fStrawHitColl->at(jclosest);
	
	//      mu2e::StrawIndex index = sh->strawIndex();

	//	straw = &tracker->getStraw(sh->strawIndex());

	//	drho = dx*straw->getDirection().x()+dy*straw->getDirection().y();

	rs = sqrt(step->position().x()*step->position().x() + 
		  step->position().y()*step->position().y());

	rh = sqrt(closest_hitp->pos().x()*closest_hitp->pos().x() + 
		  closest_hitp->pos().y()*closest_hitp->pos().y());

	dr = rh-rs;

	fHist.fDr->Fill(dr);

// 	printf("MuHitDisplay::Debug_003: i,z,dz,drho,dr = %5i %10.3f %10.3f %10.3f %10.3f\n",
// 	       i,zstep,dz_min,drho,dr);
      }
      else {
	//	printf("MuHitDisplay::Debug_003: i,z : %5i %10.3f : NO CLOSEST HIT\n",i,zstep);
      }
    NEXT_STEP:;
    }

  }
    
}

using mu2e::TrackerMCCheck;
DEFINE_ART_MODULE(TrackerMCCheck);
