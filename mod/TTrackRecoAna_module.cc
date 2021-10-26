/*
  Analyzes events to determine what the natural logarithm of the ratio
  of the sums of the squares of the residuals of pairs of doublets in
  muons and electrons are. Created in order to help determine how to
  reject muons with momentums ~105MeV/c which are reconstructed as
  having high likelihoods of being electrons.
*/

#include "mod/TTrackRecoAna_module.hh"

#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Selector.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/EDAnalyzer.h"  //Where are these files?
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"

#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"

#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/VirtualDetector.hh"

// #include "TrackerGeom/inc/Tracker.hh"

#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/PIDProduct.hh"

#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/StepPointMC.hh"
// #include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"

#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/KalmanTrack/KalHit.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"

#include "BTrkData/inc/TrkStrawHit.hh"
#include "RecoDataProducts/inc/KalRepPtrCollection.hh"

#include "Stntuple/print/TAnaDump.hh"

// C++ includes.
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>

//Histogram includes
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TEnv.h"
#include "TSystem.h"
#include "TROOT.h"


using namespace std;

namespace mu2e {
  
  TTrackRecoAna::TTrackRecoAna(fhicl::ParameterSet const& pset) :
    EDAnalyzer(pset)
    , fDiagLevel     (pset.get<int>        ("diagLevel"            ))
    , fTrkPatRecLabel(pset.get<std::string>("trkPatRecLabel"       ))
  { 
    fEvent = 0;   
  }


//-----------------------------------------------------------------------------
  void TTrackRecoAna::BookHistograms() {
    art::ServiceHandle<art::TFileService> tfs;

    for (int i=0; i<fNHistSets; i++) {
      
      fHist[i].fNDoublets   = tfs->make<TH1F>(Form("ndoub_%i"  ,i),Form("N(Doublets) [%i]"     ,i),40, 0., 40.);
      fHist[i].fNHits       = tfs->make<TH1F>(Form("nhits_%i"  ,i),Form("Nhits [%i]"           ,i),125, 0., 125.);
      fHist[i].fChi2T       = tfs->make<TH1F>(Form("chi2t_%i",i),
					 Form("Chi2 Track/NDof [%i]",i),200,0.,20);

      fHist[i].fChi2D       = tfs->make<TH1F>(Form("chi2d_%i",i),Form("Chi2D/ND [%i]",i),200,0.,20);

      fHist[i].fSlope       = tfs->make<TH1F>(Form("sl_%i",i),Form("Slope [%i]",i),200,-1.e-3,1.e-3);
      
    }
  }
  

//-----------------------------------------------------------------------------
  void TTrackRecoAna::FillHistograms(Hist_t* Hist) {

    Hist->fNDoublets->Fill(fEle.fNDoublets);
    Hist->fNHits->Fill(fEle.fNHits);

    Hist->fChi2T->Fill(fEle.fChi2T);
    Hist->fChi2D->Fill(fEle.fChi2);
    Hist->fSlope->Fill(fEle.fSlope);
  }

//-----------------------------------------------------------------------------
  void TTrackRecoAna::beginJob() {
    BookHistograms();
  }
  
//-----------------------------------------------------------------------------
  void TTrackRecoAna::beginRun(const art::Run& run) {
  }

 
//______________________________________________________________________________
  TTrackRecoAna::~TTrackRecoAna() {
//     delete fEle;
//     delete fMuo;
  }



//-----------------------------------------------------------------------------
  int  TTrackRecoAna::MakeDoublets(const KalRep* Track, DoubletRecoPar_t* Tp) {

    int      resgood[500], sector[500], device[500], layer[500], number[500], iambig[500], active[500];
    int      dn, nhits(0);
    double   hitres, hiterr;
   
    const TrkHitVector* hot_list = &Track->hitVector();
//-----------------------------------------------------------------------------
// count number of hits on the track
//-----------------------------------------------------------------------------
    const mu2e::TrkStrawHit* hit;
    mu2e::Straw*             straw;

    for (auto it=hot_list->begin(); it<hot_list->end(); it++) {
      hit   = (const mu2e::TrkStrawHit*) (*it);
      straw = (mu2e::Straw*) &hit->straw();
      
      sector [nhits] = straw->id().getPanel();
      device [nhits] = straw->id().getPlane();
      layer  [nhits] = straw->id().getLayer();
      number [nhits] = straw->id().getStraw();
      iambig [nhits] = hit->ambig();
      active [nhits] = hit->isActive();
      resgood[nhits] = hit->resid(hitres,hiterr,0); // keep the hit because it is a part of a doublet 
      
      Tp->fStrawHit[nhits] = hit;
      nhits += 1;
    }

    Tp->fTrack = Track;
    Tp->fNHits = nhits;
    Tp->fChi2T = Track->chisq()/(Track->nActive()-4.9999);

    int   hitsactive = Track->nActive();
    float costh      = Track->momentum().cosTheta();
    float fitcons    = Track->chisqConsistency().consistency();

    float h1_fltlen  = Track->firstHit()->kalHit()->hit()->fltLen() - 10;
    float hn_fltlen  = Track->lastHit ()->kalHit()->hit()->fltLen() - 10;
    float entlen     = std::min(h1_fltlen,hn_fltlen);
    
    BbrVectorErr      momerr = Track->momentumErr(entlen);
    CLHEP::Hep3Vector fitmom = Track->momentum(entlen);
    CLHEP::Hep3Vector momdir = fitmom.unit();
    
    CLHEP::HepVector momvec(3);
    for (int i=0; i<3; i++) momvec[i] = momdir[i];

    float fitmom_err = sqrt(momerr.covMatrix().similarity(momvec));
    //    float t0         = Track->t0().t0();
    float t0err      = Track->t0().t0Err();

    Tp->fIDWord = 0;

    if (hitsactive <  25                            ) Tp->fIDWord |= 0x0001;
    if (fitcons    <= 0.002                         ) Tp->fIDWord |= 0x0002;
    if ((costh     >= (1/(sqrt(2)))) || (costh<=0.5)) Tp->fIDWord |= 0x0004;
    if (fitmom_err >= 0.25                          ) Tp->fIDWord |= 0x0008;
    if (t0err      >= 0.9                           ) Tp->fIDWord |= 0x0010;

    for(int i=0; i<nhits; i++) { 
      
      Tp->fDoublet[i].fHit  [0] = Tp->fStrawHit[i];
      Tp->fDoublet[i].fIndex[0] = i;
      Tp->fDoublet[i].fNHits    = 1;

      for(int j=0; j<i; j++) { // this loop checks current hit against all previous hits in list
	dn = number[i]-number[j];

	if ((device [i] == device[j]) &&
	    (sector [i] == sector[j]) &&
	    (layer  [i] != layer [j]) && 
	    (abs(dn)    <= 2        ) && 
	    (iambig [i] != iambig[j]) &&
	    (active [i] == 1) && (active [j] == 1) && 
	    (resgood[i] == 1) && (resgood[j] == 1)    ) {   //doublet conditions

	  Tp->fDoublet[i].fHit  [Tp->fDoublet[i].fNHits] = Tp->fStrawHit[j];
	  Tp->fDoublet[i].fIndex[Tp->fDoublet[i].fNHits] = j;
	  Tp->fDoublet[i].fNHits         += 1;
	}
      }
    }

    Tp->fNDoublets    = 0;
    Tp->fNTriplets    = 0;
    Tp->fNQuadruplets = 0;

    for (int i=0; i<nhits; i++) {
      if (Tp->fDoublet[i].fNHits > 1) {
	Tp->fDoubletPtr[Tp->fNDoublets] =  &Tp->fDoublet[i];
	Tp->fNDoublets += 1;
	if (Tp->fDoublet[i].fNHits > 2) {
	  printf(" got a triplet hit number i = %5i! \n",i);
	}
      }
    }
//-----------------------------------------------------------------------------
// calculate chi2 and slope ds/dz
//-----------------------------------------------------------------------------
    CLHEP::Hep3Vector p0, p1;
    double            s12, z, r0, r1;
    Doublet_t*        d;

    double            sx(0), sy(0), sx2(0), sxy(0), sy2(0), sigxx(0), sigxy(0); //, sigyy(0);
    double            xmean, ymean, x2mean, xymean; //, y2mean;

    Tp->fChi2 = 0;
    for (int i=0; i<Tp->fNDoublets; i++) {
      d = Tp->fDoubletPtr[i];
      d->fHit[0]->resid(r0,hiterr,0);
      d->fHit[1]->resid(r1,hiterr,0);
      s12 = r0*d->fHit[0]->ambig()+r1*d->fHit[1]->ambig();
	
      d->fHit[0]->hitPosition(p0);
      d->fHit[1]->hitPosition(p1);
	
      z = (p0.z()+p1.z())/2.;
					// assume 100 microns resolution
      Tp->fChi2 += s12*s12/2/(0.1*0.1);

      sx  += z;
      sy  += s12;
      sx2 += z*z;
      sy2 += s12*s12;
      sxy += s12*z;
    }

    if (Tp->fNDoublets < 2) {
      Tp->fSlope = -1.e6;
    }
    else {
      xmean    = sx /Tp->fNDoublets;
      ymean    = sy /Tp->fNDoublets;
      x2mean   = sx2/Tp->fNDoublets;
      xymean   = sxy/Tp->fNDoublets;
      //      y2mean   = sy2/Tp->fNDoublets;
      
      sigxx = x2mean-xmean*xmean;
      sigxy = xymean-xmean*ymean;
      //      sigyy = y2mean-ymean*ymean;
      
      Tp->fChi2  = Tp->fChi2/Tp->fNDoublets;
      Tp->fSlope = sigxy/sigxx;
    }
//-----------------------------------------------------------------------------
// diagnostics printout: sum of doublet residuals vs Z
//-----------------------------------------------------------------------------
    if (fDiagLevel > 0) {
      printf(" print doublets: \n");

      for (int i=0; i<Tp->fNDoublets; i++) {
	d = Tp->fDoubletPtr[i];
	d->fHit[0]->resid(r0,hiterr,0);
	d->fHit[1]->resid(r1,hiterr,0);
	s12 = r0*d->fHit[0]->ambig()+r1*d->fHit[1]->ambig();
	
	d->fHit[0]->hitPosition(p0);
	d->fHit[1]->hitPosition(p1);
	
	z = (p0.z()+p1.z())/2.;
	
	printf("doublet: i= %3i z = %10.3f sum(res*ambig) = %10.3f \n",i,z,s12);

      }
      printf(" chi2/N = %10.3f, slope = %12.5e\n",Tp->fChi2, Tp->fSlope);
    }

    return 0;
  }

//-----------------------------------------------------------------------------
  void TTrackRecoAna::doubletmaker(const KalRep* ele_Trk) {
        
    art::Handle<mu2e::KalRepPtrCollection> eleKrepsHandle;
    art::Handle<mu2e::KalRepPtrCollection> muoKrepsHandle;
    
//     int gen_index;
//     int generatorID = 28;
//     //    int simID = 1;
//     int pdgID = 11;
//     //    float  t0true = -1.e-12;
//     static char   g4_module_label  [100] = "g4run";

//     //-------------------------------------------------------------------
//     mu2e::GeomHandle<mu2e::VirtualDetector> vdg; 
//     if (vdg->nDet() > 0) {
//       art::Handle<mu2e::StepPointMCCollection> vdhits;
//       fEvent->getByLabel(g4_module_label,"virtualdetector",vdhits);
//       int nvdhits = vdhits->size();
//       for (int i=0; i<nvdhits; i++) {
// 	const mu2e::StepPointMC* hit = &(*vdhits)[i];
// 	mu2e::VirtualDetectorId vdid(hit->volumeId());
// 	if (vdid.isTrackerMid()) {
// 	  art::Ptr<mu2e::SimParticle> const& simptr = hit->simParticle();
// 	  const mu2e::SimParticle* sim  = simptr.operator ->();
// 	  if (sim->fromGenerator()) gen_index = sim->genParticle()->generatorId().id();
// 	  else                         gen_index = -1;
// 	  int pdg_id = sim->pdgId();
// 	  //	  int sim_id = sim->id().asInt();
// 	  if (/*(sim_id == simID) && */(gen_index == generatorID) && (pdg_id == pdgID)) {
// 	    // time folding and Adding ghost hits to properly treat boundary conditions with folding, see docdb-3425
// 	    // double hitTimeUnfolded = _toff.timeWithOffsetsApplied(*hit);
// 	    // t0true                 = fmod(hitTimeUnfolded,_mbtime);
// 	  }
// 	}
//       }  
//     }
//-----------------------------------------------------------------------------
// make doublets
//-----------------------------------------------------------------------------
  }

//-------------------------------------------------------------------------
  void TTrackRecoAna::analyze(const art::Event& Evt) {

    const char* __oname = "TTrackRecoAna::analyze";

    fEvent = &Evt;
    
    art::Handle<mu2e::KalRepPtrCollection> dem_handle;
    
    art::Selector  s_dem(art::ProcessNameSelector("")                && 
			 art::ModuleLabelSelector(fTrkPatRecLabel)    );
    fEvent->get(s_dem,dem_handle);
//-----------------------------------------------------------------------------
// make sure collection exists
//-----------------------------------------------------------------------------
    if (! dem_handle.isValid()) {
      printf("%s: no TrkPatRecDem KalRepPtrCollection, BAIL OUT\n",__oname);
      return;
    }
    
    const mu2e::KalRepPtrCollection* list_of_ele_tracks = dem_handle.product();
    
    int n_ele_trk = list_of_ele_tracks->size();

    const KalRep *trk;

    printf("TTrackRecoAna event number: %10i \n", fEvent->event());
//-----------------------------------------------------------------------------
// at this point it is guaranteed that bot, electron and muon, hypotheses
// resulted in a found track
//-----------------------------------------------------------------------------
    TAnaDump* d = TAnaDump::Instance();

    for (int i=0; i<n_ele_trk; i++) {
      trk = list_of_ele_tracks->at(i).get();

      int  rc;
      rc = MakeDoublets(trk, &fEle);
      if (rc < 0) return;
//-----------------------------------------------------------------------------
// ANALYSIS:
// 1. all tracks
//-----------------------------------------------------------------------------  
      FillHistograms(&fHist[0]);
//-----------------------------------------------------------------------------
// 2. Set C tracks
//-----------------------------------------------------------------------------  
      if (fEle.fIDWord == 0) {
	FillHistograms(&fHist[1]);
      }

      if (fDiagLevel > 0) {
	d->printKalRep(trk,"data+hits");
      }
    }
    return;
  }
  //-------------------------------------------------------------------------
}

DEFINE_ART_MODULE(mu2e::TTrackRecoAna) ;
