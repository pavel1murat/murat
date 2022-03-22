//////////////////////////////////////////////////////////////////////////////
// framework
//
// parameter defaults: CalPatRec/fcl/prolog.fcl
// this module doesn't do reconstruction
// on input, it takes a list  of ComboHitFlags flags and eveluates performance 
// of the delta electron tagging
//
// hit type = 0 : proton
//            1 : ele 0 < p < 20
//            2 : ele 20 < p < 90
//            3 : ele 90 < p < 110
//            4 : everything else
//////////////////////////////////////////////////////////////////////////////
#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"

#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
// conditions
#include "Offline/ConditionsService/inc/ConditionsHandle.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
// root 
#include "TMath.h"
#include "TH1F.h"
#include "TH1.h"
#include "TTree.h"
#include "TH2F.h"
#include "TVector2.h"
// data
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/Mu2eUtilities/inc/SimParticleTimeOffset.hh"

#include <algorithm>
#include <cmath>
#include "CLHEP/Vector/ThreeVector.h"
#include "Offline/MCDataProducts/inc/StrawGasStep.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
// ------------------------------------------------------------------------------
// #include "murat/mod/TZClusterFinder_types.hh"

using namespace std; 
using CLHEP::Hep3Vector;

namespace mu2e {
  
  class TZClusterFinderAna : public art::EDAnalyzer {
  
  public:
    enum { kNStations      = 20 };
    enum { kNFaces         =  4 };
    enum { kNPanelsPerFace =  3 };
    enum { kc2_cut         = 50 };
    enum { kPerpRes        = 10 };
    enum { kMaxDxy         = 50 };
    enum { kmax_gap        = 1  };

    enum {
      kNEventHistSets    =  10,
      kNComboHitHistSets = 100,
      kNMcHistSets       = 200
    };

  protected:

    struct McHist_t {
      TH1F*  fPDGCode;
      TH1F*  fMom;
      TH1F*  fNHits;
      TH1F*  fNHitsDelta;
      TH1F*  fFractReco;
      TH2F*  fFractRecoVsNHits;
    };

    struct ComboHitHist_t {
      TH1F*  fType;
      TH1F*  fTime;
      TH1F*  fMom;			// momentum of the particle which produced the hit
      TH1F*  fEnergyDep;
      TH1F*  fDeltaT;
      TH1F*  fPDGCode;
    };

    struct EventHist_t {
      TH1F*  fEventNumber;
      TH1F*  fRunNumber;
      TH1F*  fNSecondHits;
      TH1F*  fNSecondHitsT;
      TH1F*  fNMc;
      TH1F*  fPDGCode;
      TH1F*  fNMcHits;
      TH1F*  fNHitsDeltaT;
      TH1F*  fNHitsDeltaR;
    };

    struct Hist_t {
      EventHist_t*    fEvent   [kNEventHistSets   ];
      ComboHitHist_t* fComboHit[kNComboHitHistSets];
      McHist_t*       fMc      [kNMcHistSets      ];
    }; 

    Hist_t  _hist;

//-----------------------------------------------------------------------------
// diagnostics structures
//-----------------------------------------------------------------------------
    struct McPart_t {
      int   fFirstStation;
      int   fLastStation;
      int   fNHitsDelta;	    // number of hits associated with all reconstructed delta electrons 
      int   fTime;                  // lowest out of the hit times

      const SimParticle*               fSim;
      std::vector<const ComboHit*>     fListOfHits;
      //      std::vector<const ComboHitFlag*> fListOfFlags;

      McPart_t(const SimParticle* Sim = NULL) { 
	fSim          = Sim; 
	fFirstStation = 999;
	fLastStation  = -1;
	fNHitsDelta   = 0;
	fTime         = 1.e6;
      }

      ~McPart_t() { 
      }

      int NHits()      const { return fListOfHits.size(); }
      int NHitsDelta() const { return fNHitsDelta;        }

      float Momentum() const { 
	float px = fSim->startMomentum().px();
	float py = fSim->startMomentum().py();
	float pz = fSim->startMomentum().pz();
	return sqrt(px*px+py*py+pz*pz);
      }

      float Time() const { return fTime; }

    }; 

    struct McHitInfo_t {
      const  McPart_t*     fMc;
      int                  fType;   // 0:p, 1:ele p<20, 2:ele 20<p<80  3:ele 100<p<110 4:everything else
      StrawHitFlag         fFlag;
    };

//-----------------------------------------------------------------------------
// NStations stations, 4-1=3 faces (for hit w/ lower z), 3 panels (for hit w/ lower z)
// 2017-07-27 P.Murat: the 2nd dimension should be 3, right? 
//-----------------------------------------------------------------------------
    std::vector<McPart_t*>      _list_of_mc_particles; // list_of_particles with hits in the tracker
    std::vector<McHitInfo_t>    _list_of_mc_hit_info ; // for each straw hit, pointer to the MC info
//-----------------------------------------------------------------------------
// talk-to parameters
//-----------------------------------------------------------------------------
    art::InputTag             _chCollTag;
    art::InputTag             _mcdigiCollTag;
    int                       _debugLevel;
    int                       _diagLevel;
    int                       _printElectrons;         //
    int                       _printElectronsMinNHits;
    float                     _printElectronsMaxFReco;
    //    float                     _maxElectronHitEnergy;
//-----------------------------------------------------------------------------
// cache of event or geometry objects
//-----------------------------------------------------------------------------
    const ComboHitCollection*                   _chcol;
    const StrawDigiMCCollection*                _mcdigis;

    const Tracker*                              _tracker;
    int                                         _eventNum;
    int                                         _nsh;

    int                                         fNHitsDeltaTot;
    int                                         fNHitsDeltaReco;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  public:
    explicit     TZClusterFinderAna(fhicl::ParameterSet const&);
    virtual      ~TZClusterFinderAna();

    void bookEventHistograms   (EventHist_t*    Hist, int HistSet, art::TFileDirectory* Dir);
    void bookComboHitHistograms(ComboHitHist_t* Hist, int HistSet, art::TFileDirectory* Dir);
    void bookMcHistograms      (McHist_t*       Hist, int HistSet, art::TFileDirectory* Dir);

    void bookHistograms();

    void fillEventHistograms   (EventHist_t*    Hist);
    void fillMcHistograms      (McHist_t*       Hist, McPart_t* Mc );
    void fillComboHitHistograms(ComboHitHist_t* Hist, const ComboHit* Hit, McHitInfo_t* McHitInfo);

    void fillHistograms();

    void OrderHits();
    void debug();
//-----------------------------------------------------------------------------
// overloaded methods of the base class
//-----------------------------------------------------------------------------
    virtual void beginJob()                     override;
    virtual void beginRun(const art::Run& ARun) override;
    virtual void analyze (const art::Event& e ) override;

//-----------------------------------------------------------------------------
// do these need to be private ?
//-----------------------------------------------------------------------------
    bool      findData     (const art::Event&  Evt);
    McPart_t* findParticle (const SimParticle* Sim);
    int       initMcDiag      ();
    int       associateMcTruth();
  };

//-----------------------------------------------------------------------------
  TZClusterFinderAna::TZClusterFinderAna(fhicl::ParameterSet const& pset): 
    art::EDAnalyzer(pset), 
    _chCollTag             (pset.get<string>       ("chCollTag"             )),
    _mcdigiCollTag         (pset.get<art::InputTag>("mcDigiCollTag"         )),
    _debugLevel            (pset.get<int>          ("debugLevel"            )),
    _diagLevel             (pset.get<int>          ("diagLevel"             )),
    _printElectrons        (pset.get<int>          ("printElectrons"        )),
    _printElectronsMinNHits(pset.get<int>          ("printElectronsMinNHits")),
    _printElectronsMaxFReco(pset.get<float>        ("printElectronsMaxFReco"))
    //    _maxElectronHitEnergy  (pset.get<float>        ("maxElectronHitEnergy"         ))
  {
  }

  TZClusterFinderAna::~TZClusterFinderAna() {
  }


//-----------------------------------------------------------------------------
  void TZClusterFinderAna::bookEventHistograms(EventHist_t* Hist, int HistSet, art::TFileDirectory* Dir) {
    Hist->fEventNumber     = Dir->make<TH1F>(Form("event_%02i", HistSet), "Event Number", 100, 0., 100000.);

    Hist->fNSecondHits     = Dir->make<TH1F>(Form("nhit2_%02i" ,HistSet), "N(second hits)", 100, 0., 100.);
    Hist->fNSecondHitsT    = Dir->make<TH1F>(Form("nhit2t_%02i",HistSet), "N(second hits) w/ time selection", 100, 0., 100.);

    Hist->fNMc             = Dir->make<TH1F>(Form("nmc_%02i"      ,HistSet), "N(MC particles)", 100, 0., 1000.);
    Hist->fPDGCode         = Dir->make<TH1F>(Form("pdg_code_%02i" ,HistSet), "PDG Code"       , 100, 0., 100.);
    Hist->fNMcHits         = Dir->make<TH1F>(Form("n_mc_hits_%02i",HistSet), "N(MC hits)"     , 250, 0., 500.);
    Hist->fNHitsDeltaT     = Dir->make<TH1F>(Form("n_delta_ht_%02i" ,HistSet), "N(delta hits T)", 500, 0., 5000.);
    Hist->fNHitsDeltaR     = Dir->make<TH1F>(Form("n_delta_hr_%02i" ,HistSet), "N(delta hits R)", 500, 0., 5000.);
  }

//-----------------------------------------------------------------------------
  void TZClusterFinderAna::bookMcHistograms(McHist_t* Hist, int HistSet, art::TFileDirectory* Dir) {

    Hist->fPDGCode    = Dir->make<TH1F>("pdg"  , "PDG code"        , 500, -250., 250.);
    Hist->fMom        = Dir->make<TH1F>("mom"  , "momentum"        , 200, 0., 200.);
    Hist->fNHits      = Dir->make<TH1F>("nhits", "N(hits)"         , 200, 0., 200.);
    Hist->fNHitsDelta = Dir->make<TH1F>("nhitsr", "N(hits reco)"   , 200, 0., 200.);
    Hist->fFractReco  = Dir->make<TH1F>("fractr", "NR/N"           , 100, 0.,   1.);

    Hist->fFractRecoVsNHits = Dir->make<TH2F>("freco_vs_nhits", "F(Reco) vs nhits", 100, 0., 200.,100,0,1);
  }


//-----------------------------------------------------------------------------
  void TZClusterFinderAna::bookComboHitHistograms(ComboHitHist_t* Hist, int HistSet, art::TFileDirectory* Dir) {

    Hist->fType      = Dir->make<TH1F>("type", "Hit type"        , 10, 0., 10.);
    Hist->fTime      = Dir->make<TH1F>("time", "time"            , 400, 0., 2000.);
    Hist->fMom       = Dir->make<TH1F>("mom" , "Momentum"        , 200, 0., 400.);
    Hist->fEnergyDep = Dir->make<TH1F>("edep", "edep"            , 200, 0., 2e-2);
    Hist->fDeltaT    = Dir->make<TH1F>("dt"  , "DeltaT"          , 200, -10,10);
    Hist->fPDGCode   = Dir->make<TH1F>("pdg" , "PDG code"        , 2000, -10000,10000);
  }


//-----------------------------------------------------------------------------
  void TZClusterFinderAna::bookHistograms() {
//-----------------------------------------------------------------------------
// book event histograms
//-----------------------------------------------------------------------------
    art::ServiceHandle<art::TFileService> tfs;
    char   folder_name[200];

    TH1::AddDirectory(0);

    int book_event_histset[kNEventHistSets];
    for (int i=0; i<kNEventHistSets; i++) book_event_histset[i] = 0;

    book_event_histset[ 0] = 1;		// all events

    for (int i=0; i<kNEventHistSets; i++) {
      if (book_event_histset[i] != 0) {
	sprintf(folder_name,"evt_%i",i);
	art::TFileDirectory tfdir = tfs->mkdir(folder_name);

	_hist.fEvent[i] = new EventHist_t;
	bookEventHistograms(_hist.fEvent[i],i,&tfdir);
      }
    }
//-----------------------------------------------------------------------------
// book straw hit histograms
//-----------------------------------------------------------------------------
    int book_straw_hit_histset[kNComboHitHistSets];
    for (int i=0; i<kNComboHitHistSets; i++) book_straw_hit_histset[i] = 0;

    book_straw_hit_histset[  0] = 1;		// all hits
    book_straw_hit_histset[  1] = 1;		// all hits t>500
    book_straw_hit_histset[  2] = 1;		// all hits t>500 
    book_straw_hit_histset[  3] = 1;		// all hits t>500 edepOK
    book_straw_hit_histset[  4] = 1;		// all hits t>500 edepOK non-delta
    book_straw_hit_histset[  5] = 1;		// all hits t>500 edepOK non-delta MC=delta (type 1)
    book_straw_hit_histset[  6] = 1;		// all hits t>500 edepOK delta
    book_straw_hit_histset[  7] = 1;		// all hits t>500 edepOK delta MC=delta (type 1)

    book_straw_hit_histset[ 10] = 1;		// all hits type=0
    book_straw_hit_histset[ 20] = 1;		// all hits type=1
    book_straw_hit_histset[ 30] = 1;		// all hits type=2
    book_straw_hit_histset[ 40] = 1;		// all hits type=3
    book_straw_hit_histset[ 50] = 1;		// all hits type=3

    for (int i=0; i<kNComboHitHistSets; i++) {
      if (book_straw_hit_histset[i] != 0) {
	sprintf(folder_name,"sh_%i",i);
	art::TFileDirectory tfdir = tfs->mkdir(folder_name);

	_hist.fComboHit[i] = new ComboHitHist_t;
	bookComboHitHistograms(_hist.fComboHit[i],i,&tfdir);
      }
    }
//-----------------------------------------------------------------------------
// book MC histograms
//-----------------------------------------------------------------------------
    int book_mc_histset[kNMcHistSets];
    for (int i=0; i<kNMcHistSets; i++) book_mc_histset[i] = 0;

    book_mc_histset[  0] = 1;		// all particles
    book_mc_histset[  1] = 1;		// electrons
    book_mc_histset[  2] = 1;		// electrons fTime > 550
    book_mc_histset[  3] = 1;		// electrons fTime > 550 with last>first
    book_mc_histset[  4] = 1;		// electrons fTime > 550 with last>first and 6+ hits
    book_mc_histset[  5] = 1;		// electrons fTime > 550 with last>first, 5+ hits, and reco delta
    book_mc_histset[  6] = 1;		// electrons fTime > 550 with last>first, 5+ hits and p < 20

    book_mc_histset[100] = 1;		// electrons fTime > 550 and 20 < p < 80 MeV/c
    book_mc_histset[101] = 1;		// electrons fTime > 550 and p > 80 MeV/c

    for (int i=0; i<kNMcHistSets; i++) {
      if (book_mc_histset[i] != 0) {
	sprintf(folder_name,"mc_%i",i);
	art::TFileDirectory tfdir = tfs->mkdir(folder_name);

	_hist.fMc[i] = new McHist_t;
	bookMcHistograms(_hist.fMc[i],i,&tfdir);
      }
    }
  }


//-----------------------------------------------------------------------------
  void TZClusterFinderAna::beginJob() {
    bookHistograms();
  }


 
//----Get data------------------------------------------------------------------------------------------------
  void TZClusterFinderAna::beginRun(const art::Run& aRun) {

    mu2e::GeomHandle<mu2e::Tracker> ttHandle;
    _tracker = ttHandle.get();

  }




//-----------------------------------------------------------------------------
  void  TZClusterFinderAna::fillEventHistograms(EventHist_t* Hist) {
    Hist->fEventNumber->Fill(_eventNum);

//-----------------------------------------------------------------------------
// fill MC particle histograms
//-----------------------------------------------------------------------------
    int nmc = _list_of_mc_particles.size();
    for (int i=0; i<nmc; i++) {
      McPart_t* mc = _list_of_mc_particles.at(i);
      int n_mc_hits = mc->fListOfHits.size();

      Hist->fNMcHits->Fill(n_mc_hits);
      Hist->fPDGCode->Fill(mc->fSim->pdgId());
    }

    Hist->fNMc->Fill(nmc);
    Hist->fNHitsDeltaT->Fill(fNHitsDeltaTot);
    Hist->fNHitsDeltaR->Fill(fNHitsDeltaReco);
  }


//-----------------------------------------------------------------------------
  void  TZClusterFinderAna::fillComboHitHistograms(ComboHitHist_t* Hist, const ComboHit* Hit, McHitInfo_t* McHitInfo) {

    const McPart_t* mc = McHitInfo->fMc;

    Hist->fType->Fill(McHitInfo->fType);
    Hist->fTime->Fill(Hit->time());
    Hist->fMom->Fill(mc->Momentum());
    Hist->fEnergyDep->Fill(Hit->energyDep());
    // Hist->fDeltaT->Fill(Hit->dt());
    Hist->fPDGCode->Fill(mc->fSim->pdgId());
  }

//-----------------------------------------------------------------------------
  void TZClusterFinderAna::fillMcHistograms(McHist_t* Hist, McPart_t* Mc) {
    const SimParticle* sim = Mc->fSim;
    float mom = Mc->Momentum();

    Hist->fPDGCode->Fill(sim->pdgId());
    Hist->fMom->Fill(mom);
    Hist->fNHits->Fill(Mc->NHits());
    Hist->fNHitsDelta->Fill(Mc->fNHitsDelta);

    float freco = Mc->fNHitsDelta/(Mc->NHits()+1.e-4);

    Hist->fFractReco->Fill(freco);

    Hist->fFractRecoVsNHits->Fill(Mc->NHits(),freco);
  }

//-----------------------------------------------------------------------------
// fill histograms
//-----------------------------------------------------------------------------
  void  TZClusterFinderAna::fillHistograms() {
//-----------------------------------------------------------------------------
// event histograms - just one set
//-----------------------------------------------------------------------------
    fillEventHistograms(_hist.fEvent[0]);
//-----------------------------------------------------------------------------
// straw hit histograms
//-----------------------------------------------------------------------------
    for (int i=0; i<_nsh; i++) {
      const ComboHit* ch          = &_chcol->at(i);
      McHitInfo_t*    mc_hit_info = &_list_of_mc_hit_info.at(i);
      
      fillComboHitHistograms(_hist.fComboHit[0],ch,mc_hit_info);
      if (ch->time() > 500) {
	fillComboHitHistograms(_hist.fComboHit[1],ch,mc_hit_info);
	
	const StrawHitFlag* flag = &ch->flag();

	fillComboHitHistograms(_hist.fComboHit[2],ch,mc_hit_info);
	int edepOK = flag->hasAllProperties(StrawHitFlag::energysel);
	if (edepOK) {
	  fillComboHitHistograms(_hist.fComboHit[3],ch,mc_hit_info);
	  int delta = flag->hasAllProperties(StrawHitFlag::bkg);
	  if (! delta) {
//-----------------------------------------------------------------------------
// ComboHit SET 4: hits not marked as delta electron hits
//          SET 5: hits of low energy electrons not marked as delta electron hits
//-----------------------------------------------------------------------------
	    fillComboHitHistograms(_hist.fComboHit[4],ch,mc_hit_info);
	    if (mc_hit_info->fType == 1) { // low-energy electrons
	      fillComboHitHistograms(_hist.fComboHit[5],ch,mc_hit_info);
	    }
	  }
	  else {
//-----------------------------------------------------------------------------
// ComboHit SET 6: hits marked as delta electron hits
//          SET 7: hits of low energy electrons marked as such
//-----------------------------------------------------------------------------
	    fillComboHitHistograms(_hist.fComboHit[6],ch,mc_hit_info);
	    if (mc_hit_info->fType == 1) { // low-energy electrons
	      fillComboHitHistograms(_hist.fComboHit[7],ch,mc_hit_info);
	    }
	  }
	}
      }

      if (mc_hit_info->fType == 0) fillComboHitHistograms(_hist.fComboHit[10],ch,mc_hit_info);
      if (mc_hit_info->fType == 1) fillComboHitHistograms(_hist.fComboHit[20],ch,mc_hit_info);
      if (mc_hit_info->fType == 2) fillComboHitHistograms(_hist.fComboHit[30],ch,mc_hit_info);
      if (mc_hit_info->fType == 3) fillComboHitHistograms(_hist.fComboHit[40],ch,mc_hit_info);
      if (mc_hit_info->fType == 4) fillComboHitHistograms(_hist.fComboHit[50],ch,mc_hit_info);
    }
//-----------------------------------------------------------------------------
// fill MC histograms
// for each delta electron, need to check which fraction of its hits has not been
// Associated with found DeltaCandidate's
//-----------------------------------------------------------------------------
    int nmc = _list_of_mc_particles.size();
    
    for (int i=0; i<nmc; i++) {
      McPart_t* mc = _list_of_mc_particles.at(i);
      const SimParticle* sim = mc->fSim;
      
      fillMcHistograms(_hist.fMc[0],mc);
//-----------------------------------------------------------------------------
// set 1: electrons
//-----------------------------------------------------------------------------
      if (sim->pdgId() == 11) {
	fillMcHistograms(_hist.fMc[1],mc);
	if (mc->Time() > 550) {
	  fillMcHistograms(_hist.fMc[2],mc);
	  if (mc->fLastStation > mc->fFirstStation) {
	    fillMcHistograms(_hist.fMc[3],mc);
	    if (mc->NHits() >= 5) {
	      fillMcHistograms(_hist.fMc[4],mc);

	      if (mc->Momentum() < 20) fillMcHistograms(_hist.fMc[6],mc);
	    }
//-----------------------------------------------------------------------------
// a closer look at misreconstructed delta electrons
//-----------------------------------------------------------------------------
	    float fr = mc->fNHitsDelta/(mc->NHits()+1.e-3);

	    if ((mc->Momentum() < 5) && (mc->Time() > 550) && (mc->NHits() > 40) && (fr < 0.5)) {
	      printf(" event: %6i missed delta: sim.id = %10li mom = %10.3f time= %9.3f nhits = %3i nhits(delta): %3i first: %2i last: %2i",
		     _eventNum,
		     sim->id().asInt(), mc->Momentum(), mc->Time(), 
		     mc->NHits(), mc->fNHitsDelta, 
		     mc->fFirstStation, mc->fLastStation);
	      printf(" fraction: %6.3f\n",fr);
	    }
	  }

	  if ((mc->Momentum() > 20) && (mc->Momentum() < 80)) fillMcHistograms(_hist.fMc[100],mc);
	  if (mc->Momentum()  > 80)                           fillMcHistograms(_hist.fMc[101],mc);
	}
      }
    }
  }
  
//-----------------------------------------------------------------------------
  TZClusterFinderAna::McPart_t* TZClusterFinderAna::findParticle(const SimParticle* Sim) {
    McPart_t* found(0);

    int n = _list_of_mc_particles.size();
    
    for (int i=0; i<n; i++) {
      McPart_t* mc = _list_of_mc_particles.at(i);
      if (mc->fSim == Sim) {
	found = mc;
	break;
      }
    }

    return found;
  }

//-----------------------------------------------------------------------------
// create a list of MC particles with hits in the tracker - hope, it is shorter
// than the list of all particles
//-----------------------------------------------------------------------------
  int TZClusterFinderAna::initMcDiag() {

    int n = _list_of_mc_particles.size();
    for (int i=0; i<n; i++) {
      McPart_t* p = _list_of_mc_particles.at(i);
      delete p;
    }

    _list_of_mc_particles.clear();

    _list_of_mc_hit_info.clear();

    _list_of_mc_hit_info.resize(_nsh);

    int delta_nhits_tot = 0;

    for (int i=0; i<_nsh; i++) {
      const ComboHit*             ch  = &_chcol->at(i);
      const mu2e::StrawDigiMC* mcdigi = &_mcdigis->at(i);

      const mu2e::StrawGasStep   *stmc;
      if (mcdigi->wireEndTime(mu2e::StrawEnd::cal) < mcdigi->wireEndTime(mu2e::StrawEnd::hv)) {
	stmc = mcdigi->strawGasStep(mu2e::StrawEnd::cal).get();
      }
      else {
	stmc = mcdigi->strawGasStep(mu2e::StrawEnd::hv ).get();
      }

      const mu2e::SimParticle* sim = &(*stmc->simParticle());

//-----------------------------------------------------------------------------
// search if this particle has already been registered
//-----------------------------------------------------------------------------
      McPart_t* mc = findParticle(sim);

      if (mc == NULL) {
					// add new particle
	mc = new McPart_t(sim);
	_list_of_mc_particles.push_back(mc);
      }
      mc->fListOfHits.push_back(ch);

      StrawId   shid  = ch->strawId();
      const Straw& straw = _tracker->getStraw(shid);
      int station = straw.id().getStation();
      if (station < mc->fFirstStation) mc->fFirstStation = station;
      if (station > mc->fLastStation ) mc->fLastStation  = station;

      if (ch->time() < mc->fTime) mc->fTime = ch->time();

      McHitInfo_t* mc_hit_info = &_list_of_mc_hit_info.at(i);

      mc_hit_info->fMc   = mc;
      mc_hit_info->fFlag = (StrawHitFlag) ch->flag();

      int pdg_id = mc->fSim->pdgId();

      if      (pdg_id == 2212)   mc_hit_info->fType = 0;
      else if (pdg_id == 11  ) { 
	float mom = mc->Momentum();
	if      (mom <  20)      {
	  mc_hit_info->fType = 1;
	  delta_nhits_tot++;
	}
	else if (mom <  90)      mc_hit_info->fType = 2;
	else if (mom < 110)      mc_hit_info->fType = 3;
      }
      else                       mc_hit_info->fType = 4;

      // int flagged_as_delta = shf->hasAnyProperty(ComboHitFlag::bkg);
	
      // if (flagged_as_delta) fNHitsDeltaReco++;
    }

    return 0;
  }

//-----------------------------------------------------------------------------
// form a list of MC particles with hits in the tracker
//-----------------------------------------------------------------------------
  int TZClusterFinderAna::associateMcTruth() {
//-----------------------------------------------------------------------------
// for each MC electron calculate the number of reconstructed hits
//-----------------------------------------------------------------------------
    int nmc    = _list_of_mc_particles.size();

    StrawHitFlag        deltamask(StrawHitFlag::bkg);

    for (int i=0; i<nmc; i++) {
      McPart_t* mc    = _list_of_mc_particles.at(i);
      mc->fNHitsDelta = 0;
//-----------------------------------------------------------------------------
// loop over the hits of MC delta electron and calculate fraction of them 
// which have been tagged as the delta electron hits
//-----------------------------------------------------------------------------
      // int nh = mc->fListOfHits.size();
      // for (int ih=0; ih<nh; ih++) {
      // 	const StrawHitFlag* flag = mc->fListOfFlags.at(ih);

      // 	int flagged_as_delta = flag->hasAnyProperty(deltamask);

      // 	if (flagged_as_delta) mc->fNHitsDelta += 1;
      // }

      int pdg_id = mc->fSim->pdgId();
      
      if (pdg_id == 11  ) { 
	float mom = mc->Momentum();
	if (mom < 20)      {
//-----------------------------------------------------------------------------
// call this "a delta electron"
//-----------------------------------------------------------------------------
	  fNHitsDeltaTot += mc->NHits();
	}
      }
    }

    //    printf("TZClusterFinderAna::associateMcTruth: fNHitsDeltaTot = %5i\n",fNHitsDeltaTot);

    return 0;
  }

//-----------------------------------------------------------------------------
bool TZClusterFinderAna::findData(const art::Event& Evt) {
    _chcol    = 0;
    _mcdigis  = 0;

    auto shH  = Evt.getValidHandle<ComboHitCollection>(_chCollTag);
    _chcol    = shH.product();
    _nsh      = _chcol->size();

    auto mcdH = Evt.getValidHandle<StrawDigiMCCollection>(_mcdigiCollTag);
    _mcdigis  = mcdH.product();

    return (_chcol != 0) && (_nsh > 0) && (_mcdigis != 0) ;     
  }


//-----------------------------------------------------------------------------
  void TZClusterFinderAna::analyze(const art::Event& Event) {

    _eventNum = Event.event();
    if (_debugLevel) printf(">>> TZClusterFinderAna::produce  event number: %10i\n",_eventNum);  
//-----------------------------------------------------------------------------
// process event
//-----------------------------------------------------------------------------
    if (! findData(Event)) {
      throw cet::exception("RECO")<<"mu2e::TZClusterFinderAna_module::produce: data missing or incomplete"<< endl;
    }

    fNHitsDeltaTot  = 0;
    fNHitsDeltaReco = 0;

    initMcDiag();
    associateMcTruth();
//-----------------------------------------------------------------------------
// in the end of event processing fill diagnostic histograms
//-----------------------------------------------------------------------------
    fillHistograms();

    if (_debugLevel    > 0) debug();
  }

//-----------------------------------------------------------------------------
// debugLevel > 0: print seeds
//-----------------------------------------------------------------------------
  void TZClusterFinderAna::debug() {
//-----------------------------------------------------------------------------
// print MC electrons
//-----------------------------------------------------------------------------
    if (_printElectrons) {
      int nmc = _list_of_mc_particles.size();

      for (int i=0; i<nmc; i++) {
	McPart_t* mc = _list_of_mc_particles.at(i);
	const SimParticle* sim = mc->fSim;

	if ((sim->pdgId() == 11) && (mc->Time() > 550) && (mc->NHits()  >= _printElectronsMinNHits)) {

	  float fr = mc->fNHitsDelta/(mc->NHits()+1.e-3);

	  if (fr < _printElectronsMaxFReco) {

	    printf(" electron: sim.id = %10li mom = %10.3f time= %9.3f nhits = %3i nhits(delta): %3i first: %2i last: %2i",
		   sim->id().asInt(), mc->Momentum(), mc->Time(), 
		   mc->NHits(), 
		   mc->fNHitsDelta, 
		   mc->fFirstStation, mc->fLastStation);
	    printf(" fraction: %6.3f\n",fr);
	  }
	}
      }
    }
  }


// Part of the magic that makes this class a module.
DEFINE_ART_MODULE(TZClusterFinderAna)

}
