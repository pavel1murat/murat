///////////////////////////////////////////////////////////////////////////////
// $Id: GenpHist_module.cc,v 1.11 2014/10/02 17:15:09 murat Exp $
// $Author: murat $
// $Date: 2014/10/02 17:15:09 $
//
// .fcl file to use: murat/test/genp_hist.fcl
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

#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/StrawGasStep.hh"
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"

#include "Offline/RecoDataProducts/inc/ComboHit.hh"

#include "Stntuple/mod/StntupleModule.hh"

#include "TH1F.h"

using namespace std;
using CLHEP::Hep3Vector;

namespace mu2e {

  class GenpHist : public StntupleModule {
  public:
    struct Config {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag>   genpCollTag {Name("genpCollTag"), Comment("GenParticle  coll tag") };
      fhicl::Atom<art::InputTag>   sgsCollTag  {Name("sgsCollTag" ), Comment("StrawGasStep coll tag") };
      fhicl::Atom<art::InputTag>   sdmcCollTag {Name("sdmcCollTag"), Comment("StrawDigiMC  coll tag") };
      fhicl::Atom<art::InputTag>   sschCollTag {Name("sschCollTag"), Comment("SSCH         coll tag") };
    };
//-----------------------------------------------------------------------------
// Module labels 
//-----------------------------------------------------------------------------
    art::InputTag      _genpCollTag;
    art::InputTag      _sgsCollTag;
    art::InputTag      _sdmcCollTag;
    art::InputTag      _sschCollTag;

    const GenParticleCollection*  _genpColl;
    const StrawGasStepCollection* _sgsColl;
    const StrawDigiMCCollection*  _sdmcColl;
    const ComboHitCollection*     _sschColl;

    struct Data_t {
      const art::Event* event;
      int               nsgs;              // number of straw gas steps
      int               nsdmc;              // number of straw gas steps
      int               nssch;             // N single-straw combo hits
    } _data;

//-----------------------------------------------------------------------------
// histogramming
//-----------------------------------------------------------------------------
    enum { kMaxEvtHistSets  = 100 };
    enum { kMaxGenpHistSets = 100 };
    enum { kMaxSgsHistSets  = 100 };
    enum { kMaxSdmcHistSets = 100 };
    enum { kMaxSschHistSets = 100 };

    struct EvtHist_t {
      TH1F*   evt[2];
      TH1F*   nsgs;
      TH1F*   nsdmc;
      TH1F*   nssch;
    };

    struct GenpHist_t {
      TH1F*   energy [2];
      TH1F*   mom    [2];//to be filled
      TH1F*   pdgId;
      TH1F*   costh;
      TH1F*   phi;
    };

    struct SgsHist_t {
      TH1F*   strawId;
      TH1F*   edep;
      TH1F*   time;
    };

    struct SdmcHist_t {
      TH1F*   strawId;
      TH1F*   esum;
      TH1F*   ctime;  // "early end" cluster time
      TH1F*   wtime;  // "early end" wire time
    };

    struct SschHist_t {
      TH1F*   strawId;
      TH1F*   edep;
      TH1F*   time;
    };

    struct Hist_t {
      EvtHist_t     evt [kMaxEvtHistSets ];
      GenpHist_t    genp[kMaxGenpHistSets];
      SgsHist_t     sgs [kMaxSgsHistSets ];
      SdmcHist_t    sdmc[kMaxSdmcHistSets];
      SschHist_t    ssch[kMaxSschHistSets];
    } _hist;

  public:
    explicit GenpHist(fhicl::ParameterSet const& pset);
    virtual ~GenpHist();

    void     bookEvtHistograms (art::TFileDirectory* Dir, EvtHist_t*  Hist);
    void     bookGenpHistograms(art::TFileDirectory* Dir, GenpHist_t* Hist);
    void     bookSgsHistograms (art::TFileDirectory* Dir, SgsHist_t*  Hist);
    void     bookSdmcHistograms(art::TFileDirectory* Dir, SdmcHist_t* Hist);
    void     bookSschHistograms(art::TFileDirectory* Dir, SschHist_t* Hist);
    void     bookHistograms();

    void     fillEvtHistograms (EvtHist_t*  Hist);
    void     fillGenpHistograms(GenpHist_t* Hist, const GenParticle*  Genp);
    void     fillSgsHistograms (SgsHist_t*  Hist, const StrawGasStep* Sgs );
    void     fillSdmcHistograms(SdmcHist_t* Hist, const StrawDigiMC*  Sdmc);
    void     fillSschHistograms(SschHist_t* Hist, const ComboHit*     Ssch);
    void     fillHistograms();

    void     getData(const art::Event& rEvent);
//-----------------------------------------------------------------------------
// overloaded virtual methods of the base class
//-----------------------------------------------------------------------------
    virtual void     beginJob()                         override;
    virtual void     beginRun(const art::Run& )         override;
    virtual void     endJob  ()                         override;
    virtual void     analyze (const art::Event& rEvent) override;
  };


//-----------------------------------------------------------------------------
  GenpHist::GenpHist(fhicl::ParameterSet const& pset): 
    StntupleModule (pset,"GenpHist")
    ,_genpCollTag (pset.get<art::InputTag>("genpCollTag"))
    ,_sgsCollTag  (pset.get<art::InputTag>("sgsCollTag" ))
    ,_sdmcCollTag (pset.get<art::InputTag>("sdmcCollTag"))
    ,_sschCollTag (pset.get<art::InputTag>("sschCollTag"))
  {
  }

//-----------------------------------------------------------------------------
  GenpHist::~GenpHist() { 
  }

//-----------------------------------------------------------------------------
  void GenpHist::endJob() {
    art::ServiceHandle<art::TFileService> tfs;
  }


//-----------------------------------------------------------------------------
  void GenpHist::bookEvtHistograms(art::TFileDirectory* Dir, EvtHist_t* Hist) {
    Hist->evt[0] = Dir->make<TH1F>("evt_0", "Event number"      , 1000,   0,   1e4);
    Hist->evt[1] = Dir->make<TH1F>("evt_1", "Event number"      , 1000,   0,   1e6);
    Hist->nsgs   = Dir->make<TH1F>("nsgs" , "N(straw gas steps)", 1000,   0,   1000);
    Hist->nsdmc  = Dir->make<TH1F>("nsdmc", "N(straw digi MCs )", 1000,   0,   1000);
    Hist->nssch  = Dir->make<TH1F>("nssch", "N(ss combo hits)"  , 1000,   0,   1000);
  }

//-----------------------------------------------------------------------------
  void GenpHist::bookGenpHistograms(art::TFileDirectory* Dir, GenpHist_t* Hist) {
    Hist->energy[0] = Dir->make<TH1F>("e_0"  , "Energy[0]"   , 2400,   0.0,   120);
    Hist->energy[1] = Dir->make<TH1F>("e_1"  , "Energy[1]"   , 2400,   0.0,  1200);
    Hist->mom   [0] = Dir->make<TH1F>("p_0"  , "Momentum"    ,  500,   0.0,  1000);
    Hist->mom   [1] = Dir->make<TH1F>("p_1"  , "Momentum"    ,  500,   0.0,  1000); //want to know z momentum as well
    Hist->pdgId     = Dir->make<TH1F>("pdgId", "PDG ID"      ,  500,  -250,   250);
    Hist->costh     = Dir->make<TH1F>("costh", "cos(th)"     ,  200,  -1  ,   1);
    Hist->phi       = Dir->make<TH1F>("phi"  , "phi"         ,  315,  -3.15,   3.15);
  }

//-----------------------------------------------------------------------------
  void GenpHist::bookSgsHistograms(art::TFileDirectory* Dir, SgsHist_t* Hist) {
    Hist->strawId = Dir->make<TH1F>("str_id", "straw ID"     , 2500,   0,   25000);
    Hist->edep    = Dir->make<TH1F>("edep"  , "E(dep)"       , 500 ,   0,   0.1);
    Hist->time    = Dir->make<TH1F>("time"  , "time"         , 200 ,   0,   2000);
  }

//-----------------------------------------------------------------------------
  void GenpHist::bookSdmcHistograms(art::TFileDirectory* Dir, SdmcHist_t* Hist) {
    Hist->strawId = Dir->make<TH1F>("str_id", "straw ID"              , 2500,   0,   25000);
    Hist->esum    = Dir->make<TH1F>("esum"  , "E(sum)"                ,  500,   0,   0.1);
    Hist->ctime   = Dir->make<TH1F>("ctime" , "early end cluster time",  200,   0,   2000);
    Hist->wtime   = Dir->make<TH1F>("wtime" , "early end wire time"   ,  200,   0,   2000);
  }

//-----------------------------------------------------------------------------
  void GenpHist::bookSschHistograms(art::TFileDirectory* Dir, SschHist_t* Hist) {
    Hist->strawId = Dir->make<TH1F>("str_id", "straw ID"     , 2500,   0,   25000);
    Hist->edep    = Dir->make<TH1F>("edep"  , "E(dep)"       , 500 ,   0,   0.1);
    Hist->time    = Dir->make<TH1F>("time"  , "time"         , 200 ,   0,   2000);
  }

//-----------------------------------------------------------------------------
  void GenpHist::bookHistograms() {

    art::ServiceHandle<art::TFileService> tfs;

//-----------------------------------------------------------------------------
// single-straw combo hits
//-----------------------------------------------------------------------------
    int book_ssch_histograms[kMaxSschHistSets];
    for (int i=0; i<kMaxSschHistSets; i++)  book_ssch_histograms[i] = 0;

    book_ssch_histograms[0] = 1;     // all

    for (int i=0; i<kMaxSschHistSets; i++) {
      if (book_ssch_histograms[i] == 0) continue;

      art::TFileDirectory dir = tfs->mkdir(Form("ssch_%02i",i));
      bookSschHistograms(&dir,&_hist.ssch[i]);
    }
//-----------------------------------------------------------------------------
// event-level histograms
//-----------------------------------------------------------------------------
    int book_evt_histograms[kMaxEvtHistSets];
    for (int i=0; i<kMaxEvtHistSets; i++)  book_evt_histograms[i] = 0;

    book_evt_histograms[0] = 1;     // all

    for (int i=0; i<kMaxEvtHistSets; i++) {
      if (book_evt_histograms[i] == 0) continue;

      art::TFileDirectory dir = tfs->mkdir(Form("evt_%02i",i));
      bookEvtHistograms(&dir,&_hist.evt[i]);
    }
//-----------------------------------------------------------------------------
// GenParticle histograms
//-----------------------------------------------------------------------------
    int book_genp_histograms[kMaxGenpHistSets];
    for (int i=0; i<kMaxGenpHistSets; i++)  book_genp_histograms[i] = 0;

    book_genp_histograms[0] = 1;     // all
    book_genp_histograms[1] = 1;     // electrons

    for (int i=0; i<kMaxGenpHistSets; i++) {
      if (book_genp_histograms[i] == 0) continue;

      art::TFileDirectory dir = tfs->mkdir(Form("genp_%02i",i));
      bookGenpHistograms(&dir,&_hist.genp[i]);
    }
//-----------------------------------------------------------------------------
// straw gas step-level histograms
//-----------------------------------------------------------------------------
    int book_sgs_histograms[kMaxSgsHistSets];
    for (int i=0; i<kMaxSgsHistSets; i++)  book_sgs_histograms[i] = 0;

    book_sgs_histograms[0] = 1;     // all

    for (int i=0; i<kMaxSgsHistSets; i++) {
      if (book_sgs_histograms[i] == 0) continue;

      art::TFileDirectory dir = tfs->mkdir(Form("sgs_%02i",i));
      bookSgsHistograms(&dir,&_hist.sgs[i]);
    }
//-----------------------------------------------------------------------------
// StrawDigiMC-level histograms
//-----------------------------------------------------------------------------
    int book_sdmc_histograms[kMaxSdmcHistSets];
    for (int i=0; i<kMaxSdmcHistSets; i++)  book_sdmc_histograms[i] = 0;

    book_sdmc_histograms[0] = 1;     // all

    for (int i=0; i<kMaxSdmcHistSets; i++) {
      if (book_sdmc_histograms[i] == 0) continue;

      art::TFileDirectory dir = tfs->mkdir(Form("sdmc_%02i",i));
      bookSdmcHistograms(&dir,&_hist.sdmc[i]);
    }
  }

//-----------------------------------------------------------------------------
  void GenpHist::beginJob() {
    bookHistograms();
  }

//-----------------------------------------------------------------------------
  void GenpHist::beginRun(const art::Run& ) {
  }


//-----------------------------------------------------------------------------
  void GenpHist::fillSschHistograms(SschHist_t* Hist, const ComboHit* Ssch) {
    Hist->strawId->Fill(Ssch->strawId().asUint16());
    Hist->edep->Fill(Ssch->energyDep());
    Hist->time->Fill(Ssch->endTime(Ssch->earlyEnd()));
  }

//-----------------------------------------------------------------------------
  void GenpHist::fillEvtHistograms(EvtHist_t* Hist) {

    Hist->evt[0]->Fill(_data.event->event());
    Hist->evt[1]->Fill(_data.event->event());
    Hist->nsgs->Fill(_data.nsgs);
    Hist->nsdmc->Fill(_data.nsdmc);
  }

//-----------------------------------------------------------------------------
  void GenpHist::fillGenpHistograms(GenpHist_t* Hist, const GenParticle* Genp) {

    Hist->pdgId->Fill(Genp->pdgId());
    
    double e = Genp->momentum().e();
    double p = Genp->momentum().vect().mag();
    
    Hist->energy[0]->Fill(e);
    Hist->energy[1]->Fill(e);
    Hist->mom   [0]->Fill(p);
    Hist->mom   [1]->Fill(p);

    
    double costh(-100.), phi(-100.);

    if (p > 0) { 
      costh = Genp->momentum().z()/p;
      phi   = atan2(Genp->momentum().y(),Genp->momentum().x());
    }

    Hist->costh->Fill(costh);
    Hist->phi->Fill(phi);
  }


//-----------------------------------------------------------------------------
  void GenpHist::fillSgsHistograms(SgsHist_t* Hist, const StrawGasStep* Sgs) {
    Hist->strawId->Fill(Sgs->strawId().asUint16());
    Hist->edep->Fill(Sgs->ionizingEdep());
    Hist->time->Fill(Sgs->time());
  }

//-----------------------------------------------------------------------------
  void GenpHist::fillSdmcHistograms(SdmcHist_t* Hist, const StrawDigiMC* Sdmc) {
    Hist->strawId->Fill(Sdmc->strawId().asUint16());
    Hist->esum->Fill(Sdmc->energySum());
    StrawEnd early_end = Sdmc->earlyEnd();
    Hist->ctime->Fill(Sdmc->clusterTime(early_end));
    Hist->wtime->Fill(Sdmc->wireEndTime(early_end));
  }

//-----------------------------------------------------------------------------
  void GenpHist::fillHistograms() {

    fillEvtHistograms(&_hist.evt[0]);

    if (_genpColl) {
      for (const GenParticle& genp : *_genpColl) {
        fillGenpHistograms(&_hist.genp[0],&genp);
        
        if (abs(genp.pdgId()) == 11) fillGenpHistograms(&_hist.genp[1],&genp);
      }
    }

    if (_sgsColl) {
      for (const StrawGasStep& sgs : *_sgsColl) {
        fillSgsHistograms(&_hist.sgs[0],&sgs);
      }
    }

    if (_sdmcColl) {
      for (const StrawDigiMC& sdmc : *_sdmcColl) {
        fillSdmcHistograms(&_hist.sdmc[0],&sdmc);
      }
    }

    if (_sschColl) {
      for (const ComboHit& ssch : *_sschColl) {
        fillSschHistograms(&_hist.ssch[0],&ssch);
      }
    }
  }

//-----------------------------------------------------------------------------
  void GenpHist::getData(const art::Event& rEvent) {
    const char* oname = "GenpHist::getData";
    bool ok;

    _data.event = &rEvent;
    _data.nsgs  = 0;

    art::Handle<mu2e::ComboHitCollection> sschcH;
    ok = rEvent.getByLabel(_sschCollTag,sschcH); 
    if (ok) {
      _sschColl   = sschcH.product();
      _data.nssch = _sschColl->size();
    }
    else {
      _sschColl  = nullptr;
      mf::LogWarning(oname) << " WARNING in " << oname << ":" << __LINE__ 
                            << ": ComboHitCollection:" 
                            << _sschCollTag.encode().data() << " NOT FOUND";
    }

    art::Handle<mu2e::StrawGasStepCollection> sgscH;
    ok = rEvent.getByLabel(_sgsCollTag,sgscH); 
    if (ok) {
      _sgsColl   = sgscH.product();
      _data.nsgs = _sgsColl->size();
    }
    else {
      _sgsColl  = nullptr;
      mf::LogWarning(oname) << " WARNING in " << oname << ":" << __LINE__ 
                            << ": StrawGasStepCollection:" 
                            << _sgsCollTag.encode().data() << " NOT FOUND";
    }

    art::Handle<mu2e::StrawDigiMCCollection> sdmccH;
    ok = rEvent.getByLabel(_sdmcCollTag,sdmccH); 
    if (ok) {
      _sdmcColl   = sdmccH.product();
      _data.nsdmc = _sdmcColl->size();
    }
    else {
      _sdmcColl  = nullptr;
      mf::LogWarning(oname) << " WARNING in " << oname << ":" << __LINE__ 
                            << ": StrawDigiMCCollection:" 
                            << _sdmcCollTag.encode().data() << " NOT FOUND";
    }

    art::Handle<mu2e::GenParticleCollection> genpcH;
    ok = rEvent.getByLabel(_genpCollTag,genpcH); 
    if (ok) _genpColl = genpcH.product();
    else {
      _genpColl = nullptr;
      mf::LogWarning(oname) << " WARNING in " << oname << ":" << __LINE__ 
                            << ": GenParticleCollection:" 
                            << _genpCollTag.encode().data() << " NOT FOUND";
    }
  }

//-----------------------------------------------------------------------------
  void GenpHist::analyze(const art::Event& rEvent) {

    getData(rEvent);

    fillHistograms();
  }
}

using mu2e::GenpHist;
DEFINE_ART_MODULE(GenpHist)
