/////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////////////
#include <string>
#include <map>
#include <sstream>

// art includes.

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/OptionalTable.h"

#include "cetlib_except/exception.h"

#include "art_root_io/TFileService.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/StrawGasStep.hh"
#include "Offline/MCDataProducts/inc/CaloShowerStep.hh"
#include "Offline/MCDataProducts/inc/CrvStep.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"

#include "TH1F.h"

#include <map>

namespace mu2e {

  class DetStepAna : public art::EDAnalyzer {
    public:

    struct EvtHist_t {
      TH1F* event;
      TH1F* ns_trk;
      TH1F* ns_cal;
      TH1F* ns_crv;
    };
    
    struct TrkHist_t {
      TH1F* mom;
      TH1F* time;
      TH1F* edep;
    };
    
    struct CalHist_t {
      TH1F* mom;
      TH1F* time;
      TH1F* edep;
    };
    
    struct CrvHist_t {
      TH1F* mom;
      TH1F* time;
      TH1F* edep;
    };
    
    struct Hist_t {
      EvtHist_t evt;
      TrkHist_t trk;
      CalHist_t cal;
      CrvHist_t crv;
    } _hist;

    struct Config {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Sequence<art::InputTag> sgsCollTags { Name("sgsCollTags"), Comment("StrawGasStep   coll tags") };
      fhicl::Sequence<art::InputTag> cssCollTags { Name("cssCollTags"), Comment("CaloShowerStep coll tags") };
      fhicl::Sequence<art::InputTag> crvCollTags { Name("crvCollTags"), Comment("CrvStep coll tags") };
    };

    using    Parameters = art::EDAnalyzer::Table<Config>;
    explicit DetStepAna(const Parameters& conf);
    void     analyze(const art::Event& event) override;
    void     endJob() override;
    
    void      bookEvtHistograms  (EvtHist_t*   Hist, int HistSet, art::TFileDirectory* Dir);
    void      bookTrkHistograms  (TrkHist_t*   Hist, int HistSet, art::TFileDirectory* Dir);
    void      bookCalHistograms  (CalHist_t*   Hist, int HistSet, art::TFileDirectory* Dir);
    void      bookCrvHistograms  (CrvHist_t*   Hist, int HistSet, art::TFileDirectory* Dir);
    void      bookHistograms();

    private:
      std::vector<art::InputTag> _sgsCollTags, _cssCollTags, _crvCollTags;
  };


//================================================================
  DetStepAna::DetStepAna(const Parameters& conf) : art::EDAnalyzer{conf}
    {
      for(const auto& tag : conf().sgsCollTags()) { _sgsCollTags.emplace_back(tag); consumes<StrawGasStepCollection>(tag); }
      for(const auto& tag : conf().cssCollTags()) { _cssCollTags.emplace_back(tag); consumes<CaloShowerStepCollection>(tag); }
      for(const auto& tag : conf().crvCollTags()) { _crvCollTags.emplace_back(tag); consumes<CrvStepCollection>(tag); }
      
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory tfdir = tfs->mkdir("DetStepAna");

      TH1::AddDirectory(0);
      
      art::TFileDirectory d0 = tfdir.mkdir("evt");
      bookEvtHistograms(&_hist.evt,0,&d0);
      
      art::TFileDirectory d1 = tfdir.mkdir("trk");
      bookTrkHistograms(&_hist.trk,0,&d1);
      
      art::TFileDirectory d2 = tfdir.mkdir("cal");
      bookCalHistograms(&_hist.cal,0,&d2);
      
      art::TFileDirectory d3 = tfdir.mkdir("crv");
      bookCrvHistograms(&_hist.crv,0,&d3);
      
    }

//-----------------------------------------------------------------------------
  void DetStepAna::bookEvtHistograms(EvtHist_t* Hist, int HistSet, art::TFileDirectory* Dir) {

    Hist->event   = Dir->make<TH1F>("event" , "Event Number",  100, 0., 100000.);
    Hist->ns_trk  = Dir->make<TH1F>("ns_trk", "N(trk steps)", 1000, -0.5,   999.5);
    Hist->ns_cal  = Dir->make<TH1F>("ns_cal", "N(cal steps)", 1000, -0.5,   999.5);
    Hist->ns_crv  = Dir->make<TH1F>("ns_crv", "N(crv steps)", 1000, -0.5,   999.5);
  }

//-----------------------------------------------------------------------------
  void DetStepAna::bookTrkHistograms(TrkHist_t* Hist, int HistSet, art::TFileDirectory* Dir) {

    Hist->mom   = Dir->make<TH1F>("mom" , "mom" ,  500, 0.,    500.);
    Hist->time  = Dir->make<TH1F>("time", "time", 2000, 0., 100000.);
    Hist->edep  = Dir->make<TH1F>("edep", "edep", 1000, 0.,     0.1);
  }


//-----------------------------------------------------------------------------
  void DetStepAna::bookCalHistograms(CalHist_t* Hist, int HistSet, art::TFileDirectory* Dir) {

    Hist->mom   = Dir->make<TH1F>("mom" , "mom" ,  500, 0.,    500.);
    Hist->time  = Dir->make<TH1F>("time", "time", 2000, 0., 100000.);
    Hist->edep  = Dir->make<TH1F>("edep", "edep",  250, 0.,     50.);
  }


//-----------------------------------------------------------------------------
  void DetStepAna::bookCrvHistograms(CrvHist_t* Hist, int HistSet, art::TFileDirectory* Dir) {

    Hist->mom   = Dir->make<TH1F>("mom" , "mom" ,  500, 0.,    500.);
    Hist->time  = Dir->make<TH1F>("time", "time", 2000, 0., 100000.);
    Hist->edep  = Dir->make<TH1F>("edep", "edep", 1000, 0.,   1000.);
  }

  
//-----------------------------------------------------------------------------  
  void DetStepAna::analyze(const art::Event& event) {
    
    double mbtime = GlobalConstantsHandle<PhysicsParams>()->getNominalDRPeriod();
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
    int ns_trk = 0;
    for(const auto& tag : _sgsCollTags) {
      auto sgscolH = event.getValidHandle<StrawGasStepCollection>(tag);
      for(const auto& sgs : *sgscolH ) {
        double mom  = sgs.momentum().R();
        double edep = sgs.ionizingEdep();
        double time = sgs.time();

        _hist.trk.mom->Fill(mom);
        _hist.trk.edep->Fill(edep);
        _hist.trk.time->Fill(time);
        ns_trk++;
      }
    }
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
    int ns_cal = 0;
    for(const auto& tag : _cssCollTags) {
      auto csscolH = event.getValidHandle<CaloShowerStepCollection>(tag);
      for(const auto& css : *csscolH ) {
        double mom  = css.momentumIn();
        double edep = css.energyDepBirks();
        double time = css.time();
        _hist.cal.mom->Fill(mom);
        _hist.cal.edep->Fill(edep);
        _hist.cal.time->Fill(time);
        ns_cal++;
      }
    }
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
    int ns_crv = 0;
    for(const auto& tag : _crvCollTags) {
      auto crvscolH = event.getValidHandle<CrvStepCollection>(tag);
      for(const auto& crvs : *crvscolH ) {
        double mom = crvs.startMom().R();
        double edep = crvs.visibleEDep();
        double time = fmod(crvs.startTime(),mbtime);
        _hist.cal.mom->Fill(mom);
        _hist.cal.edep->Fill(edep);
        _hist.cal.time->Fill(time);
        ns_crv++;
      }
    }
//-----------------------------------------------------------------------------
// fill event histograms
//-----------------------------------------------------------------------------
    int evt    = event.event();
    
    _hist.evt.event->Fill(evt);
    if (ns_trk > 0) _hist.evt.ns_trk->Fill(ns_trk);
    if (ns_cal > 0) _hist.evt.ns_cal->Fill(ns_cal);
    if (ns_crv > 0) _hist.evt.ns_crv->Fill(ns_crv);
  }

  void DetStepAna::endJob() {
    // mf::LogInfo("Summary")
    //   <<"DetStepAna_module: passed "
    //   <<nPassed_<<" / "<<nEvt_<<" events\n";
  }

}

DEFINE_ART_MODULE(mu2e::DetStepAna)
