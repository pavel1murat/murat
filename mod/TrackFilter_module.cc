
// ======================================================================
//
// TrackFilter_module: allows filtering on the TrkQual
// when given a cut efficiency
//
// ======================================================================

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

#include "art/Framework/Principal/Handle.h"

#include <iostream>
#include <string>

#include "Offline/RecoDataProducts/inc/KalSeed.hh"

#include "TH1.h"

using namespace std;

namespace mu2e {

  class TrackFilter : public art::EDFilter {
  public:
    struct Config {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<art::InputTag> trkCollTag{Name("trkCollTag"), Comment("Track Collection Tag")};
      fhicl::Atom<double>        pMin      {Name("pMin"   )   , Comment("pMin, MeV/c"         )};
      fhicl::Atom<double>        nTrkMin   {Name("ntrkMin")   , Comment("min N(trk)"          )};
    };

    using Parameters = art::EDFilter::Table<Config>;
    explicit TrackFilter(const Parameters& conf);

    enum { kMaxEvtHistSets  = 100 };
    enum { kMaxTrkHistSets  = 100 };

    struct EvtHist_t {
      TH1F*   evt[2];
      TH1F*   ntrk;
      TH1F*   ntrk_good;
    };

    struct TrkHist_t {
      TH1F*   p;
    };

    struct Hist_t {
      EvtHist_t     evt[kMaxEvtHistSets];
      TrkHist_t     trk[kMaxTrkHistSets];
    } _hist;
//-----------------------------------------------------------------------------
// data
//-----------------------------------------------------------------------------
    struct TrackPar_t {
      const KalSegment* kseg;
      float             p;
    };

    struct Data_t {
      const art::Event* event;
      int               ntrk;
      int               ntrk_good;
      TrackPar_t        tp[100];               // assume less 100 tracks
      
    } _data;

  private:
    virtual void beginJob()                  override;
    virtual bool beginRun(art::Run&   run  ) override;
    virtual bool endRun  (art::Run&   run  ) override;
    virtual bool filter  (art::Event& event) override;

    void     bookEvtHistograms(art::TFileDirectory* Dir, EvtHist_t* Hist);
    void     bookTrkHistograms(art::TFileDirectory* Dir, TrkHist_t* Hist);
    void     bookHistograms();

    void     fillEvtHistograms(EvtHist_t* Hist);
    void     fillTrkHistograms(TrkHist_t* Hist, const KalSeed* Kseg, TrackPar_t* Tp);
    void     fillHistograms();

    art::InputTag                  _trkCollTag;
    int                            _nTrkMin;
    float                          _pMin;

    const mu2e::KalSeedCollection* _ksColl;
    bool                           _passed;
  };



//-----------------------------------------------------------------------------
  TrackFilter::TrackFilter(const Parameters& conf): 
    art::EDFilter{conf}
    , _trkCollTag(conf().trkCollTag())
    , _nTrkMin   (conf().nTrkMin()   )
    , _pMin      (conf().pMin()      )
    {
    }

//-----------------------------------------------------------------------------
  void TrackFilter::beginJob() {
    bookHistograms();
  }

//-----------------------------------------------------------------------------
  bool TrackFilter::beginRun(art::Run& run) {
    return true;
  }

//-----------------------------------------------------------------------------
  bool TrackFilter::endRun(art::Run& run) {
    return true;
  }

//-----------------------------------------------------------------------------
  void TrackFilter::bookEvtHistograms(art::TFileDirectory* Dir, EvtHist_t* Hist) {
    Hist->evt[0]    = Dir->make<TH1F>("evt_0", "Event number"      , 1000,   0,   1e4);
    Hist->evt[1]    = Dir->make<TH1F>("evt_1", "Event number"      , 1000,   0,   1e6);
    Hist->ntrk      = Dir->make<TH1F>("ntrk" , "N(tracks)"         ,   10,   0,    10);
    Hist->ntrk_good = Dir->make<TH1F>("ntrkg", "N(straw digi MCs )",   10,   0,    10);
  }

//-----------------------------------------------------------------------------
  void TrackFilter::bookTrkHistograms(art::TFileDirectory* Dir, TrkHist_t* Hist) {
    Hist->p     = Dir->make<TH1F>("p" , "P(track)"         ,   500,   0,   500);
  }

//-----------------------------------------------------------------------------
  void TrackFilter::bookHistograms() {
    art::ServiceHandle<art::TFileService> tfs;
//-----------------------------------------------------------------------------
// event-level histograms
//-----------------------------------------------------------------------------
    int book_evt_histograms[kMaxEvtHistSets];
    for (int i=0; i<kMaxEvtHistSets; i++)  book_evt_histograms[i] = 0;

    book_evt_histograms[0] = 1;     // all
    book_evt_histograms[1] = 1;     // passed

    for (int i=0; i<kMaxEvtHistSets; i++) {
      if (book_evt_histograms[i] == 0) continue;

      art::TFileDirectory dir = tfs->mkdir(Form("evt_%02i",i));
      bookEvtHistograms(&dir,&_hist.evt[i]);
    }

//-----------------------------------------------------------------------------
// track-level histograms
//-----------------------------------------------------------------------------
    int book_trk_histograms[kMaxTrkHistSets];
    for (int i=0; i<kMaxTrkHistSets; i++)  book_trk_histograms[i] = 0;

    book_trk_histograms[0] = 1;     // all

    for (int i=0; i<kMaxTrkHistSets; i++) {
      if (book_trk_histograms[i] == 0) continue;

      art::TFileDirectory dir = tfs->mkdir(Form("trk_%02i",i));
      bookTrkHistograms(&dir,&_hist.trk[i]);
    }


  }

//-----------------------------------------------------------------------------
  void TrackFilter::fillEvtHistograms(EvtHist_t* Hist) {
    Hist->evt[0]->Fill(_data.event->event());
    Hist->evt[1]->Fill(_data.event->event());
    Hist->ntrk->Fill(_data.ntrk);
    Hist->ntrk_good->Fill(_data.ntrk_good);
  }

//-----------------------------------------------------------------------------
  void TrackFilter::fillTrkHistograms(TrkHist_t* Hist, const KalSeed* Ks, TrackPar_t* Tp) {
    Hist->p->Fill(Tp->p);
  }

//-----------------------------------------------------------------------------
  void TrackFilter::fillHistograms() {

    fillEvtHistograms(&_hist.evt[0]);

    if (_passed) fillEvtHistograms(&_hist.evt[1]);

    for (int i=0; i<_data.ntrk; i++) {
      const KalSeed* ks = &_ksColl->at(i);
      fillTrkHistograms(&_hist.trk[0],ks,&_data.tp[i]);
    }
  }

//-----------------------------------------------------------------------------
  bool TrackFilter::filter(art::Event& event) {
    
    art::Handle<mu2e::KalSeedCollection> kscH;
    event.getByLabel(_trkCollTag,kscH);
    
    _ksColl         = nullptr;

    _data.event     = &event;
    _data.ntrk      = 0;
    _data.ntrk_good = 0;

    if (kscH.isValid())    { 
      _ksColl    = kscH.product();
      _data.ntrk = _ksColl->size();
    }

//-----------------------------------------------------------------------------
// look for the segment with minimal Z
//-----------------------------------------------------------------------------
    for (int i=0; i < _data.ntrk ; i++) {
      const mu2e::KalSeed* ks = &_ksColl->at(i);
      
      const mu2e::KalSegment *kseg(nullptr);

      double zmin(1.e6);

      for(auto const& ks : ks->segments() ) {
        double z = ks.position3().z();
        if (z < zmin) {
          kseg = &ks;
          zmin = z  ;  
        }
      }

      float p = kseg->mom();
      _data.tp[i].p    = p;
      _data.tp[i].kseg = kseg;

      if (p >= _pMin) _data.ntrk_good++;
    }

    _passed = (_data.ntrk_good >= _nTrkMin);
    fillHistograms();

    return _passed;
  }
}

using mu2e::TrackFilter;
DEFINE_ART_MODULE(TrackFilter)
