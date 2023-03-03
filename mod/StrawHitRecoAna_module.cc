//////////////////////////////////////////////////////////////////////////////
// framework
//
// parameter defaults: CalPatRec/fcl/prolog.fcl
// this module doesn't do reconstruction
// on input, it takes a list  of StrawHitFlags flags and evaluates performance
// of the delta electron tagging
//////////////////////////////////////////////////////////////////////////////
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "art/Framework/Core/EDAnalyzer.h"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "art_root_io/TFileService.h"
// conditions
#include "Offline/ConditionsService/inc/ConditionsHandle.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"

#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom3.h"

#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitPosition.hh"
#include "Offline/RecoDataProducts/inc/StereoHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"

#include <algorithm>
#include <cmath>
#include "CLHEP/Vector/ThreeVector.h"
#include "Offline/MCDataProducts/inc/StrawGasStep.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"

#include "Offline/CalPatRec/inc/HlPrint.hh"

using namespace std;
using CLHEP::Hep3Vector;

namespace mu2e {

  class StrawHitRecoAna : public art::EDAnalyzer {

    struct Config {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag>    shCollTag              {Name("shCollTag"             ), Comment("SComboHit collection tag"         ) };
      fhicl::Atom<art::InputTag>    chCollTag              {Name("chCollTag"             ), Comment("ComboHit collection tag"          ) };
      fhicl::Atom<art::InputTag>    sschCollTag            {Name("sschCollTag"           ), Comment("SS ComboHit collection tag"       ) };
      fhicl::Atom<art::InputTag>    chfCollTag             {Name("chfCollTag"            ), Comment("ComboHit flag collection tag"     ) };
      fhicl::Atom<art::InputTag>    shfCollTag             {Name("shfCollTag"            ), Comment("SSC Hit  flag collection tag"     ) };
      fhicl::Atom<art::InputTag>    sdmcCollTag            {Name("sdmcCollTag"           ), Comment("StrawDigiMC collection Name"      ) };

      fhicl::Atom<float>            scale2                 {Name("scale2"                ), Comment("scale for 2-hitters "             ) };
      fhicl::Atom<float>            scale3                 {Name("scale3"                ), Comment("scale for 3-hitters "             ) };

      fhicl::Atom<int>              debugLevel             {Name("debugLevel"            ), Comment("debug level"                      ) };
      fhicl::Atom<int>              diagLevel              {Name("diagLevel"             ), Comment("diag level"                       ) };
      fhicl::Atom<int>              printComboHits         {Name("printComboHits"        ), Comment("if 1, print combo hits"           ) };
      fhicl::Atom<int>              printSingleComboHits   {Name("printSingleComboHits"  ), Comment("if 1, print single straw ComboH"  ) };
    };

  public:
    enum { kNStations      = 20 };
    enum { kNFaces         =  4 };
    enum { kNPanelsPerFace =  3 };

    enum {
      kNEventHistSets    =  10,
      kNComboHitHistSets = 100,
    };

    struct ComboHitPar_t {
      int    nsh;
      float  werr;
      float  sigw;
      float  sigw_over_werr;
      float  sigws;
      float  sigws_over_werr;

      float  ss_wres[10];
    };

  protected:

    struct ComboHitHist_t {
      TH1F* fNsh;
      TH1F* fWres;
      TH1F* fSigw;
      TH1F* fSigwOverWres;
      TH1F* fSigws;
      TH1F* fSigwsOverWres;
      TH1F* fSSWres;
    };

    struct EventHist_t {
      TH1F*  fEventNumber;
      TH1F*  fNsh;
      TH1F*  fNch;
      TH1F*  fNssch;
    };

    struct Hist_t {
      EventHist_t*    fEvent   [kNEventHistSets   ];
      ComboHitHist_t* fComboHit[kNComboHitHistSets];
    };

    Hist_t  _hist;

//-----------------------------------------------------------------------------
// talk-to parameters
//-----------------------------------------------------------------------------
    art::InputTag                  _shCollTag;              // straw hits        by "makeSH"
    art::InputTag                  _chCollTag;
    art::InputTag                  _sschCollTag;            // single straw combohits by "makeSH"
    art::InputTag                  _shfCollTag;             // 1-straw  ComboHit flags
    art::InputTag                  _chfCollTag;             // combined ComboHit flags
    art::InputTag                  _sdmcCollTag;
    float                          _scale2;
    float                          _scale3;
    int                            _debugLevel;
    int                            _diagLevel;
    int                            _printComboHits;
    int                            _printSingleComboHits;
//-----------------------------------------------------------------------------
// cache of event or geometry objects
//-----------------------------------------------------------------------------
    const ComboHitCollection*      _chColl;
    const ComboHitCollection*      _sschColl; // one straw hit per combo hit
    const StrawHitCollection*      _shColl;
    const StrawHitFlagCollection*  _chfColl;
    const StrawDigiMCCollection*   _sdmcColl;

    const Tracker*                 _tracker;

    int                            _eventNum;
    int                            _nSingleSH;  // true number of straw hits
    int                            _nComboHits;
    int                            _nStrawHits; // not the total number of straw hits, but the
                                                // number of straw hits from combohits

    int                            _n2;
    int                            _n3;

    HlPrint*                       _hlp;

    TRandom3                       _trn;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  public:
    explicit  StrawHitRecoAna(const art::EDAnalyzer::Table<Config>& config);
    virtual   ~StrawHitRecoAna();

    int       associateMcTruth();

    void      bookEventHistograms   (EventHist_t*    Hist, int HistSet, art::TFileDirectory* Dir);
    void      bookComboHitHistograms(ComboHitHist_t* Hist, int HistSet, art::TFileDirectory* Dir);
    void      bookHistograms();

    void      debug();

    void      fillEventHistograms   (EventHist_t*    Hist);
    void      fillComboHitHistograms(ComboHitHist_t* Hist, const ComboHit* Hit, ComboHitPar_t* Chp);
    void      fillHistograms();

    bool      findData     (const art::Event&  Evt);
//-----------------------------------------------------------------------------
// overloaded methods of the base class
//-----------------------------------------------------------------------------
    virtual void beginJob();
    virtual void beginRun(const art::Run&   r);
    virtual void endRun  (const art::Run&   r);
    virtual void analyze (const art::Event& e);
  };

//-----------------------------------------------------------------------------
  StrawHitRecoAna::StrawHitRecoAna(const art::EDAnalyzer::Table<Config>& config):
    art::EDAnalyzer(config),
    _shCollTag             (config().shCollTag  ()      ),
    _chCollTag             (config().chCollTag  ()      ),
    _sschCollTag           (config().sschCollTag()      ),
    _shfCollTag            (config().shfCollTag ()      ),
    _chfCollTag            (config().chfCollTag ()      ),
    _sdmcCollTag           (config().sdmcCollTag()      ),
    _scale2                (config().scale2     ()      ),
    _scale3                (config().scale3     ()      ),
    _debugLevel            (config().debugLevel ()      ),
    _diagLevel             (config().diagLevel  ()      ),
    _printComboHits        (config().printComboHits        ()),
    _printSingleComboHits  (config().printSingleComboHits  ())
  {
    _hlp = HlPrint::Instance();

  }

  StrawHitRecoAna::~StrawHitRecoAna() {
  }

//-----------------------------------------------------------------------------
  void StrawHitRecoAna::bookEventHistograms(EventHist_t* Hist, int HistSet, art::TFileDirectory* Dir) {

    Hist->fEventNumber   = Dir->make<TH1F>(Form("event_%02i", HistSet), "Event Number", 100, 0., 100000.);
    Hist->fNch           = Dir->make<TH1F>(Form("nch_%02i"  , HistSet), "N(combo hits)"  ,1000, 0., 10000.);
    Hist->fNsh           = Dir->make<TH1F>(Form("nsh_%02i"  , HistSet), "N(straw hits)"  ,1000, 0., 10000.);
    Hist->fNssch         = Dir->make<TH1F>(Form("nssch_%02i", HistSet), "N(1-straw hits)",1000, 0., 10000.);
  }

//-----------------------------------------------------------------------------
  void StrawHitRecoAna::bookComboHitHistograms(ComboHitHist_t* Hist, int HistSet, art::TFileDirectory* Dir) {

    Hist->fNsh           = Dir->make<TH1F>("nsh"            , "N(straw hits)",   10, 0.,  10.);
    Hist->fSSWres        = Dir->make<TH1F>("ss_wres"        , "ss Wres"      ,  200, 0., 200.);
    Hist->fWres          = Dir->make<TH1F>("wres"           , "Wres"         ,  200, 0., 200.);
    Hist->fSigw          = Dir->make<TH1F>("sigw"           , "sigw"         ,  200, 0., 200.);
    Hist->fSigwOverWres  = Dir->make<TH1F>("sigw_over_wres" , "sigw/wres"    ,  200, 0.,  10.);
    Hist->fSigws         = Dir->make<TH1F>("sigws"          , "sigws"         ,  200, 0., 200.);
    Hist->fSigwsOverWres = Dir->make<TH1F>("sigws_over_wres", "sigw/wres"    ,  200, 0.,  10.);
  }

//-----------------------------------------------------------------------------
  void StrawHitRecoAna::bookHistograms() {
//-----------------------------------------------------------------------------
// book event histograms
//-----------------------------------------------------------------------------
    art::ServiceHandle<art::TFileService> tfs;
    char   folder_name[200];

    TH1::AddDirectory(0);

    int book_event_histset[kNEventHistSets];
    for (int i=0; i<kNEventHistSets; i++) book_event_histset[i] = 0;

    book_event_histset[ 0] = 1;                // all events

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
    int book_combo_hit_histset[kNComboHitHistSets];
    for (int i=0; i<kNComboHitHistSets; i++) book_combo_hit_histset[i] = 0;

    book_combo_hit_histset[  0] = 1;                // all
    book_combo_hit_histset[  1] = 1;                // nsh>1
    book_combo_hit_histset[  2] = 1;                // nsh= 2
    book_combo_hit_histset[  3] = 1;                // nsh= 3

    book_combo_hit_histset[ 12] = 1;                // MC nsh= 2
    book_combo_hit_histset[ 13] = 1;                // MC nsh= 3

    for (int i=0; i<kNComboHitHistSets; i++) {
      if (book_combo_hit_histset[i] != 0) {
        sprintf(folder_name,"sh_%i",i);
        art::TFileDirectory tfdir = tfs->mkdir(folder_name);

        _hist.fComboHit[i] = new ComboHitHist_t;
        bookComboHitHistograms(_hist.fComboHit[i],i,&tfdir);
      }
    }
  }

//-----------------------------------------------------------------------------
  void StrawHitRecoAna::beginJob() {
    bookHistograms();
  }

//----------------------------------------------------------------------------------------------------
  void StrawHitRecoAna::beginRun(const art::Run& R) {

    mu2e::GeomHandle<mu2e::Tracker> ttHandle;
    _tracker = ttHandle.get();
  }

//----------------------------------------------------------------------------------------------------
  void StrawHitRecoAna::endRun(const art::Run& R) {

    ComboHitPar_t chp;

    TRandom3 trn;

    TH1F* hrn = _hist.fComboHit[0]->fSSWres;

    float sig[10];
//-----------------------------------------------------------------------------
// total number of two-sh combo hits
//-----------------------------------------------------------------------------
    int nh = 2;

    _n2 = _hist.fComboHit[2]->fNsh->GetEntries();

    for (int i=0; i<_n2; i++) {
      float sy(0), sy2(0), sw(0);
      for (int k=0; k<nh; k++) {
                                        // resolution
        sig[k] = hrn->GetRandom();

        float y = trn.Gaus(0,sig[k]);
        
        sy  += y;
        sy2 += y*y;
        sw  += 1./(sig[k]*sig[k]);
      }
    
      double w2m = sy2/nh;
      double wm  = sy /nh;

      chp.werr           = sqrt(1/sw);
      chp.sigw           = sqrt(w2m-wm*wm)*sqrt(nh/(nh-1.));
      chp.sigw_over_werr = chp.sigw/chp.werr;
      
      chp.sigws           = chp.sigw*_scale2;
      chp.sigws_over_werr = chp.sigws/chp.werr;
      
      fillComboHitHistograms(_hist.fComboHit[12],nullptr, &chp);  // all
    } 
//-----------------------------------------------------------------------------
// the same for nsh=3 hits
//-----------------------------------------------------------------------------
    nh = 3;

    _n3 = _hist.fComboHit[3]->fNsh->GetEntries();

    for (int i=0; i<_n3; i++) {
      float sy(0), sy2(0), sw(0);
      for (int k=0; k<nh; k++) {
                                        // resolution
        sig[k] = hrn->GetRandom();
        
        float y = trn.Gaus(0,sig[k]);
        
        sy  += y;
        sy2 += y*y;
        sw  += 1./(sig[k]*sig[k]);
      }
        
      double w2m = sy2/nh;
      double wm  = sy /nh;

      chp.werr           = sqrt(1/sw);
      chp.sigw           = sqrt(w2m-wm*wm)*sqrt(nh/(nh-1.));
      chp.sigw_over_werr = chp.sigw/chp.werr;

      chp.sigws           = chp.sigw*_scale2;
      chp.sigws_over_werr = chp.sigws/chp.werr;
      
      fillComboHitHistograms(_hist.fComboHit[13],nullptr, &chp);  // all
    }
  }

//-----------------------------------------------------------------------------
  void  StrawHitRecoAna::fillEventHistograms(EventHist_t* Hist) {
    Hist->fEventNumber->Fill(_eventNum);

    Hist->fNch->Fill (_nComboHits);
    Hist->fNsh->Fill (_nStrawHits);
    Hist->fNssch->Fill(_nSingleSH );
  }

//-----------------------------------------------------------------------------
  void  StrawHitRecoAna::fillComboHitHistograms(ComboHitHist_t* Hist, const ComboHit* Ch, ComboHitPar_t* Chp) {

    Hist->fNsh->Fill(Chp->nsh);
    Hist->fWres->Fill(Chp->werr);
    Hist->fSigw->Fill(Chp->sigw);
    Hist->fSigwOverWres->Fill(Chp->sigw_over_werr);
    Hist->fSigws->Fill(Chp->sigws);
    Hist->fSigwsOverWres->Fill(Chp->sigws_over_werr);

    for (int i=0; i<Chp->nsh; i++) {
      Hist->fSSWres->Fill(Chp->ss_wres[i]);
    }
  }

//-----------------------------------------------------------------------------
// fill histograms
//-----------------------------------------------------------------------------
  void  StrawHitRecoAna::fillHistograms() {
//-----------------------------------------------------------------------------
// event histograms - just one set
//-----------------------------------------------------------------------------
    fillEventHistograms(_hist.fEvent[0]);
//-----------------------------------------------------------------------------
// straw hit histograms, 
//-----------------------------------------------------------------------------
    ComboHitPar_t chp;

    _n2 = 0; 
    _n3 = 0;


    for (int i=0; i<_nComboHits; i++) {
      const ComboHit* ch = &_chColl->at(i);
      int nsh = ch->nStrawHits();

      float sw(0), sw2(0);
      for (int ish=0; ish<nsh; ish++) {
        int ind = ch->index(ish);

        const ComboHit* sh = &_sschColl->at(ind);

        float wp = sh->wireDist();
        sw      += wp;
        sw2     += wp*wp;

        chp.ss_wres[ish] = sh->wireRes();
      }

      chp.nsh             = ch->nStrawHits();
      chp.werr            = ch->wireRes();
//-----------------------------------------------------------------------------
// resolution estimated based on the spread of the hits along the wire
//-----------------------------------------------------------------------------
      chp.sigw            = sqrt(sw2/nsh-sw*sw/nsh/nsh)*sqrt(nsh/(nsh-1+1.e-12));
      chp.sigw_over_werr  = chp.sigw/ch->wireRes();
      chp.sigws           = chp.sigw;
      chp.sigws_over_werr = chp.sigw_over_werr;

      fillComboHitHistograms(_hist.fComboHit[0],ch,&chp);  // all

      if (nsh > 1) fillComboHitHistograms(_hist.fComboHit[1],ch,&chp);
      if (nsh == 2) { 
        _n2 += 1;
        fillComboHitHistograms(_hist.fComboHit[2],ch,&chp);
      }

      if (nsh == 3) { 
        _n3 += 1;
        fillComboHitHistograms(_hist.fComboHit[3],ch,&chp);
      }
    }
  }

//-----------------------------------------------------------------------------
bool StrawHitRecoAna::findData(const art::Event& Evt) {
    _chColl     = nullptr;
    _sdmcColl   = nullptr;
    _chfColl    = nullptr;

    auto chcH   = Evt.getValidHandle<ComboHitCollection>(_chCollTag);
    _chColl     = chcH.product();
    _nComboHits = _chColl->size();

    auto shcH   = Evt.getValidHandle<StrawHitCollection>(_shCollTag);
    _shColl     = shcH.product();
    _nSingleSH  = _shColl->size();

    auto sschcH = Evt.getValidHandle<ComboHitCollection>(_sschCollTag);
    _sschColl   = sschcH.product();

    auto chfcH  = Evt.getValidHandle<StrawHitFlagCollection>(_chfCollTag);
    _chfColl    = chfcH.product();

    auto sdmccH = Evt.getValidHandle<StrawDigiMCCollection>(_sdmcCollTag);
    _sdmcColl   = sdmccH.product();

    return (_chColl != 0) && (_chfColl != 0) && (_sdmcColl != 0) ;
  }

//-----------------------------------------------------------------------------
  void StrawHitRecoAna::analyze(const art::Event& uEvent) {

    _eventNum = uEvent.event();
    if (_debugLevel) {
      printf("* >>> StrawHitRecoAna::%s event number: %10i\n",__func__,_eventNum);
    }

    _hlp->SetEvent(&uEvent);
//-----------------------------------------------------------------------------
// process event
//-----------------------------------------------------------------------------
    if (! findData(uEvent)) {
      throw cet::exception("RECO")
        << "mu2e::StrawHitRecoAna_module::produce: missing data" << endl;
    }
//-----------------------------------------------------------------------------
// in the end of event processing fill histograms
//-----------------------------------------------------------------------------
    fillHistograms  ();

    if (_debugLevel    > 0) debug();
  }

//-----------------------------------------------------------------------------
// debugLevel > 0: print seeds
//-----------------------------------------------------------------------------
  void StrawHitRecoAna::debug() {

    if (_printComboHits) {
      printf("* ComboHits \n");
//-----------------------------------------------------------------------------
// print ComboHits
//-----------------------------------------------------------------------------
      _hlp->printComboHitCollection(_chCollTag.encode().data(),
                                    _chfCollTag.encode().data(),
                                    _sdmcCollTag.encode().data());
    }

    if (_printSingleComboHits) {
      printf("* Single straw ComboHits tag:  %s\n",_shCollTag.encode().data());
//-----------------------------------------------------------------------------
// print ComboHits
//-----------------------------------------------------------------------------
      _hlp->printComboHitCollection(_shCollTag.encode().data(),
                                    _shfCollTag.encode().data(),
                                    _sdmcCollTag.encode().data());
    }
  }

// Part of the magic that makes this class a module.
DEFINE_ART_MODULE(StrawHitRecoAna)

}
