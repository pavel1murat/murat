//////////////////////////////////////////////////////////////////////////////
// test timing performance of sorting the combohit colection in time
// clones from anoother module, so comments could be wrong
//////////////////////////////////////////////////////////////////////////////
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art_root_io/TFileService.h"

#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"

#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
// conditions
#include "Offline/ConditionsService/inc/ConditionsHandle.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
// data
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"

#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/MCDataProducts/inc/StrawGasStep.hh"

#include "Offline/CalPatRec/inc/HlPrint.hh"

// diagnostics

// #include "Offline/Mu2eUtilities/inc/ModuleHistToolBase.hh"
// #include "art/Utilities/make_tool.h"

#include <algorithm>
#include <cmath>

using namespace std;

namespace mu2e {

  using CLHEP::Hep3Vector;

  class ComboHitSorter: public art::EDAnalyzer {
  public:
    struct Config {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag>     chCollTag  {Name("chCollTag"  ), Comment("ComboHit collection Name"    ) };
      fhicl::Atom<art::InputTag>     chfCollTag {Name("chfCollTag" ), Comment("StrawHitFlag collection Name") };
      fhicl::Atom<art::InputTag>     sdmcCollTag{Name("sdmcCollTag"), Comment("StrawDigiMC collection Name" ) };
      // fhicl::Atom<int>               simID      {Name("simID"      ), Comment("Selected sim particle ID"    ) };
      fhicl::Atom<int>               debugLevel {Name("debugLevel" ), Comment("debug level"                 ) };
    };
//-----------------------------------------------------------------------------
// talk-to parameters: input collections and algorithm parameters
//-----------------------------------------------------------------------------
    art::InputTag   _chfCollTag;
    art::InputTag   _chCollTag;
    art::InputTag   _sdmcCollTag;
    int             _simID;
    int             _debugLevel;
//-----------------------------------------------------------------------------
// cache event/geometry objects
//-----------------------------------------------------------------------------
    const ComboHitCollection*      _chColl  ;
    const StrawHitFlagCollection*  _chfColl ;  // combo hit flags
    const StrawDigiMCCollection*   _sdmcColl;

    int                            _nComboHits;
    std::vector<const ComboHit*>   _v;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  public:
    explicit ComboHitSorter(fhicl::ParameterSet const&);

  private:

    bool         findData     (const art::Event&  Evt);
//-----------------------------------------------------------------------------
// overloaded methods of the module class
//-----------------------------------------------------------------------------
    void beginJob()                          override;
    void beginRun(const art::Run&   ARun   ) override;
    void analyze (const art::Event& AnEvent) override;
  };

//-----------------------------------------------------------------------------
  ComboHitSorter::ComboHitSorter(fhicl::ParameterSet const& pset):
    art::EDAnalyzer{pset},
    _chfCollTag            (pset.get<string>("chfCollTag" )),
    _chCollTag             (pset.get<string>("chCollTag"  )),
    _sdmcCollTag           (pset.get<string>("sdmcCollTag")),
    _debugLevel            (pset.get<int>   ("debugLevel" ))
  {
    consumes<ComboHitCollection>    (_chCollTag  );
    consumes<StrawHitFlagCollection>(_chfCollTag );
    consumes<StrawDigiMCCollection> (_sdmcCollTag);

    // produces<ComboHitCollection>();
    // produces<StrawHitFlagCollection>();
  }

  //-----------------------------------------------------------------------------
  void ComboHitSorter::beginJob() {
    if (_debugLevel > 0) {
    }
  }

//-----------------------------------------------------------------------------
// create a Z-ordered map of the tracker
//-----------------------------------------------------------------------------
  void ComboHitSorter::beginRun(const art::Run& aRun) {
  }


//-----------------------------------------------------------------------------
  bool ComboHitSorter::findData(const art::Event& Evt) {

    auto chcH   = Evt.getValidHandle<ComboHitCollection>(_chCollTag);
    _chColl     = chcH.product();

    auto chfcH  = Evt.getValidHandle<StrawHitFlagCollection>(_chfCollTag);
    _chfColl    = chfcH.product();

    bool res = (_chColl != nullptr) and (_chfColl != nullptr) ;
    if (_debugLevel > 0) {
      auto sdmccH = Evt.getValidHandle<StrawDigiMCCollection>(_sdmcCollTag);
      _sdmcColl   = sdmccH.product();

      res  = res and (_sdmcColl != nullptr);

    }
    return  res;
  }

//-----------------------------------------------------------------------------
// Custom comparator to sort in descending order
//-----------------------------------------------------------------------------
  bool comparator(const ComboHit*& a, const ComboHit*& b) {
    return a->correctedTime() < b->correctedTime();
  }


//-----------------------------------------------------------------------------
  void ComboHitSorter::analyze(const art::Event& Event) {

    if (_debugLevel) printf(">>> ComboHitSorter::produce  event number: %10i\n",Event.event());
//-----------------------------------------------------------------------------
// process event
//-----------------------------------------------------------------------------
    if (! findData(Event)) {
      const char* message = "mu2e::ComboHitSorter_module::produce: data missing or incomplete";
      throw cet::exception("RECO")<< message << endl;
    }

    _nComboHits = _chColl->size();

    _v.clear();
    // _v.resize(_nComboHits);

    for (int i=0; i<_nComboHits; i++) {
      _v.push_back(&(*_chColl)[i]); //  = &(*_chColl)[i]; //  
    }

    std::sort(_v.begin(), _v.end(), comparator);
//-----------------------------------------------------------------------------
// debug, if requested
//-----------------------------------------------------------------------------
    if (_debugLevel > 0) {

      HlPrint* hlp = HlPrint::Instance();
      hlp->SetEvent(&Event);

      hlp->printComboHit(0,0,"banner",0,0);

      const ComboHit* ch0 = &(*_chColl)[0];

      for (int i=0; i<_nComboHits; i++) {
        int ind    = _v[i]-ch0;

        const StrawDigiMC*  sdmc = &_sdmcColl->at(ind);
        const StrawGasStep* sgs  = sdmc->earlyStrawGasStep().get();

        int flags = *((int*) &(*_chfColl)[i]);

        hlp->printComboHit(_v[i],sgs,"data",ind,flags);
      }
    }
  }
}
//-----------------------------------------------------------------------------
// magic that makes this class a module.
//-----------------------------------------------------------------------------
DEFINE_ART_MODULE(mu2e::ComboHitSorter)
//-----------------------------------------------------------------------------
// done
//-----------------------------------------------------------------------------
