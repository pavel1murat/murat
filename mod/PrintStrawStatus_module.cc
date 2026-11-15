///////////////////////////////////////////////////////////////////////////////
// print the long-term and the short-term status of every tracker straw, once per run
///////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"

#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/DataProducts/inc/StrawStatus.hh"
#include "Offline/DbService/inc/DbHandle.hh"
#include "Offline/DbTables/inc/TrkElementStatus.hh"

#include <cstdint>
#include <format>
#include <iostream>
#include <map>
#include <memory>
#include <string>

namespace mu2e {

class PrintStrawStatus : public art::EDAnalyzer {
public:
  struct Config {
    using Name    = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Atom<bool> onlyFlagged{Name("onlyFlagged"), Comment("print only the straws with a status flag")};
  };
  using Parameters = art::EDAnalyzer::Table<Config>;

  explicit PrintStrawStatus(Parameters const& config);

  void analyze(art::Event const& event) override;

private:
  void printRun(art::Event const& event);

  bool                                           onlyFlagged_;
  unsigned int                                   runNumber_{0};
  bool                                           haveRun_{false};
  std::unique_ptr<DbHandle<TrkStrawStatusLong>>  _tssl_p;
  std::unique_ptr<DbHandle<TrkStrawStatusShort>> _tsss_p;
};

//-----------------------------------------------------------------------------
PrintStrawStatus::PrintStrawStatus(Parameters const& config)
  : art::EDAnalyzer(config)
  , onlyFlagged_(config().onlyFlagged())
  , _tssl_p(std::make_unique<DbHandle<TrkStrawStatusLong>>())
  , _tsss_p(std::make_unique<DbHandle<TrkStrawStatusShort>>())
{
}

//-----------------------------------------------------------------------------
void PrintStrawStatus::analyze(art::Event const& event) {
  if (haveRun_ and event.run() == runNumber_) return;
  runNumber_ = event.run();
  haveRun_   = true;
  printRun(event);
}

//-----------------------------------------------------------------------------
void PrintStrawStatus::printRun(art::Event const& event) {
  auto const& longTable  = _tssl_p->get(event.id());
  auto const& shortTable = _tsss_p->get(event.id());

  // status of the straws in the tables, keyed by the StrawId
  std::map<uint16_t,std::string> longStatus;
  std::map<uint16_t,std::string> shortStatus;
  for (auto const& row : longTable.rows())  longStatus [row.id().asUint16()] = row.status().stringRep();
  for (auto const& row : shortTable.rows()) shortStatus[row.id().asUint16()] = row.status().stringRep();

  std::cout << std::format("straw status for run {}: long-term rows {}, short-term rows {}\n",
                           runNumber_,longTable.rows().size(),shortTable.rows().size());

  for (int plane = 0; plane < StrawId::_nplanes; ++plane) {
    for (int panel = 0; panel < StrawId::_npanels; ++panel) {
      for (int straw = 0; straw < StrawId::_nstraws; ++straw) {
        StrawId sid(plane,panel,straw);
        auto lit = longStatus.find(sid.asUint16());
        auto sit = shortStatus.find(sid.asUint16());
        std::string lstat = (lit == longStatus.end()) ? "" : lit->second;
        std::string sstat = (sit == shortStatus.end()) ? "" : sit->second;

        if (onlyFlagged_ and lstat.empty() and sstat.empty()) continue;

        std::cout << std::format("plane {:2d} panel {} straw {:2d} long: {:<40} short: {}\n",
                                 plane,panel,straw,
                                 lstat.empty() ? std::string("-") : lstat,
                                 sstat.empty() ? std::string("-") : sstat);
      }
    }
  }
}

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::PrintStrawStatus)
