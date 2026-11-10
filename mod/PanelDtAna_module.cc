// art implementation of the STNTUPLE TPanelDtAnaModule.
// Histogram directory names and the sel/slot/panel/channel hierarchy match
// the original module so existing ROOT-based plotting scripts can be reused.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Name.h"
#include "fhiclcpp/types/Comment.h"

#include "Offline/DataProducts/inc/StrawEnd.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"

#include "murat/alg/RunInfo_t.hh"

#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "TRACE/tracemf.h"
#define  TRACE_NAME "PanelDtAna_module"
namespace mu2e {

class PanelDtAna : public art::EDAnalyzer {
public:
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    fhicl::Atom<art::InputTag> shCollTag{Name("shCollTag"), Comment("StrawHit collection")};
    fhicl::Atom<art::InputTag> chCollTag{Name("chCollTag"), Comment("ComboHit collection")};
    fhicl::Atom<art::InputTag> tcCollTag{Name("tcCollTag"), Comment("TimeCluster collection")};
    //    fhicl::Atom<art::InputTag> trackCollTag{Name("trackCollTag"), Comment("KalSeed collection")};
  };
  using Parameters = art::EDAnalyzer::Table<Config>;

  static constexpr int nSlots = 18;
  static constexpr int nPanelsPerSlot = 12;
  static constexpr int nChannels = 96;
  static constexpr int nPlanes = 36;
  static constexpr int nPanels = 216;

  struct ChannelHist {
    TH1F *tdc0{}, *dt10{}, *edep{};
  };
  
  struct PanelHist {
    std::array<ChannelHist, nChannels> ch{};
    TH1F *occup{}, *edep{};
  };
  
  struct SlotHist {
    std::array<PanelHist, nPanelsPerSlot> panel{};
    TH1F* edep{};
  };

  struct SelHist {
    TH1F *occup{};
    TH2F *occup2d{}, *strawVsPanel{}, *panelVsEvent{};
    std::array<SlotHist, nSlots> slot{};
  };
  
  struct TcHist {
    TH1F *nch{}, *ngh{}, *nplanes{}, *npanels{}, *maxNhPanel{}, *maxNhPlane{};
    TH2F *edepVsPanel{}, *dtVsPanel{};
    std::array<std::array<TH1F*, nPlanes>, nPlanes> dtPlane{};
    std::array<TH2F*, nPanels> edep_vs_ewt_panel{};
  };
  
  struct EventHist {
    TH1F *nch{}, *nsh{}, *ntc{}, *ntrk{}, *nshVsEvent{};
  };
  
  struct Hist {
    EventHist event;
    std::array<SelHist,2> sel{};
    std::array<TcHist ,2> tc {};
  } _hist;


  int fRunNumber{-1};

  murat::RunInfo_t fRunInfo;

  explicit PanelDtAna(const Parameters& config);
  void beginJob() override;
  void beginRun(const art::Run& run) override;
  void analyze(const art::Event& event) override;

private:
  art::InputTag _shTag, _chTag, _tcTag, _trackTag;
  void bookHistograms();
  void bookSelection(SelHist& hist, int selection, art::TFileDirectory& dir);
  void bookSlot(SlotHist& hist, int slot, art::TFileDirectory& dir);
  void bookPanel(PanelHist& hist, int slot, int pnl12, art::TFileDirectory& dir);
  void bookChannel(ChannelHist& hist, int slot, int plane, int panel, int channel,
                   art::TFileDirectory& dir);
  void bookTimeCluster(TcHist& hist, int selection, art::TFileDirectory& dir);
  void fillSelection(SelHist& hist, int selection, int slot, int pnl12,
                     int plane, int panel, int channel, float edep,
                     float tdc0, float dt10, double eventNumber);
};

PanelDtAna::PanelDtAna(const Parameters& config)
  : art::EDAnalyzer{config},
    _shTag{config().shCollTag()}, _chTag{config().chCollTag()},
    _tcTag{config().tcCollTag()}
    // , _trackTag{config().trackCollTag()}
{
  consumes<StrawHitCollection>(_shTag);
  consumes<ComboHitCollection>(_chTag);
  consumes<TimeClusterCollection>(_tcTag);
  // consumes<KalSeedCollection>(_trackTag);
  TH1::AddDirectory(false);
}

//-----------------------------------------------------------------------------
  void PanelDtAna::beginRun(const art::Run& run) {

    int rn = run.run();
    TLOG(TLVL_INFO) << std::format(" BEGIN RUN entered: rn:{}",rn);
    
    if (rn != fRunNumber) {
      fRunNumber = rn;
      fRunInfo.GetFromRunDb(rn,0);
      bookHistograms();
    }
  }

  //-----------------------------------------------------------------------------
  void PanelDtAna::beginJob() {
  }

void PanelDtAna::bookHistograms() {
  art::ServiceHandle<art::TFileService> tfs;
  auto evt = tfs->mkdir("evt_00");
  _hist.event.nch = evt.make<TH1F>("nch", "N(combo hits)", 200, 0, 2000);
  _hist.event.nsh = evt.make<TH1F>("nsh", "N(straw hits)", 200, 0, 2000);
  _hist.event.ntc = evt.make<TH1F>("ntc", "N(time clusters)", 100, 0, 100);
  //  _hist.event.ntrk = evt.make<TH1F>("ntrk", "N(tracks)", 10, 0, 10);
  
  int nb  = 10000;
  int bin  = fRunInfo.max_ewtag/nb+1;
  _hist.event.nshVsEvent = evt.make<TH1F>("nsh_vs_ewt", "N(straw hits) vs event",nb,0,nb*bin);

  for (int i = 0; i < 2; ++i) {
    auto selectionDir = tfs->mkdir(Form("sel_%02d", i));
    bookSelection(_hist.sel[i], i, selectionDir);
  }
  for (int i = 0; i < 2; ++i) {
    auto tcDir = tfs->mkdir(Form("tc_%02d", i));
    bookTimeCluster(_hist.tc[i], i, tcDir);
  }
}

void PanelDtAna::bookSelection(SelHist& h, int selection, art::TFileDirectory& dir) {
  const std::string prefix = "selection " + std::to_string(selection);
  h.occup2d = dir.make<TH2F>("occup_2d", (prefix + " : occupancy, panel vs ch").c_str(),
                             100, 0, 100, nPanels, 0, nPanels);
  h.occup = dir.make<TH1F>("occup", (prefix + " : per panel occupancy").c_str(), nPanels, 0, nPanels);
  h.strawVsPanel = dir.make<TH2F>("straw_vs_pnl", (prefix + " : straw vs panel occupancy").c_str(),
                                  nPanels, 0, nPanels, nChannels, 0, nChannels);
  
  int nb  = 1000;
  int bin  = fRunInfo.max_ewtag/nb+1;
  h.panelVsEvent = dir.make<TH2F>("pnl_vs_ewt", (prefix + " : panel vs event").c_str(),nb, 0, nb*bin, nPanels, 0, nPanels);
  for (int slot = 0; slot < nSlots; ++slot) {
    auto slotDir = dir.mkdir(Form("slot_%02d", slot));
    bookSlot(h.slot[slot], slot, slotDir);
  }
}

//-----------------------------------------------------------------------------
void PanelDtAna::bookSlot(SlotHist& h, int slot, art::TFileDirectory& dir) {
  h.edep = dir.make<TH1F>("edep", Form("slot %02d : edep, MeV", slot), 200, -0.001, 0.009);
  for (int pnl12 = 0; pnl12 < nPanelsPerSlot; ++pnl12) {
    auto panelDir = dir.mkdir(Form("pnl_%03d", pnl12));
    bookPanel(h.panel[pnl12], slot, pnl12, panelDir);
  }
}

//-----------------------------------------------------------------------------
void PanelDtAna::bookPanel(PanelHist& h, int slot, int pnl12, art::TFileDirectory& dir) {
  const int plane = 2 * slot + pnl12 / 6;
  const int panel = pnl12 % 6;
  h.occup = dir.make<TH1F>("occup", Form("slot %02d plane %02d panel %03d : occupancy", slot, plane, panel),
                           100, 0, 100);
  h.edep = dir.make<TH1F>("edep", "energy deposition, MeV", 200, -0.001, 0.009);
  for (int channel = 0; channel < nChannels; ++channel) {
    auto channelDir = dir.mkdir(Form("ch_%02d", channel));
    bookChannel(h.ch[channel], slot, plane, panel, channel, channelDir);
  }
}

void PanelDtAna::bookChannel(ChannelHist& h, int slot, int plane, int panel, int channel,
                             art::TFileDirectory& dir) {
  const std::string prefix = Form("slot %02d plane %02d panel %02d channel %02d", slot, plane, panel, channel);
  h.tdc0 = dir.make<TH1F>("tdc0", (prefix + " : tdc0, ns").c_str(), 1500, 0, 150);
  h.dt10 = dir.make<TH1F>("dt10", (prefix + " : dt10, ns").c_str(), 1000, -50, 50);
  h.edep = dir.make<TH1F>("edep", (prefix + " : edep, MeV").c_str(), 200, -0.001, 0.009);
}


//-----------------------------------------------------------------------------
void PanelDtAna::bookTimeCluster(TcHist& h, int selection, art::TFileDirectory& dir) {
  const std::string prefix = "selection " + std::to_string(selection);
  h.nch = dir.make<TH1F>("nch", (prefix + " : N(combo hits)").c_str(), 200, 0, 200);
  h.ngh = dir.make<TH1F>("ngh", (prefix + " : N(good hits)").c_str(), 200, 0, 200);
  h.nplanes = dir.make<TH1F>("nplanes", (prefix + " : N(planes)").c_str(), 36, 0, 36);
  h.npanels = dir.make<TH1F>("npanels", (prefix + " : N(panels)").c_str(), 250, 0, 250);
  h.maxNhPlane = dir.make<TH1F>("max_nh_plane", "maximum hits per plane", 250, 0, 250);
  h.maxNhPanel = dir.make<TH1F>("max_nh_panel", "maximum hits per panel", 100, 0, 100);
  h.edepVsPanel = dir.make<TH2F>("edep_vs_pnl", "mean energy deposition vs panel", nPanels, 0, nPanels,
                                 200, -0.0001, 0.0099);
  h.dtVsPanel = dir.make<TH2F>("dt_vs_pnl", "T(panel)-T0(TC) vs panel, ns", nPanels, 0, nPanels,
                               400, -100, 100);
  for (int p2 = 1; p2 < nPlanes; ++p2) {
    for (int p1 = 0; p1 < p2; ++p1) {
      h.dtPlane[p2][p1] = dir.make<TH1F>(Form("dt_pln_%02d%02d", p2, p1),
        Form("%s : time cluster T[%02d]-T[%02d]", prefix.c_str(), p2, p1), 1000, -500, 500);
    }
  }

  int nb  = 1000;
  int bin  = fRunInfo.max_ewtag/nb+1;
  
  for (int panel = 0; panel < nPanels; ++panel) {
    h.edep_vs_ewt_panel[panel] = dir.make<TH2F>(Form("edep_vs_ewt_pnl_%03d", panel),
                                                Form("%s : panel %03d TC edep vs event", prefix.c_str(), panel),
                                                nb, 0, nb*bin, 200, -0.001, 0.009);
  }
}

void PanelDtAna::fillSelection(SelHist& h, int, int slot, int pnl12, int plane, int panel,
                               int channel, float edep, float tdc0, float dt10, double eventNumber) {
  const int offlinePanel = plane * 6 + panel;
  h.occup->Fill(offlinePanel);
  h.occup2d->Fill(channel, offlinePanel);
  h.strawVsPanel->Fill(offlinePanel, channel);
  h.panelVsEvent->Fill(eventNumber, offlinePanel);
  auto& slotHist = h.slot.at(slot);
  slotHist.edep->Fill(edep);
  auto& panelHist = slotHist.panel.at(pnl12);
  panelHist.occup->Fill(channel);
  panelHist.edep->Fill(edep);
  auto& channelHist = panelHist.ch.at(channel);
  channelHist.tdc0->Fill(tdc0);
  channelHist.dt10->Fill(dt10);
  channelHist.edep->Fill(edep);
}

void PanelDtAna::analyze(const art::Event& event) {
  const auto& strawHits = *event.getValidHandle<StrawHitCollection>(_shTag);
  const auto& comboHits = *event.getValidHandle<ComboHitCollection>(_chTag);
  const auto& timeClusters = *event.getValidHandle<TimeClusterCollection>(_tcTag);
  // const auto& tracks = *event.getValidHandle<KalSeedCollection>(_trackTag);
  const double eventNumber = event.event();

  _hist.event.nsh->Fill(strawHits.size());
  _hist.event.nch->Fill(comboHits.size());
  _hist.event.ntc->Fill(timeClusters.size());
  // _hist.event.ntrk->Fill(tracks.size());
  _hist.event.nshVsEvent->Fill(eventNumber, strawHits.size());

  for (const auto& hit : strawHits) {
    const auto sid = hit.strawId();
    const int plane = sid.plane();
    const int panel = sid.panel();
    const int channel = sid.straw();
    const int slot = plane / 2;
    const int pnl12 = 6 * (plane % 2) + panel;
    if (slot < 0 || slot >= nSlots || pnl12 < 0 || pnl12 >= nPanelsPerSlot ||
        channel < 0 || channel >= nChannels) continue;
    const float edep = hit.energyDep();
    const float tdc0 = hit.time(StrawEnd::cal);
    const float dt10 = hit.time(StrawEnd::hv) - hit.time(StrawEnd::cal);
    fillSelection(_hist.sel[0], 0, slot, pnl12, plane, panel, channel, edep, tdc0, dt10, eventNumber);
    if (edep > 0.0005f)
      fillSelection(_hist.sel[1], 1, slot, pnl12, plane, panel, channel, edep, tdc0, dt10, eventNumber);
  }

  for (const auto& tc : timeClusters) {
    std::array<int, nPlanes> nhPlane{};
    std::array<double, nPlanes> timePlane{};
    std::array<int, nPanels> nhPanel{};
    std::array<double, nPanels> timePanel{};
    std::array<double, nPanels> edepPanel{};
    int goodHits = 0;

    for (const auto hitIndex : tc._strawHitIdxs) {
      if (hitIndex >= comboHits.size()) continue;
      const auto& hit = comboHits[hitIndex];
      const auto sid = hit.strawId();
      const int plane = sid.plane();
      const int panel = sid.panel();
      if (plane < 0 || plane >= nPlanes || panel < 0 || panel >= 6) continue;
      const int panelIndex = plane * 6 + panel;
      if (panelIndex >= nPanels) continue;
      ++nhPlane[plane];
      timePlane[plane] += hit.correctedTime();
      ++nhPanel[panelIndex];
      timePanel[panelIndex] += hit.correctedTime();
      edepPanel[panelIndex] += hit.energyDep();
      if (hit.energyDep() > 0.0005f) ++goodHits;
    }

    int nPlanesWithHits = 0, nPanelsWithHits = 0, maxNhPlane = 0, maxNhPanel = 0;
    for (int i = 0; i < nPlanes; ++i) {
      if (nhPlane[i] > 0) ++nPlanesWithHits;
      maxNhPlane = std::max(maxNhPlane, nhPlane[i]);
      if (nhPlane[i] > 1) timePlane[i] /= nhPlane[i];
    }
    for (int i = 0; i < nPanels; ++i) {
      if (nhPanel[i] > 0) ++nPanelsWithHits;
      maxNhPanel = std::max(maxNhPanel, nhPanel[i]);
      if (nhPanel[i] > 1) {
        timePanel[i] /= nhPanel[i];
        edepPanel[i] /= nhPanel[i];
      }
    }
//-----------------------------------------------------------------------------
// fill time cluster histograms
//-----------------------------------------------------------------------------
    auto fillTc = [&](TcHist& h) {
      h.nch->Fill(tc._strawHitIdxs.size());
      h.ngh->Fill(goodHits);
      h.nplanes->Fill(nPlanesWithHits);
      h.npanels->Fill(nPanelsWithHits);
      h.maxNhPlane->Fill(maxNhPlane);
      h.maxNhPanel->Fill(maxNhPanel);
      for (int p2 = 1; p2 < nPlanes; ++p2) {
        if (nhPlane[p2] < 2) continue;
        for (int p1 = 0; p1 < p2; ++p1) {
          if (nhPlane[p1] >= 2) h.dtPlane[p2][p1]->Fill(timePlane[p2] - timePlane[p1]);
        }
      }
      for (int panelIndex = 0; panelIndex < nPanels; ++panelIndex) {
        if (nhPanel[panelIndex] < 2) continue;
        h.edepVsPanel->Fill(panelIndex, edepPanel[panelIndex]);
        h.edep_vs_ewt_panel[panelIndex]->Fill(eventNumber, edepPanel[panelIndex]);
        if (nPanelsWithHits > 10) h.dtVsPanel->Fill(panelIndex, timePanel[panelIndex] - tc.t0().t0());
      }
    };
    fillTc(_hist.tc[0]);
    if (tc._strawHitIdxs.size() >= 8 && goodHits >= 2) fillTc(_hist.tc[1]);
  }
}

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::PanelDtAna)
