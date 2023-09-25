// ======================================================================
//
// ======================================================================
#ifndef __murat_mod_TrkFragmentAna_hh__
#define __murat_mod_TrkFragmentAna_hh__

// ROOT includes
#include "TH1F.h"
//#include "TFolder.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/ParameterSet.h"

#ifndef __CLING__ 
#include "artdaq-core-mu2e/Overlays/FragmentType.hh"
#include "artdaq-core-mu2e/Data/TrackerFragment.hh"
#include "artdaq-core/Data/Fragment.hh"
#else 
namespace mu2e {
  class TrackerFragment;
}

namespace artdaq {
  class Fragment;
}
#endif

// Mu2e includes
#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/DataProducts/inc/TrkTypes.hh"

// #include "Stntuple/print/TAnaDump.hh"
#include "Stntuple/mod/THistModule.hh"

#include <iostream>
#include <memory>
#include <string>

// ======================================================================
namespace mu2e {

  class TrkFragmentAna : public THistModule {

  public:
    struct DtcDMAPacket_t {              // 2 16-byte words
      uint16_t       byteCount   : 16;   ///< Byte count of current block
      uint8_t        subsystemID :  4;   ///< Hop count
      uint16_t       packetType  :  4;   ///< Packet type               DTC_PacketType
      uint16_t       ROCID       :  3;   ///< Link identifier of packet DTC_LinkID
      uint16_t       unused      :  4;
      bool           valid       :  1;   ///< Whether the DTC believes the packet to be valid
    };

    struct DtcDataHeaderPacket_t : public DtcDMAPacket_t {  // 8 16-byte words in tota
      uint16_t            nPackets     : 11;
      uint16_t            unused       :  5;
      uint16_t            eventTag[3];             // DTC_EventWindowTag
      uint8_t             status;                  // DTC_DataStatus
      uint8_t             version;
      uint8_t             DTCID;
      uint8_t             EVBMode;
    };

    struct DtcDataBlock_t : public DtcDataHeaderPacket_t {
      uint16_t            hitData[10000];
    };

    DtcDataBlock_t*  _trkFragment;
    int              _diagLevel;

    art::InputTag    _trkfCollTag;

    struct EventHist_t {
      TH1F*         nbtot;
      TH1F*         nfrag;
    };

    //   TH1F*         _hTrkNFragment;
    //   TH1F*         _hTrkStrawId;
    //   TH1F*         _hTrkTDC[4];
    //   TH1F*         _hTrkTOT;
    //   TH1F*         _hTrkPMP;
    //   TH1F*         _hTrkMeanADC;
    //   TH1F*         _hTrkMaxADC;
    //   TH1F*         _hTrkWfSize;
    // };

    struct FragmentHist_t {
      TH1F*         nbytes;
      TH1F*         npackets;
      TH1F*         nhits;
      TH1F*         valid;
    };
                                        // assume one panel
    struct ChannelHist_t {
      TH1F*         h_time[2];
      TH1F*         h_nwf;
      TH1F*         h_pmp;
    };
                                        // assume one panel
    struct WaveformHist_t {
      ChannelHist_t     channel[100];
      TH1F*             h_wf   [100][10];
    };

    struct Hist_t {
      EventHist_t       event;
      FragmentHist_t    frag;
    } _Hist;

    int _nwf[100];

    // struct Config {
    //   fhicl::Atom<int>               diagLevel{fhicl::Name("diagLevel"), fhicl::Comment("diagnostic level")};
    //   fhicl::Atom<int>               parseTRK {fhicl::Name("parseTRK" ), fhicl::Comment("parseTRK"        )};
    //   fhicl::Atom<art::InputTag>     trkTag   {fhicl::Name("trkTag"   ), fhicl::Comment("trkTag"          )};
    //   fhicl::Table<TAnaDump::Config> tAnaDump {fhicl::Name("tAnaDump" ), fhichl::Comment("Diag plugin"    )};
    // };
    
    explicit TrkFragmentAna(fhicl::ParameterSet const& pset);
    // explicit TrkFragmentAna(const art::EDAnalyzer::Table<Config>& config);
    virtual ~TrkFragmentAna() {}
    
    virtual void beginJob() override;
    virtual void endJob  () override;
    
    virtual void analyze        (const art::Event& e) override;
    //    void         analyze_tracker(const mu2e::TrackerFragment& cc);
    void         analyze_fragment(const artdaq::Fragment* Fragment,FragmentHist_t* Hist);
  };
}
#endif
