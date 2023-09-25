
#include "murat/mod/TrkFragmentAna_module.hh"

namespace mu2e {

//-----------------------------------------------------------------------------
unsigned int reverseBits(unsigned int num) {
  unsigned int numOfBits = 10; // sizeof(num) * 8; // Number of bits in an unsigned int

  unsigned int reversedNum = 0;
  for (unsigned int i = 0; i < numOfBits; ++i) {
    if ((num & (1 << i)) != 0)
      reversedNum |= 1 << ((numOfBits - 1) - i);
  }

  return reversedNum;
}

// ======================================================================

  TrkFragmentAna::TrkFragmentAna(fhicl::ParameterSet const& PSet) : 
    THistModule     (PSet                  ,"TrkFragmentAna"),
    _diagLevel      (PSet.get<int>         ("diagLevel")), 
    _trkFragmentsTag(PSet.get<std::string> ("trkTag")) {
  }

//--------------------------------------------------------------------------------
// book histograms
//--------------------------------------------------------------------------------
void TrkFragmentAna::beginJob() {
  art::ServiceHandle<art::TFileService> tfs;

  art::TFileDirectory top_dir  = tfs->mkdir("trk");
  art::TFileDirectory frag_dir = tfs->mkdir("trk/frag_0");

  _Hist.event.nbtot       = top_dir.make<TH1F>("nbtot"   , "nbytes total" ,10000, 0., 10000.);
  _Hist.event.nfrag       = top_dir.make<TH1F>("nfrag"   , "n fragments"  ,  100, 0.,   100.);

  _Hist.frag.nbytes       = frag_dir.make<TH1F>("nbytes"  , "n bytes"     ,10000, 0., 10000.);
  _Hist.frag.npackets     = frag_dir.make<TH1F>("npackets", "n packets"   , 1000, 0.,  1000.);
  _Hist.frag.nhits        = frag_dir.make<TH1F>("nhits"   , "n hits"      ,  300, 0.,   300.);
  _Hist.frag.valid        = frag_dir.make<TH1F>("valid"   , "valid"       ,    2, 0.,     2.);

  // _hTrkStrawId   = trkDir.make<TH1F>("h_sid"     , " strawId", 20000, 0., 20000.);
  // _hTrkTDC[0]    = trkDir.make<TH1F>("h_tdc0"    , "TDC0; TDC[0]", 264, 0., 264000);
  // _hTrkTDC[1]    = trkDir.make<TH1F>("h_tdc1"    , "TDC1; TDC[1]", 264, 0., 264000);
  // _hTrkTDC[2]    = trkDir.make<TH1F>("h_tdcm"    , "(TDC[0]+TDC[1])/2",2000, 0., 20000.);
  // _hTrkTDC[3]    = trkDir.make<TH1F>("h_dtdc"    , "TDC[1]-TDC[0]", 220, -100., 10000.);
  // _hTrkTOT       = trkDir.make<TH1F>("h_tot"     , "(TOT[0]+TOT[1])/2", 100, 0., 200.);
  // _hTrkPMP       = trkDir.make<TH1F>("h_pmp"     , "trk fragment average PMP; PMP", 100, 0., 200.);
  // _hTrkMeanADC   = trkDir.make<TH1F>("h_mean_adc", "trk fragment Mean ADC; <ADC>", 250, 0., 2500.);
  // _hTrkMaxADC    = trkDir.make<TH1F>("h_max_adc" , "trk fragment Max ADC; Max_ADC", 250, 0., 2500.);
  // _hTrkWfSize    = trkDir.make<TH1F>("h_nsamples", "waveform N(samples)",  20, 0., 20.);

  printf("[mu2e::TrkFragmentAna] pointer to the module: 0x%8p\n",(void*) this);

  // for (int i=0; i<100; i++) {
  //   _Hist.channel[i].h_time[0] = trkDir.make<TH1F>(Form("h_time_ch_%02i_0",i),Form("Channel %02i time[0]",i),250, 0., 25000.);
  //   _Hist.channel[i].h_time[1] = trkDir.make<TH1F>(Form("h_time_ch_%02i_1",i),Form("Channel %02i time[1]",i),250, 0., 25000.);
  //   _Hist.channel[i].h_nwf     = trkDir.make<TH1F>(Form("h_nwf_ch_%02i"   ,i),Form("Channel %02i nwf"    ,i), 50,  0.,    50.);
  //   _Hist.channel[i].h_pmp     = trkDir.make<TH1F>(Form("h_pmp_ch_%02i"   ,i),Form("Channel %02i pmp"    ,i), 100, 0.,  500.);

  //   for (int j=0; j<10; j++) {
  //     _Hist.h_wf [i][j] = trkDir.make<TH1F>(Form("h_wf_ch_%02i_%i",i,j),Form("Channel [%02i][%i] waveform",i,j),20, 0.,20.);
  //   }
  // }
}

//-----------------------------------------------------------------------------
void TrkFragmentAna::endJob() {}

//--------------------------------------------------------------------------------
// assume that we only have tracker fragment(s)
//-----------------------------------------------------------------------------
void TrkFragmentAna::analyze(const art::Event& event) {
  //art::EventNumber_t eventNumber = event.event();

  int    nbtot = 0;
  int    nfrag = 0;

  for (int i=0; i<100; i++) _nwf[i] = 0;

  auto handle = event.getValidHandle<std::vector<artdaq::Fragment> >(_trkFragmentsTag);
//-----------------------------------------------------------------------------
// calculate the fragment size manually - big thank you to designers (:
//----------------------------------------------------------------------------- 
  for (/* auto */ const artdaq::Fragment& frag : *handle) {
    ushort* buf = (ushort*) (frag.dataBegin());
    int nbytes  = buf[0];
    nfrag      += 1;
    nbtot      += nbytes;

    analyze_fragment(&frag,&_Hist.frag);

    if (_diagLevel > 2) {
//-----------------------------------------------------------------------------
// print fragments in HEX, for the tracker, the data has to be in 2-byte words
//-----------------------------------------------------------------------------
      int nw      = nbytes/2;
      int loc     = 0;
      
      for (int i=0; i<nw; i++) {
        if (loc == 0) printf(" 0x%08x: ",i*2);
        
        ushort  word = buf[i];
        printf("0x%04x ",word);
        
        loc += 1;
        if (loc == 8) {
          printf("\n");
          loc = 0;
        }
      }
      
      if (loc != 0) printf("\n");
    }
  }

  _Hist.event.nbtot->Fill(nbtot);
  _Hist.event.nfrag->Fill(nfrag);

  if (_diagLevel > 1) {
    printf(" Run : %5i subrun: %5i event: %8i nfrag: %3i nbytes: %5i\n", 
           event.run(),event.subRun(),event.event(), nfrag, nbtot);
  }
//-----------------------------------------------------------------------------
// go into interactive mode, 
// fInteractiveMode = 0 : do not stop
// fInteractiveMode = 1 : stop after each event (event display mode)
// fInteractiveMode = 2 : stop only in the end of run, till '.q' is pressed
//-----------------------------------------------------------------------------
  TModule::analyze(event);
}

//-----------------------------------------------------------------------------
// void TrkFragmentAna::analyze_tracker(const mu2e::TrackerFragment& Fragment) {
  void TrkFragmentAna::analyze_fragment(const artdaq::Fragment* Fragment, FragmentHist_t* Hist) {


  int nbytes = Fragment->size();

  short* data = (short*) Fragment->dataBegin();

  DtcDataHeaderPacket_t* dh = (DtcDataHeaderPacket_t*) (data+0x20);

  int npackets = dh->nPackets;

  Hist->nbytes->Fill(nbytes);
  Hist->npackets->Fill(npackets);

//   for (size_t curBlockIdx = 0; curBlockIdx < Fragment.block_count(); curBlockIdx++) {

//     auto block = Fragment.dataAtBlockIndex(curBlockIdx);
    
//     if (block == nullptr) {
//       mf::LogError("TrkFragmentAna") << "Unable to retrieve header from block " 
//                                      << curBlockIdx << "!"
//                                      << std::endl;
//       continue;
//     }
//     auto hdr = block->GetHeader();
// //-----------------------------------------------------------------------------
// // hdr is the DTC data packer header
// //-----------------------------------------------------------------------------
//     if (hdr->GetPacketCount() > 0) {

//       auto trkDataVec = Fragment.GetTrackerData(curBlockIdx);
//       if (trkDataVec.empty()) {
//         mf::LogError("TrkFragmentAna") << "Error retrieving Tracker data from DataBlock "
//                                     << curBlockIdx << "! Aborting processing of this block!";
//         continue;
//       }
// //-----------------------------------------------------------------------------
// // trkDataVec is a list of hits
// //-----------------------------------------------------------------------------
//       for (auto& trkDataPair : trkDataVec) {
//         mu2e::StrawId sid(trkDataPair.first->StrawIndex);
//         mu2e::TrkTypes::TDCValues tdc = {trkDataPair.first->TDC0(), trkDataPair.first->TDC1()};
//         mu2e::TrkTypes::TOTValues tot = {trkDataPair.first->TOT0, trkDataPair.first->TOT1};
//         mu2e::TrkTypes::ADCValue  pmp = trkDataPair.first->PMP;
//         int sum{0};
//         unsigned short maxadc{0};
//         for (auto adc : trkDataPair.second) {
//           sum += adc;
//           maxadc = std::max(maxadc, adc);
//         }

//         _hTrkStrawId->Fill(sid.asUint16());
//         _hTrkTDC[0]->Fill(tdc[0]);
//         _hTrkTDC[1]->Fill(tdc[1]);
//         _hTrkTDC[2]->Fill((tdc[0] + tdc[1]) / 2.);
//         _hTrkTDC[3]->Fill(tdc[1] - tdc[0]);
//         _hTrkTOT->Fill((tot[0] + tot[1]) / 2.);
//         _hTrkPMP->Fill(pmp);
//         int mean = (trkDataPair.second.size() != 0) ? sum / trkDataPair.second.size() : -1.;
//         _hTrkMeanADC->Fill(mean);
//         _hTrkMaxADC->Fill(maxadc);
// //-----------------------------------------------------------------------------
// // apparently, the waveform is the second..
// //-----------------------------------------------------------------------------
//         int nsamples = trkDataPair.second.size();

//         _hTrkWfSize->Fill(nsamples);

//         int ch = sid.straw();
// //-----------------------------------------------------------------------------
// // deal with the bug in the channel numbering
// //-----------------------------------------------------------------------------
//         if (ch > 128) ch = ch-128;

//         float t0 = trkDataPair.first->TDC0()*5/256.;
//         float t1 = trkDataPair.first->TDC1()*5/256.;

//         _Hist.channel[ch].h_time[0]->Fill(t0);
//         _Hist.channel[ch].h_time[1]->Fill(t1);
//         _Hist.channel[ch].h_pmp->Fill(pmp);

//         int nw = _nwf[ch];

//         _Hist.h_wf[ch][nw]->Reset();
//         for (int i=0; i<nsamples; i++) {
//           uint adc = trkDataPair.second.at(i);

//           uint adc1 = reverseBits(adc);

//           _Hist.h_wf[ch][nw]->Fill(i+1,adc1);
//         }

//         _nwf[ch] += 1;
//       }

//       for (int ch=0; ch<100; ch++) {
//         _Hist.channel[ch].h_nwf->Fill(_nwf[ch]);
//       }
//     }
//   }

  // Fragment.ClearUpgradedPackets();
}

} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::TrkFragmentAna)
