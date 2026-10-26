///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <format>

#include "Stntuple/obj/TCrvDigi.hh"
#include "murat/gui/TEvdCrvCounter.hh"
#include "TCanvas.h"
#include "TH1F.h"

#include "TRACE/trace.h"
#define TRACE_NAME "TEvdCrvCounter"


ClassImp(murat::TEvdCrvCounter)

namespace murat {

  TCanvas* TEvdCrvCounter::fgCanvasWf = nullptr;

//-----------------------------------------------------------------------------
TEvdCrvCounter::TEvdCrvCounter(const char* Name, TGeoShape* Shape, TGeoMedium* Medium):
  TGeoVolume(Name,Shape,Medium) {
  // can Ralf make more than one pulse out of a waveform ?
  fListOfHits = new TObjArray();
  fListOfHits->SetOwner(kFALSE);

  fListOfDigis = new TObjArray();
  fListOfDigis->SetOwner(kFALSE);

                                        // the counter color will depend on whether the bar has hits
  SetLineColor(kCyan + 1);
  SetFillColor(kCyan + 1);
  SetTransparency(0);

  fListOfHistograms = new TObjArray();
}

//-----------------------------------------------------------------------------
  void TEvdCrvCounter::AddRecoPulse(TCrvRecoPulse* Hit) {
    fListOfHits->Add(Hit);
    
    SetLineColor(kRed + 1);
    SetFillColor(kRed + 1);
  }

//-----------------------------------------------------------------------------
  void TEvdCrvCounter::AddDigi(TCrvDigi* Digi) {
    fListOfDigis->Add(Digi);
    
    SetLineColor(kRed + 1);
    SetFillColor(kRed + 1);
  }

//-----------------------------------------------------------------------------
  void TEvdCrvCounter::PlotWaveforms() {

    TLOG(TLVL_DEBUG) << "START";

    if (fgCanvasWf) delete fgCanvasWf;
    std::string cname = std::format("c_sbid_{:04}",fSbid);
    fgCanvasWf = new TCanvas(cname.data(),cname.data(),1200,800);

    fgCanvasWf->Divide(2,2);

    fListOfHistograms->Delete();

    int first_hist[4] = {1,1,1,1};

    int np = fListOfDigis->GetEntriesFast();
    for (int ip=0; ip<np; ++ip) {
      TCrvDigi* digi = (TCrvDigi*) fListOfDigis->UncheckedAt(ip);
      int sipm = digi->Sipm();
      TLOG(TLVL_DEBUG+1) << std::format("ip:{} sipm{}",ip,sipm);
      std::string hname = std::format("h_{:04}_{}_{:02}",digi->Sbid(),sipm,ip);
      // CRV waveform has 12 samples
      int ns = digi->Ns();
      TH1F* h = new TH1F(hname.data(),hname.data(),ns,0,ns);
      fListOfHistograms->Add(h);
      for (int is=0; is<ns; is++) {
        int adc = digi->Adc()[is];
        h->SetBinContent(is+1,adc);
        h->SetBinError  (is+1,  0);
      }
      TLOG(TLVL_DEBUG+1) << std::format("first_hist[sipm]:",first_hist[sipm]);
      fgCanvasWf->cd(sipm+1);
      if (first_hist[sipm] == 1) {
        h->Draw("");
        first_hist[sipm] = 0;
      }
      else {
        h->Draw("same");
      }
    }

    fgCanvasWf->Modified();
    fgCanvasWf->Update();

    TLOG(TLVL_DEBUG) << "END";
  }

//-----------------------------------------------------------------------------
void TEvdCrvCounter::Clear(Option_t* Opt) {
  fListOfHits->Clear();
  fListOfDigis->Clear();
  SetLineColor(kCyan + 1);
  SetFillColor(kCyan + 1);
  fListOfHistograms->Delete();
}

//-----------------------------------------------------------------------------
void TEvdCrvCounter::Print(Option_t* Opt) const {
  std::cout << std::format("TEvdCrvCounter::Print emoe\n");
}

//-----------------------------------------------------------------------------
  void TEvdCrvCounter::PrintRecoPulses() {
    bool banner_printed(false);

    TIter next(fListOfHits);
    while (auto* pulse = (TCrvRecoPulse*) next()) {
      if (not banner_printed) {
        pulse->Print("banner");
        banner_printed = true;
      }
      pulse->Print("data");
    }
  }

}
