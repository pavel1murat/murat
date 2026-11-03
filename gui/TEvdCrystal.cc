///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <format>
#include "TCanvas.h"
#include "TH1F.h"

#include "Stntuple/obj/TCaloDigiBlock.hh"

#include "murat/gui/TEvdCrystal.hh"
#include "murat/gui/TEvdManager.hh"

#include "TRACE/tracemf.h"
#define TRACE_NAME "TEvdCrystal"

ClassImp(murat::TEvdCrystal)

namespace murat {

  TCanvas* TEvdCrystal::fgCanvasWf = nullptr;

//-----------------------------------------------------------------------------
  TEvdCrystal::TEvdCrystal(int Index, TGeoShape* Shape, TGeoMedium* Medium):
    TGeoVolume(Form("crystal_%04i",Index),Shape,Medium) {
    fCid              = Index;
    fListOfHits       = new TObjArray();
    fEDep             = 0;
    fListOfHistograms = new TObjArray();
  }
  
//-----------------------------------------------------------------------------
  void TEvdCrystal::AddHit(TCaloHit* Hit) {
    fListOfHits->Add(Hit);
    fEDep += Hit->EDep();

    SetLineColor(kRed + 1);
    SetFillColor(kRed + 1);
  }

//-----------------------------------------------------------------------------
void TEvdCrystal::Clear(Option_t* Opt) {
  fListOfHits->Clear();
  fEDep = 0;
                                        // crystal color will depend on whether the crystal has hits
  SetLineColor(kOrange + 1);
  SetFillColor(kOrange + 1);
  SetTransparency(0);
}

//-----------------------------------------------------------------------------
  void TEvdCrystal::PlotWaveforms() {

    TLOG(TLVL_DEBUG) << "START";

    auto* vm = TEvdManager::Instance();

    TCaloDigiBlock* calo_digi_block = (TCaloDigiBlock*) vm->FindDataBlock("CaloDigiBlock");

    if (fgCanvasWf) delete fgCanvasWf;
    std::string cname = std::format("c_calo_cid_{:04}",fCid);
    fgCanvasWf = new TCanvas(cname.data(),cname.data(),1200,600);

    fgCanvasWf->Divide(2,1);

    fListOfHistograms->Delete();

    int first_hist[2] = {0,0};

    int nh = fListOfHits->GetEntriesFast();
    for (int ih=0; ih<nh; ++ih) {
      TCaloHit* hit = Hit(ih);
      for (int sipm=0; sipm<2; sipm++) {
        int icrd = hit->CrdIndex(sipm);
        if (icrd < 0)                    continue ;
        TLOG(TLVL_DEBUG+1) << std::format("ih:{} sipm{}",ih,sipm);
        
        TCaloDigi* digi = calo_digi_block->Digi(icrd);
        
        int chid = 2*fCid+sipm;
        
        std::string hname = std::format("h_{:04}_{:02}",chid,ih);
        // calo waveform has a variable number of samples, book histogram with 100 bins
        int ns = digi->Ns();
        TH1F* h = new TH1F(hname.data(),hname.data(),100,0,100);
        h->SetMinimum(2000.);
        h->SetMaximum(4100.);
        fListOfHistograms->Add(h);
        for (int is=0; is<ns; is++) {
          int adc = digi->Wf()[is];
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
    }

    fgCanvasWf->Modified();
    fgCanvasWf->Update();
    
    TLOG(TLVL_DEBUG) << "END";
  }

//-----------------------------------------------------------------------------
// there are two SiPMs
//-----------------------------------------------------------------------------
void TEvdCrystal::PrintHits() {

  int nhits = NHits();
  std::cout << std::format("N(hits):{}\n",nhits);

  if (nhits > 0) {
    bool header_printed(false);
    for (int i=0; i<nhits; i++) {
      TCaloHit* hit = Hit(i);
      if (not header_printed) {
        hit->Print("banner");
        header_printed = true;
      }
      hit->Print("data");
    }
  }
}

}
