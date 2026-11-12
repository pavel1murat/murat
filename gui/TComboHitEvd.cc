//

#include <iostream>
#include <format>

#include "TH1F.h"
#include "TCanvas.h"
#include "TGeoMedium.h"
#include "TGeoTube.h"

#include "Stntuple/obj/TStrawDigi.hh"
#include "Stntuple/obj/TComboHit.hh"
#include "murat/gui/TComboHitEvd.hh"
#include "murat/gui/TEvdManager.hh"

#include "TRACE/tracemf.h"
#define TRACE_NAME "TComboHitEvd"

ClassImp(murat::TComboHitEvd)

namespace {
//-----------------------------------------------------------------------------
// unique name of a combo hit volume
//-----------------------------------------------------------------------------
  std::string comboHitVolumeName() {
    static int n(0);
    return std::format("combohit_{}",n++);
  }
}

namespace murat {
//-----------------------------------------------------------------------------
  TComboHitEvd::TComboHitEvd(TComboHit* Ch, double RMin, double RMax, TGeoMedium* Medium):
  TGeoVolume(comboHitVolumeName().c_str(),new TGeoTube(RMin,RMax,Ch->Ures()),Medium),
  fCh(Ch)
{
  SetLineColor(kRed);
  SetFillColor(kRed);
  SetTransparency(0);
}

//-----------------------------------------------------------------------------
void TComboHitEvd::ShowPanel() {
  if (fCh == nullptr) return;
  TEvdManager::Instance()->OpenPanelView(fCh->Sid().getPlane(), fCh->Sid().getPanel());
}

//-----------------------------------------------------------------------------
void TComboHitEvd::PrintComboHit() {
  //const char oname[] = "TComboHitEvd::Paint";

  fCh->Print("banner+data");
}

//-----------------------------------------------------------------------------
  void TComboHitEvd::PlotWaveform() {
  //const char oname[] = "TComboHitEvd::Paint";

    TLOG(TLVL_DEBUG) << "START";

    if (fSd == nullptr) {
      TLOG(TLVL_ERROR) << std::format("fSd = nullptr, BAIL OUT");
      return ;
    }

    std::string cname = std::format("c_wf_{:02}_{:1}",fCh->Sid().plane(),fCh->Sid().panel());
    if (fCanvas) delete fCanvas;
    fCanvas = new TCanvas(cname.data(),cname.data(),800,600);

    // not limited by time, can go slowly
    if (fHist != nullptr) delete fHist;

    std::string hname = std::format("h_wf_{:02}_{:1}",fCh->Sid().plane(),fCh->Sid().panel());
    int ns = fSd->Ns();
    fHist  = new TH1F(hname.data(),hname.data(),ns,0,ns);

    for (int is=0; is<ns; is++) {
      int adc = fSd->Adc()[is];
      fHist->SetBinContent(is+1,adc);
      fHist->SetBinError  (is+1,  0);
    }

    fHist->Draw("");
    
    fCanvas->Modified();
    fCanvas->Update();

    TLOG(TLVL_DEBUG) << "END";
  }

}
