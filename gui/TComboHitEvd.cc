//

#include <iostream>
#include <format>

#include "TH1F.h"
#include "TCanvas.h"

#include "Stntuple/obj/TStrawDigi.hh"
#include "Stntuple/obj/TComboHit.hh"
#include "murat/gui/TComboHitEvd.hh"
#include "murat/gui/TEvdManager.hh"

#include "TRACE/tracemf.h"
#define TRACE_NAME "TComboHitEvd"

ClassImp(murat::TComboHitEvd)

namespace murat {
//-----------------------------------------------------------------------------
  TComboHitEvd::TComboHitEvd(TComboHit* Ch, int TrackFlag) : TPolyLine3D(2) { 
  float tracker_pos[3] = {-3904.,    0.0, 22935.0 };

  double dx = Ch->Ures()*Ch->Ux();
  double dy = Ch->Ures()*Ch->Uy();

  int line_width(2);
  
  if (TrackFlag == 1) {
                                        // small and thick error bars
    dx         = 5.*Ch->Ux();
    dy         = 5.*Ch->Uy();
    line_width = 6;
  }
  
  SetPoint(0, Ch->X()+tracker_pos[0]-dx, Ch->Y()+tracker_pos[1]-dy, Ch->Z()+tracker_pos[2]);
  SetPoint(1, Ch->X()+tracker_pos[0]+dx, Ch->Y()+tracker_pos[1]+dy, Ch->Z()+tracker_pos[2]);
  SetLineWidth(line_width);

  fCh = Ch;
}

//-----------------------------------------------------------------------------
TComboHitEvd::TComboHitEvd(float X, float Y, float Z, float T, float EDep, int Pdg, float P)
  : TPolyLine3D(2) 
{ 


  fX    = X;
  fY    = Y;
  fZ    = Z;
  fT    = T;
  fEDep = EDep;
  fPdg  = Pdg;
  fP    = P;
  
  int   style(0), color(0);
  float size(0.);
  
  if      (Pdg == 11) {
    if    (P    > 20  ) { style = 20; size = 0.8; color = kRed; }
    else                { style = 20; size = 0.4; color = kRed+2; }
  }
  else if (Pdg ==  -11) { style = 24; size = 0.8; color = kBlue;   } 
  else if (Pdg ==   13) { style = 20; size = 0.8; color = kGreen+2;} 
  else if (Pdg ==  -13) { style = 20; size = 0.8; color = kGreen-2;} 
  else if (Pdg == 2212) { style = 20; size = 1.0; color = kBlue+2; } 
  else                  { style = 20; size = 1.0; color = kMagenta;} 
  
  fXYMarker.SetX(fX);
  fXYMarker.SetY(fY);
  fXYMarker.SetMarkerStyle(style);
  fXYMarker.SetMarkerSize (size );
  fXYMarker.SetMarkerColor(color);
  
  fTZMarker.SetX(fZ);
  fTZMarker.SetY(fT);
  fTZMarker.SetMarkerStyle(style);
  fTZMarker.SetMarkerSize (size );
  fTZMarker.SetMarkerColor(color);
}

// //-----------------------------------------------------------------------------
// void TComboHitEvd::Paint(Option_t* Option) {
//   // const char oname[] = "TComboHitEvd::Paint";

//   // //  int   iv;

//   // auto* vm = murat::TEvdManager::Instance();

//   // int view = vm->GetCurrentView()->Type();

//   // if      (view == vm->GetViewID("xy")) PaintXY (Option);
//   // else if (view == vm->GetViewID("tz")) PaintTZ (Option);
//   // else {
//   //   printf("[%s] >>> ERROR: unknown view: %i, DO NOTHING\n",oname,view);
//   // }

//   gPad->Modified();
// }

// //-----------------------------------------------------------------------------
// void TComboHitEvd::Draw(Option_t* Option) {
//   AppendPad(Option);
// }

// //-----------------------------------------------------------------------------
// void TComboHitEvd::Paint(Option_t* Option) {
//   //const char oname[] = "TComboHitEvd::Paint";

//   TPolyLine3D::Paint(Option);
// }

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
