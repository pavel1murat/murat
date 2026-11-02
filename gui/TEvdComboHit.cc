//

#include "TPad.h"
#include "Stntuple/obj/TComboHit.hh"
#include "murat/gui/TEvdComboHit.hh"
#include "murat/gui/TEvdManager.hh"

ClassImp(murat::TEvdComboHit)

namespace murat {
//-----------------------------------------------------------------------------
  TEvdComboHit::TEvdComboHit(TComboHit* Ch, int TrackFlag) : TObject(),
                                                             fPolyLine(2) { 
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
  
  fPolyLine.SetPoint(0, Ch->X()+tracker_pos[0]-dx, Ch->Y()+tracker_pos[1]-dy, Ch->Z()+tracker_pos[2]);
  fPolyLine.SetPoint(1, Ch->X()+tracker_pos[0]+dx, Ch->Y()+tracker_pos[1]+dy, Ch->Z()+tracker_pos[2]);
  fPolyLine.SetLineWidth(line_width);

  fCh = Ch;
}

//-----------------------------------------------------------------------------
TEvdComboHit::TEvdComboHit(float X, float Y, float Z, float T, float EDep, int Pdg, float P) { 


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
// void TEvdComboHit::Paint(Option_t* Option) {
//   // const char oname[] = "TEvdComboHit::Paint";

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

//-----------------------------------------------------------------------------
void TEvdComboHit::Draw(Option_t* Option) {
  AppendPad(Option);
}

//-----------------------------------------------------------------------------
void TEvdComboHit::Paint(Option_t* Option) {
  //const char oname[] = "TEvdComboHit::Paint";

  fPolyLine.Paint(Option);
}

//-----------------------------------------------------------------------------
void TEvdComboHit::PrintComboHit(Option_t* Option) {
  //const char oname[] = "TEvdComboHit::Paint";

  fCh->Print("banner+data");
}

}
