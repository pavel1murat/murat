//

#include "TPad.h"
#include "murat/gui/TComboHitData.hh"
#include "murat/gui/TEvdManager.hh"

ClassImp(TComboHitData)

//-----------------------------------------------------------------------------
TComboHitData::TComboHitData(float X, float Y, float Z, float T, float EDep, int Pdg, float P) { 
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

//-----------------------------------------------------------------------------
void TComboHitData::Paint(Option_t* Option) {
  const char oname[] = "TComboHitData::Paint";

  //  int   iv;

  TEvdManager* vm = TEvdManager::Instance();

  int view = vm->GetCurrentView()->Type();

  if      (view == vm->GetViewID("xy")) PaintXY (Option);
  else if (view == vm->GetViewID("tz")) PaintTZ (Option);
  else {
    printf("[%s] >>> ERROR: unknown view: %i, DO NOTHING\n",oname,view);
  }

  gPad->Modified();
}

//_____________________________________________________________________________
void TComboHitData::PaintXY(Option_t* Option) {
  fXYMarker.Paint(Option);
}

//_____________________________________________________________________________
void TComboHitData::PaintTZ(Option_t* Option) {
  fTZMarker.Paint(Option);
}


//-----------------------------------------------------------------------------
void TComboHitData::Print(Option_t* Option) const {
  //const char oname[] = "TComboHitData::Paint";

  printf("TComboHitData::Print: time = %10.3f Z=%10.3f edep=%10.6f pdg=%10i P=%10.3f\n",
	 fT,fZ,fEDep,fPdg,fP);
}

