//
#ifndef __murat_gui_TComboHitEvd_hh__
#define __murat_gui_TComboHitEvd_hh__

#include "TMarker.h"
#include "TCanvas.h"
#include "TPolyLine3D.h"

namespace mu2e {
  class ComboHit;                       // for now, prototyping interfaces
}

class TComboHit;
class TStrawDigi;

namespace murat {
//-----------------------------------------------------------------------------
class TComboHitEvd: public TPolyLine3D{

public:
  float   fX;
  float   fY;
  float   fZ;
  float   fT;
  int     fPdg;
  int     fPln;
  int     fPnl;
  int     fLay;
  int     fStr;
  float   fEDep;
  float   fP;     // partcle momentum

  TMarker fXYMarker;
  TMarker fTZMarker;

  TComboHit*   fCh    {nullptr};        // backward pointer
  TStrawDigi*  fSd    {nullptr};
  TCanvas*     fCanvas{nullptr};
  TH1F*        fHist  {nullptr};

  mu2e::ComboHit* fOfflineCh;
    
  int Flag() { 
    if (fPdg == 11) {
      if (fP > 20)         return 1;
      else                 return 2;
    }
    else if (fPdg ==  -11) return 3;
    else if (fPdg ==   13) return 4;
    else if (fPdg ==  -13) return 5;
    else if (fPdg == 2212) return 6;
    else                   return 7;
  }

  float  T()    { return fT; }
  float  Time() { return fT; }

  float  X()    { return fX; }
  float  Y()    { return fY; }
  float  Z()    { return fZ; }
  
public:
  TComboHitEvd(TComboHit* Ch, int TrackFlag = 0);
  TComboHitEvd(float X, float Y, float Z, float T, float EDep, int Pdg, float P);

  virtual ~TComboHitEvd() {}

  void SetStrawDigi(TStrawDigi* Sd) { fSd = Sd; }
//-----------------------------------------------------------------------------
// overloaded functions of TObject
//-----------------------------------------------------------------------------
//  virtual void  Draw         (Option_t* option = "") override;
  void          PrintComboHit();                                      // *MENU* 
  void          PlotWaveform ();                                      // *MENU* 

  ClassDefOverride(murat::TComboHitEvd,0)
};

}
#endif
