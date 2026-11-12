//
#ifndef __murat_gui_TComboHitEvd_hh__
#define __murat_gui_TComboHitEvd_hh__

#include "TMarker.h"
#include "TCanvas.h"
#include "TGeoVolume.h"

namespace mu2e {
  class ComboHit;                       // for now, prototyping interfaces
}

class TComboHit;
class TStrawDigi;

namespace murat {
//-----------------------------------------------------------------------------
class TComboHitEvd: public TGeoVolume {

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
  // tube along the wire of the hit: radius range [RMin,RMax], half-length the resolution along the wire
  TComboHitEvd(TComboHit* Ch, double RMin, double RMax, TGeoMedium* Medium);

  virtual ~TComboHitEvd() {}

  void SetStrawDigi(TStrawDigi* Sd) { fSd = Sd; }
//-----------------------------------------------------------------------------
// overloaded functions of TObject
//-----------------------------------------------------------------------------
//  virtual void  Draw         (Option_t* option = "") override;
  void          PrintComboHit();                                      // *MENU* 
  void          PlotWaveform ();                                      // *MENU* 
  void          ShowPanel    ();                                      // *MENU*

  ClassDefOverride(murat::TComboHitEvd,0)
};

}
#endif
