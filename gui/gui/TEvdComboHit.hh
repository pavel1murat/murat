//
#ifndef __murat_gui_TEvdComboHit_hh__
#define __murat_gui_TEvdComboHit_hh__

#include "TMarker.h"
#include "TObject.h"

#include "TObjArray.h"
#include "TPolyLine3D.h"

namespace mu2e {
  class ComboHit;                       // for now, prototyping interfaces
}

class TComboHit;

namespace murat {
//-----------------------------------------------------------------------------
class TEvdComboHit: public TObject {

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

  TPolyLine3D fPolyLine;

  TMarker fXYMarker;
  TMarker fTZMarker;

  TComboHit*   fCh;                     // backward pointer

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
  TEvdComboHit() {}
  TEvdComboHit(TComboHit* Ch, int TrackFlag = 0);
  TEvdComboHit(float X, float Y, float Z, float T, float EDep, int Pdg, float P);

  virtual ~TEvdComboHit() {}
//-----------------------------------------------------------------------------
// overloaded functions of TObject
//-----------------------------------------------------------------------------
  virtual void  Draw         (Option_t* option = "") override;
  virtual void  Paint        (Option_t* option = "") override;
  void          PrintComboHit(Option_t* option = "");               // MENU* 

  ClassDefOverride(murat::TEvdComboHit,0)
};

//-----------------------------------------------------------------------------
struct StationData_t {
public:
  TObjArray fHits;

  StationData_t() { fHits.SetOwner(1); }
  
  TEvdComboHit* GetHit  (int I) { return (TEvdComboHit*) fHits.UncheckedAt(I); }
  int           GetNHits()      { return fHits.GetEntriesFast(); }
};

//-----------------------------------------------------------------------------
struct HitData_t {
   StationData_t   fStation[18];

  HitData_t() { }

};

}
#endif
