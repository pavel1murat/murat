//
#ifndef __murat_gui_TComboHitData_hh__
#define __murat_gui_TComboHitData_hh__

#include "TMarker.h"
#include "TObject.h"

#include "TObjArray.h"

//-----------------------------------------------------------------------------
class TComboHitData: public TObject {
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
  TComboHitData() {}

  virtual ~TComboHitData() {}

  TComboHitData(float X, float Y, float Z, float T, float EDep, int Pdg, float P);

  virtual void  Paint      (Option_t* option = "");
  virtual void  PaintXY    (Option_t* option = "");
  virtual void  PaintTZ    (Option_t* option = "");

  virtual void  Print      (Option_t* option = "") const ;   // *MENU* 

  ClassDef(TComboHitData,0)
};

//-----------------------------------------------------------------------------
struct StationData_t {
public:
  TObjArray fHits;

  StationData_t() { fHits.SetOwner(1); }
  
  TComboHitData* GetHit  (int I) { return (TComboHitData*) fHits.UncheckedAt(I); }
  int            GetNHits()      { return fHits.GetEntriesFast(); }
};

//-----------------------------------------------------------------------------
struct HitData_t {
   StationData_t   fStation[18];

  HitData_t() { }

};

#endif
