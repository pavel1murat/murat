#ifndef murat_tdr_Disk_hh
#define murat_tdr_Disk_hh

#include "TVector2.h"
#include "TMath.h"
#include "TPolyLine.h"
#include "TEllipse.h"
#include "TObjArray.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TH1.h"

class HexIndex {
public:
  int  fL;
  int  fK;
  
  HexIndex() {}
  HexIndex(int L, int K) {
    fL = L;
    fK = K;
  }

  void Set(int L, int K) { fL = L; fK = K; }

  ~HexIndex() {}

};

//-----------------------------------------------------------------------------
class Disk {
public:

  double   fRMin;
  double   fRMax;
  double   fHexSize;
  double   fMinFraction;

  int      fNRings;
  int      fNCrystals;
  int      fNInsideTot;

  int      fMinRing;
  int      fMaxRing;

  int*     fFirst   ;           // index of the first crystal in the ring
  int*     fNCrystalsPerRing ;  // number of crystals in the i-th ring
  int*     fNInside ;
  
  HexIndex  fPos[6];

  static  Disk* gDisk;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  Disk();
  Disk(double RMin, double RMax, double HexSize, double Fraction = 1.);
  ~Disk();

				// for crystal number I return its hex index
  HexIndex  GetHexIndex(int I);
				// for crystal number I return its ring number
  int       GetRing (int I);
  int       GetRing (HexIndex* Index);

  double    GetRMin() { return fRMin; }
  double    GetRMax() { return fRMax; }

  int       GetNRings   () { return fNRings;    }
  int       GetNCrystals() { return fNCrystals; }

  void      GetPosition(int I          , TVector2* Pos);
  void      GetPosition(HexIndex* Index, TVector2* Pos);

  int       GetNCrystalsPerRing(int Ring) { return fNCrystalsPerRing[Ring]; }
  int       GetNInside         (int Ring) { return fNInside         [Ring]; }

  int       GetNInsideTot() { return fNInsideTot;    }

  double    GetActiveArea() { return GetCrystalArea()*GetNInsideTot(); }
  double    GetTotalArea () { return (fRMax*fRMax-fRMin*fRMin)*TMath::Pi() ; }

  double    GetRadius(int I);
  double    GetRadius(HexIndex* Index);

  double    GetCrystalArea() { return  fHexSize*fHexSize*sqrt(3.)/2.; }
  
  int       IsInside(HexIndex* Index, double *Fraction);

  void      Draw(Option_t* Opt = "");

  static void  Test1(double RMin, double RMax, double HexSize, double Fraction = 1., const char* Opt="");

  static void  OptimizeGeometry(double RMin, double RMax, double HexCrystalSize);
};


#endif
