///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef pet_TZubalPhantom
#define pet_TZubalPhantom

#include "TObject.h"
#include "TObjArray.h"

#include "TZubalVolume.hh"

//-----------------------------------------------------------------------------
class TZubalPhantom: public TObject {
public:
  int fIXMin;
  int fIXMax;
  int fIYMin;
  int fIYMax;
  int fIZMin;
  int fIZMax;

  short fVolID[128][256][256];
  char  fData [128][256][256];

  int  fNVolumes;

  TZubalVolume* fVol;

  TObjArray* fListOfVolumes;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  TZubalPhantom(const char* Filename = NULL);
  ~TZubalPhantom();

  TZubalVolume*  Volume(int I) { 	// 
    return (TZubalVolume*) fListOfVolumes->At(I); 
  }

  char&  Data (int Ix, int Iy, int Iz) { return fData [Iz][Ix][Iy]; }
  short& VolID(int Ix, int Iy, int Iz) { return fVolID[Iz][Ix][Iy]; }

  int NVolumes() { return fNVolumes; }

  int CheckSurroundingVoxels(int Ix, int Iy, int Iz, int Ic, TZubalVolume* Vol);

  int ExpandVolume(TZubalVolume* Vol);
  int Parse3DMap(int IZMin=0, int IZMax=128);

  static void Test01(const char* Filename, int IZ);
  static void Test02(int IZMin, int IZMax, int IVol);

  ClassDef(TZubalPhantom,0)
};


#endif
