///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef pet_TZubalVolume
#define pet_TZubalVolume

#include "TObject.h"
#include "TObjArray.h"

//-----------------------------------------------------------------------------
class TZubalVolume: public TObject {
public:
  int  fID;       // normally, the volume number
  int fStatus;   // 1: closed   0: updated

  int fIxMin;
  int fIxMax;
  int fIyMin;
  int fIyMax;
  int fIzMin;
  int fIzMax;

  int fColor;
  int fNVoxels;
  
  int fNAddedVoxels;

  TObjArray* fListOfVoxels;		// map of voxels

  short*    fListOfAddedVoxels;
  int       fMaxAddedVoxels;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  TZubalVolume();
  TZubalVolume(int ID, int Ix, int Iy, int Iz, int Color);
  ~TZubalVolume();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  int ID     () { return fID; }
  int Color  () { return fColor; }
  int NVoxels() { return fNVoxels; }
					// 1 if voxel is present
  int HasVoxel(int Ix, int Iy, int Iz);

  short*  ListOfAddedVoxels() { return fListOfAddedVoxels; }
  int     MaxAddedVoxels() { return fMaxAddedVoxels; }
  int     NAddedVoxels  () { return fNAddedVoxels; }

  void    ResetListOfAddedVoxels(int& Size, short*& List);
  

  int     AddVoxel(int Ix, int Iy, int Iz);

  ClassDef(TZubalVolume,0)
};


#endif
