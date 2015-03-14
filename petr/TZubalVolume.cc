///////////////////////////////////////////////////////////////////////////////
#include <cstdlib>
#include "murat/petr/TZubalVolume.hh"

#include "TArrayC.h"

ClassImp(TZubalVolume)

//-----------------------------------------------------------------------------
TZubalVolume::TZubalVolume(): TObject() { 
  fListOfVoxels      = NULL;
  fNVoxels           = 0;
  fNAddedVoxels      = 0;
  fListOfAddedVoxels = 0;
}

//-----------------------------------------------------------------------------
TZubalVolume::TZubalVolume(int ID, int Ix, int Iy, int Iz, int Color): TObject() {

  TObjArray* slice;
  TArrayC*   line;

  fID    = ID;

  fIxMin = Ix;
  fIxMax = Ix;
  fIyMin = Iy;
  fIyMax = Iy;
  fIzMin = Iz;
  fIzMax = Iz;

  fColor = Color;
					// this is a list of slice data
  fListOfVoxels = new TObjArray(128);

  // for (int i=0; i<128; i++) {
  //   printf(" TZubalVolume::TZubalVolume : slice %3i = %08x\n",i,fListOfVoxels->UncheckedAt(i));
  // }



  slice         = new TObjArray(256);
  line          = new TArrayC(256);
  for (int i=0; i<256; i++) {
    (*line)[i] = 0;
  }
  
  fListOfVoxels->AddAt((TObject*) slice,Iz);
  slice->AddAt((TObject*) line,Ix);
  (*line)[Iy] = 1;

  fNVoxels = 1;

  fMaxAddedVoxels    = 1000;
  fListOfAddedVoxels = new short[3*fMaxAddedVoxels];

  fListOfAddedVoxels[0] = Ix;
  fListOfAddedVoxels[1] = Iy;
  fListOfAddedVoxels[2] = Iz;

  fNAddedVoxels = 1;
}

//-----------------------------------------------------------------------------
TZubalVolume::~TZubalVolume() { 
  TObjArray* a; 

  if (fListOfVoxels != NULL) {
    for (int i=0; i<128; i++) {
      a = (TObjArray*) fListOfVoxels->At(i);
      a->Delete();
    }

    delete fListOfVoxels;
  }
}



//-----------------------------------------------------------------------------
int TZubalVolume::AddVoxel(int Ix, int Iy, int Iz) {

  static int nerrors = 0;

  TObjArray* slice;
  TArrayC*   line;

  if (fIxMin > Ix) fIxMin = Ix;
  if (fIxMax < Ix) fIxMax = Ix;

  if (fIyMin > Iy) fIyMin = Iy;
  if (fIyMax < Iy) fIyMax = Iy;

  if (fIzMin > Iz) fIzMin = Iz;
  if (fIzMax < Iz) fIzMax = Iz;
					// this is a list of slice data

  slice = (TObjArray*) fListOfVoxels->UncheckedAt(Iz);
  if (slice == 0) {
    slice = new TObjArray(256);
    fListOfVoxels->AddAt((TObject*) slice,Iz);
  }

  line = (TArrayC*) slice->UncheckedAt(Ix);
  if (line == 0) {
    line = new TArrayC(256);
    for (int i=0; i<256; i++) {
      (*line)[i] = 0;
    }
    slice->AddAt((TObject*) line,Ix);
  }

  if ((*line)[Iy] != 0) {
    printf("TZubalVolume::AddVoxel: VolID=%3i ERROR : re-adding voxel IX=%3i IY=%3i IZ=%3i, flag=%i\n",
	   fID,Ix,Iy,Iz,(*line)[Iy]);
    nerrors += 1;
  }

  (*line)[Iy] = 1;

  fNVoxels++;

  int loc = fNAddedVoxels*3;
  if (loc == fMaxAddedVoxels*3) {
    short* a = new short[3*fMaxAddedVoxels*2];
    memcpy(a,fListOfAddedVoxels,fMaxAddedVoxels*3*sizeof(short));
    delete fListOfAddedVoxels;
    fListOfAddedVoxels = a;
    fMaxAddedVoxels    = fMaxAddedVoxels*2;
  }

  fListOfAddedVoxels[loc  ] = Ix;
  fListOfAddedVoxels[loc+1] = Iy;
  fListOfAddedVoxels[loc+2] = Iz;

  fNAddedVoxels += 1;

  //  printf("TZubalVolume::AddVoxel: ID=%i nvoxels=%8i\n",fID,fNVoxels);

  return -nerrors;
}

//-----------------------------------------------------------------------------
int TZubalVolume::HasVoxel(int Ix, int Iy, int Iz) {

  int color = 0;

  TObjArray* slice;
  TArrayC*   line;


  slice = (TObjArray*) fListOfVoxels->UncheckedAt(Iz);

  if (slice != 0) {
    line = (TArrayC*) slice->UncheckedAt(Ix);

    if (line != 0) {
      color = (*line)[Iy];
    }
  }

  return color;
}


//-----------------------------------------------------------------------------
// assumes that previous list has been cached somewhere else
//-----------------------------------------------------------------------------
void TZubalVolume::ResetListOfAddedVoxels(int& N, short*& List) {
					// save previous list
  List               = fListOfAddedVoxels;
  N                  = fNAddedVoxels;
					// reallocate empty list
  fMaxAddedVoxels    = fNAddedVoxels;
  fListOfAddedVoxels = new short[3*fMaxAddedVoxels];
  fNAddedVoxels      = 0;
}
