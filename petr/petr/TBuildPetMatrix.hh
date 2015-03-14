///////////////////////////////////////////////////////////////////////////////
// a toy model: 
//
// 1. the imager: R ~ 8 cm, H = 6cm
//    Nphi = 100, nz = 12 ("crystal size" approximately 5mm x 5mm)
//
// 2. the phantom: ZMax= 5cm, R=5cm, voxel size 5x5x5mm
//    handling of the edge voxels may need to be improved
// 
///////////////////////////////////////////////////////////////////////////////
#include "TObject.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH2F.h"

#include "PetRecoConstants.hh"

#include "TVoxelData.hh"

#include <vector>

class TBuildPetMatrix: public TObject {
public:

  double              fZImager;		// imager half-Z
  double              fRImager;		// imager radius

  int                 fNCrystalsPhi;	// granularity in Phi
  int                 fNCrystalsZ;	// granularity in Z

  int                 fNEventsPerVoxel; 
  int                 fWriteVoxelData; // flag, 0 or 1

  double              fCrystalDz;
  double              fCrystalDPhi;
//-----------------------------------------------------------------------------
// description of the phantom
//-----------------------------------------------------------------------------
  double              fPhantomR;
  double              fPhantomMaxZ;
					// phantom voxelization
  double              fVoxelDx;
  double              fVoxelDy;
  double              fVoxelDz;

  int                 fNVoxelsX;
  int                 fNVoxelsY;
  int                 fNVoxelsZ;
					// phantom voxel map

  TVoxelData*         fVoxel[kNVoxelsX][kNVoxelsY][kNVoxelsZ];

  TVoxelData*         fVoxelBuffer;
//-----------------------------------------------------------------------------
// then at a certain point I'll need the image data - just a 2D array with the 
// number of counts
//-----------------------------------------------------------------------------
  int                 fImageData[IMAGER_PHI_SEG][IMAGER_Z_SEG];
//-----------------------------------------------------------------------------
// reconstructed map - counts per voxel (doubles)
//-----------------------------------------------------------------------------
  double              fRecoMap[kNVoxelsX][kNVoxelsY][kNVoxelsZ];
//-----------------------------------------------------------------------------
// the rest - I/O, histogramming
//-----------------------------------------------------------------------------
  TRandom3*           fRn;
  TH2F*               fHist;
					// diagnostics
  struct Hist_t {
    TH1F* fX;
    TH1F* fY;
    TH1F* fZ;
    TH1F* fTheta;
    TH1F* fPhi1;
    TH2F* fIp2VsIp1;
    TH2F* fIz2VsIz1;
  };

  Hist_t fDiagHist[2] ;
  int    fDebugFlag[100];

  TFile*              fFile;
  TTree*              fTree;
  TBranch*            fVoxelDataBranch;
//-----------------------------------------------------------------------------
//  constructors and destructor
//-----------------------------------------------------------------------------
   TBuildPetMatrix();
  ~TBuildPetMatrix();

//-----------------------------------------------------------------------------
//  accessors
//-----------------------------------------------------------------------------
  TVoxelData*  GetVoxel(int Ix, int Iy, int Iz) { return fVoxel[Ix][Iy][Iz]; }
//-----------------------------------------------------------------------------
//  setters
//-----------------------------------------------------------------------------
  void SetWriteVoxelData(int Flag)     { fWriteVoxelData = Flag; }
  void SetDebugFlag(int I, int IValue) { fDebugFlag[I] = IValue; }
//-----------------------------------------------------------------------------
//  other methods
//-----------------------------------------------------------------------------
  int  BookHist(Hist_t* Hist, int ISet);

  int  Generate(int NEventsPerVoxel);

  int  GenerateVoxelData(int Ix, int Iy, int Iz, int NEvents);

  int  PlotResults(int VoxelIx, int VoxelIy, int VoxelIz, int IPhi, int Iz);

  int  PlotIPhiIZ (int VoxelIX, int VoxelIY, int VoxelIZ, int IPhi = -1, int IZ = -1);

  int  ReadVoxelMap  (const char* Filename);

  void PrintVoxel  (int IX, int IY, int IZ);

  int  ResetDiagHist (Hist_t* Hist);

  void Validate      (int Ix = -1, int Iy = -1, int Iz = -1, int CPhi = -1, int Cz = -1);

  int  WriteVoxelData   (TVoxelData* Voxel);

  ClassDef(TBuildPetMatrix,0)

};

