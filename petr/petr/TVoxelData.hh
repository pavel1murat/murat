///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef murat_pet_TVoxelData
#define murat_pet_TVoxelData

#include "TObject.h"
#include <vector>

#include "PetRecoConstants.hh"

class TVoxelData: public TObject {
public:
  
  struct Occ_t {
    int   fIPhi;
    int   fIz;
    int   fNEntries;
  };
				// voxel itself
  int         fIx;
  int         fIy;
  int         fIz;
  
  int         fNEvents;		// N(events) generated for this voxel
  int         fNin    ;		// total within acceptance
  
  int         fMatrix[IMAGER_PHI_SEG][IMAGER_Z_SEG];	// store population of the cells
  int         fNCells[IMAGER_PHI_SEG][IMAGER_Z_SEG];	// n(cells) per [iphi][iz]

				// data, first ncells for [0][0] ...

  int         fNCellsTot;	// total number of filled cells
				// 
  std::vector<Occ_t>* fCellData[IMAGER_PHI_SEG][IMAGER_Z_SEG];  // ! 3 or 4 per[iphi1][iz1]
//-----------------------------------------------------------------------------
//  constructors and destructor
//-----------------------------------------------------------------------------
  TVoxelData();
  TVoxelData(int Ix, int Iy, int Iz);
  ~TVoxelData();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  int                  GetNCells  (int IPhi, int IZ) { return fNCells  [IPhi][IZ]; }
  std::vector<Occ_t>*  GetCellData(int IPhi, int IZ) { return fCellData[IPhi][IZ]; }
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  virtual void Copy (TVoxelData* Voxel ) const;
  virtual void Clear(Option_t*   Opt="") ;
  virtual void Print(Option_t*   Opt="") const ;

  ClassDef(TVoxelData,1)
};

#endif
