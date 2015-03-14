///////////////////////////////////////////////////////////////////////////////
// a toy model: 
//
// 1. the imager: R ~ 8 cm, H = 6cm
//    Nphi = 100, nz = 12 ("crystal size" approximately 5mm x 5mm)
//
// 2. the phantom: ZMax= 5cm, R=5cm, 
// 3. voxel size : 5x5x5mm
//
// 4. handling of the edge voxels may need to be improved
// 
///////////////////////////////////////////////////////////////////////////////

#ifndef murat_pet_PetRecoConstants
#define murat_pet_PetRecoConstants

double const IMAGER_R       =  80.;  // mm
double const IMAGER_ZMAX    =  60;   // mm

int    const IMAGER_PHI_SEG = 100;   // 2*pi*80/100 ~ 5mm
int    const IMAGER_Z_SEG   =  24;   // 5mm in Z

double const PHANTOM_ZMAX   = 50.;   // mm
double const PHANTOM_R      = 50.;   // mm

double const VOXEL_DX       = 5.;    // mm
double const VOXEL_DY       = 5.;    // mm
double const VOXEL_DZ       = 5.;    // mm

int const kNVoxelsX         = 2*PHANTOM_R/VOXEL_DX;
int const kNVoxelsY         = 2*PHANTOM_R/VOXEL_DY;
int const kNVoxelsZ         = 2*PHANTOM_ZMAX/VOXEL_DZ;
     
#endif
