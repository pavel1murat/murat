#ifndef murat_pet_constructBrainImager_hh
#define murat_pet_constructBrainImager_hh
//
// Free function to create the calorimeter.
//
// $Id: constructBrainImager.hh,v 1.2 2013/06/04 19:20:18 murat Exp $
// $Author: murat $
// $Date: 2013/06/04 19:20:18 $
//
// Original author Rob Kutschke
//
// Notes:
// 1) Arguments are:
//    1 - pointer to the mother logical volume.
//    2 - geometry file

// Mu2e includes.
#include "murat/pet/PetVolumeInfo.hh"

namespace mu2e {

  class SimpleConfig;

  PetVolumeInfo constructBrainImager( PetVolumeInfo const&  mother,SimpleConfig const& config );

}

#endif /* murat_pet_constructBrainImager_hh */
