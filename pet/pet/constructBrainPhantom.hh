#ifndef murat_BrainPhantom_hh
#define murat_BrainPhantom_hh
//
// Free function to create the calorimeter.
//
// $Id: constructBrainPhantom.hh,v 1.1 2013/11/04 21:09:32 murat Exp $
// $Author: murat $
// $Date: 2013/11/04 21:09:32 $
//
// Notes:
// 1) Arguments are:
//    1 - pointer to the mother logical volume.
//    2 - z postition of the origin of the Mu2e coordintate system in the
//        frame of the mother.
//    3 - geometry file

// Mu2e includes.
#include "murat/pet/PetVolumeInfo.hh"

namespace mu2e {

  class SimpleConfig;

  PetVolumeInfo constructBrainPhantom(PetVolumeInfo const& mother,
				      SimpleConfig  const& config);

}

#endif /* murat_constructBrainPhantom_hh */
