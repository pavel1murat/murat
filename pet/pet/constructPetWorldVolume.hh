#ifndef Mu2eG4_constructPetWorldVolume_hh
#define Mu2eG4_constructPetWorldVolume_hh
//
// Free function to construct World Mother Volume
//
// $Id: constructPetWorldVolume.hh,v 1.1 2013/06/04 19:20:18 murat Exp $
// $Author: murat $
// $Date: 2013/06/04 19:20:18 $
//
// Original author KLG
//

// Mu2e includes.
#include "murat/pet/PetVolumeInfo.hh"

namespace mu2e {

  class SimpleConfig;
  PetVolumeInfo constructPetWorldVolume(const SimpleConfig& config);

}

#endif /* Mu2eG4_constructWorldVolume_hh */
