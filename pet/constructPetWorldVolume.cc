// Free function to create world mother volume and partly fill it with
// dirt around the formal hall box.
//
// $Id: constructPetWorldVolume.cc,v 1.2 2013/11/04 21:09:32 murat Exp $
// $Author: murat $
// $Date: 2013/11/04 21:09:32 $
//
// Original author KLG based on Mu2eWorld constructDirt
// Updated by Andrei Gaponenko.

// Mu2e includes.
#include "murat/pet/constructPetWorldVolume.hh"
#include "murat/pet/PetVolumeInfo.hh"
#include "murat/pet/PetGeomHandle.hh"
#include "murat/pet/PetWorldG4.hh"

#include "murat/pet/PetG4Helper.hh"

#include "Mu2eG4/inc/MaterialFinder.hh"
#include "murat/pet/PetnestBox.hh"
#include "murat/pet/PetfinishNesting.hh"

// G4 includes
#include "G4Material.hh"
#include "G4Color.hh"
#include "G4Box.hh"

using namespace std;

namespace mu2e {

  PetVolumeInfo constructPetWorldVolume(const SimpleConfig &config) {
    // A helper class.
    MaterialFinder materialFinder(config);

    // Dimensions and material of the world.

    G4Material* worldMaterial = materialFinder.get("world.materialName");

    const bool worldBoxVisible     = config.getBool("world.boxVisible");
    const bool worldBoxSolid       = config.getBool("world.boxSolid");

    const bool forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible");
    //    const bool doSurfaceCheck      = config.getBool("g4.doSurfaceCheck");
    const bool placePV             = true;

    PetGeomHandle<PetWorldG4> world;

    PetVolumeInfo worldInfo = PetnestBox("World", world->halfLengths(),
                                   worldMaterial, 0, G4ThreeVector(),
                                   0,
                                   0, worldBoxVisible, G4Colour::Red(), worldBoxSolid,
                                   forceAuxEdgeVisible, placePV, false);
    return worldInfo;

  } // constructPetWorldVolume()

} // namespace mu2e
