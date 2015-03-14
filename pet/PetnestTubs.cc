//
// Free function to create and place a new G4Tubs, place inside a logical volume.
//
// $Id: PetnestTubs.cc,v 1.1 2013/06/04 19:20:18 murat Exp $
// $Author: murat $
// $Date: 2013/06/04 19:20:18 $
//
// Original author Rob Kutschke
//

#include <string>

// Mu2e includes
#include "murat/pet/PetnestTubs.hh"
#include "murat/pet/PetfinishNesting.hh"

// G4 includes
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4ThreeVector.hh"

using namespace std;

namespace mu2e {

  //
  // Create and place a G4Tubs inside a logical volume.
  //
  PetVolumeInfo PetnestTubs ( string const & name,
                        double const params[5],
                        G4Material* material,
                        G4RotationMatrix const* rot,
                        G4ThreeVector const & offset,
                        G4LogicalVolume* parent,
                        int copyNo,
                        bool const isVisible,
                        G4Colour const color,
                        bool const forceSolid,
                        bool const forceAuxEdgeVisible,
                        bool const placePV,
                        bool const doSurfaceCheck
                        ){


    PetVolumeInfo info;

    info.name     = name;

    info.solid    = new G4Tubs( name, params[0], params[1], params[2], params[3], params[4]  );

    PetfinishNesting(info,
                  material,
                  rot,
                  offset,
                  parent,
                  copyNo,
                  isVisible,
                  color,
                  forceSolid,
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck
                  );

    return info;

  }

  PetVolumeInfo PetnestTubs ( string const & name,
                        double const params[5],
                        G4Material* material,
                        G4RotationMatrix const* rot,
                        G4ThreeVector const & offset,
                        PetVolumeInfo const & parent,
                        int copyNo,
                        bool const isVisible,
                        G4Colour const color,
                        bool const forceSolid,
                        bool const forceAuxEdgeVisible,
                        bool const placePV,
                        bool const doSurfaceCheck
                        ){


    PetVolumeInfo info(name,offset,parent.centerInWorld);

    info.solid    = new G4Tubs( name, params[0], params[1], params[2], params[3], params[4] );

    PetfinishNesting(info,
                  material,
                  rot,
                  offset,
                  parent.logical,
                  copyNo,
                  isVisible,
                  color,
                  forceSolid,
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck
                  );

    return info;

  }

}
