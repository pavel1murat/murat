//
// Free function to create a new G4 Box, placed inside a logical volume.
//
// $Id: PetnestBox.cc,v 1.1 2013/06/04 19:20:18 murat Exp $
// $Author: murat $
// $Date: 2013/06/04 19:20:18 $
//
// Original author Rob Kutschke
//

#include <string>

// Mu2e includes
#include "murat/pet/PetnestBox.hh"
#include "murat/pet/PetfinishNesting.hh"

// G4 includes
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"

using namespace std;

namespace mu2e {

  //
  // Place a box inside a logical volume.
  //
  PetVolumeInfo PetnestBox ( string const& name,
                       double const halfDim[3],
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

    info.name    = name;

    info.solid   = new G4Box( name, halfDim[0], halfDim[1], halfDim[2] );

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

  PetVolumeInfo PetnestBox ( string const& name,
                       double const halfDim[3],
                       G4Material* material,
                       G4RotationMatrix const* rot,
                       G4ThreeVector const& offset,
                       const PetVolumeInfo& parent,
                       int copyNo,
                       bool const isVisible,
                       G4Colour const color,
                       bool const forceSolid,
                       bool const forceAuxEdgeVisible,
                       bool const placePV,
                       bool const doSurfaceCheck
                       ){

    PetVolumeInfo info(name,offset,parent.centerInWorld);

    info.solid   = new G4Box( name, halfDim[0], halfDim[1], halfDim[2] );

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
