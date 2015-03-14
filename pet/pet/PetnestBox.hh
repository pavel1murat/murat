#ifndef Mu2eG4_PetnestBox_hh
#define Mu2eG4_PetnestBox_hh
//
// Free function to create a new G4 Box, placed inside a logical volume.
//
// $Id: PetnestBox.hh,v 1.1 2013/06/04 19:20:18 murat Exp $
// $Author: murat $
// $Date: 2013/06/04 19:20:18 $
//
// Original author Rob Kutschke
//

#include <string>
#include <vector>

#include "murat/pet/PetVolumeInfo.hh"

class G4Material;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4CSGSolid;

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Colour.hh"

namespace mu2e {

  PetVolumeInfo PetnestBox ( std::string const& name,
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
                       );

  // Alternate argument list, using a vector for the half dimensions.
  inline PetVolumeInfo PetnestBox ( std::string const& name,
                              std::vector<double> const&  halfDim,
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
    return PetnestBox( name,
                    &halfDim[0],
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
  }

  // Alternate argument list (and different behavior)
  // using PetVolumeInfo object
  PetVolumeInfo PetnestBox ( std::string const& name,
                       double const halfDim[3],
                       G4Material* material,
                       G4RotationMatrix const* rot,
                       G4ThreeVector const& offset,
                       PetVolumeInfo const & parent,
                       int copyNo,
                       bool const isVisible,
                       G4Colour const color,
                       bool const forceSolid,
                       bool const forceAuxEdgeVisible,
                       bool const placePV,
                       bool const doSurfaceCheck
                       );

  // Alternate argument list, (and different behavior)
  // using PetVolumeInfo object and using a vector for the half dimensions.
  inline PetVolumeInfo PetnestBox ( std::string const& name,
                              std::vector<double> const&  halfDim,
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
    return PetnestBox( name,
                    &halfDim[0],
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
  }

}

#endif /* Mu2eG4_PetnestBox_hh */
