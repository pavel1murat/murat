#ifndef Mu2eG4_PetnestTubs_hh
#define Mu2eG4_PetnestTubs_hh
//
// Free function to create and place a new G4Tubs, place inside a logical volume.
//
// $Id: PetnestTubs.hh,v 1.1 2013/06/04 19:20:18 murat Exp $
// $Author: murat $
// $Date: 2013/06/04 19:20:18 $
//
// Original author Rob Kutschke
//

#include <string>
#include <vector>

#include "murat/pet/PetVolumeInfo.hh"
#include "GeomPrimitives/inc/TubsParams.hh"

class G4CSGSolid;
class G4LogicalVolume;
class G4Material;
class G4VPhysicalVolume;

// G4 includes
#include "G4Colour.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"


namespace mu2e {

  PetVolumeInfo PetnestTubs ( std::string const& name,
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
                        );



  // Alternate argument list, using a vector for the parameters.
  inline PetVolumeInfo PetnestTubs ( std::string const& name,
                               std::vector<double>&  params,
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
    return PetnestTubs( name,
                     &params[0],
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

  // Alternate argument list, using a TubsParams object for the parameters.
  inline PetVolumeInfo PetnestTubs ( std::string const& name,
                               TubsParams const & params,
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
    return PetnestTubs( name,
                     params.data(),
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
  // using a TubsParams & PetVolumeInfo object for the parameters.
  PetVolumeInfo PetnestTubs ( std::string const& name,
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
                        );

  // Alternate argument list (and different behavior)
  // using a TubsParams & PetVolumeInfo object for the parameters.
  inline PetVolumeInfo PetnestTubs ( std::string const& name,
                               TubsParams const & params,
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
    return PetnestTubs( name,
                     params.data(),
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

#endif /* Mu2eG4_PetnestTubs_hh */
