#ifndef Mu2eG4_finishNesting_hh
#define Mu2eG4_finishNesting_hh
//
// Free function to be used by the nest... functions
//
// $Id: PetfinishNesting.hh,v 1.1 2013/06/04 19:20:18 murat Exp $
// $Author: murat $
// $Date: 2013/06/04 19:20:18 $
//
// Original author KLG
//


//class G4RotationMatrix;
//class G4ThreeVector;
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

class PetVolumeInfo;

class G4Material;
class G4LogicalVolume;
class G4Colour;

namespace mu2e {

  void PetfinishNesting(PetVolumeInfo& info,
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
                     bool const doSurfaceCheck,
                     bool const verbose = false
                     );

}

#endif /* Mu2eG4_finishNesting_hh */
