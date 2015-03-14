#ifndef Mu2eG4_PetWorld_hh
#define Mu2eG4_PetWorld_hh
//
// Construct the Mu2e G4 world and serve information about that world.
// Note that the class inherits from Mu2eUniverse now
//
// $Id: PetWorld.hh,v 1.3 2013/11/04 21:09:32 murat Exp $
// $Author: murat $
// $Date: 2013/11/04 21:09:32 $
//
// Original author Rob Kutschke
//
// Notes
// 1) The variable _volumeInfoList, holds some pointers to information
//    about volumes.  It also holds some position information.  The
//    data member centerInWorld is not always meaningful.  If you follow
//    the trail from a given volume back to the world volume, if you
//    do not find any rotations, then this data member is meaningful.
//    If you do find a rotation, then the data member is not meaningful.
//    This data member will only be filled for some of the upper level
//    volumes.   It won't be filled for straws and crystals.  It's purpose
//    is to make it easier to break up one giant method into many smaller
//    ones, by allowing code to look up its mother volume by name.
//    The bottom line is that this is not a fully general facility and
//    must be used with care.

#include <string>
#include <memory>
#include <vector>
#include <map>

// Forward references.
class G4Material;
class G4Mag_UsualEqRhs;
class G4UserLimits;

// Mu2e includes
#include "murat/pet/PetUniverse.hh"
#include "murat/pet/PetVolumeInfo.hh"
#include "Mu2eG4/inc/FieldMgr.hh"
#include "murat/pet/PetG4Helper.hh"
#include "GeomPrimitives/inc/TubsParams.hh"

//G4 includes
#include "G4String.hh"
#include "G4Colour.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4VisAttributes.hh"

namespace mu2e {
					// Forward references within mu2e namespace.
  class SimpleConfig;
  class PetSensitiveDetectorHelper;

  class PetWorld : public PetUniverse {
  public:
					// note no ownership passing

    explicit PetWorld(PetSensitiveDetectorHelper *sdHelper) : sdHelper_(sdHelper) {}

    // Construct everything. The non-const return type is eventually required 
    // by G4VUserDetectorConstruction::Construct();

    virtual G4VPhysicalVolume * construct() override;

  private:

    // Do all of the work.
    G4VPhysicalVolume * constructWorld();

    // Break the big task into many smaller ones.

    void constructStepLimiters();

    void instantiateSensitiveDetectors();

    // Field managers for the different regions of magnetic field.
    // These have a lifetime equal to that of the G4 geometry.

    PetSensitiveDetectorHelper *sdHelper_; // Non-owning
  };
}

#endif
