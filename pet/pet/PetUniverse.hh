#ifndef Mu2eG4_PetUniverse_hh
#define Mu2eG4_PetUniverse_hh
//
// (Pure virtual) Umbrela for the the Mu2e G4 world classes 
//
// $Id: PetUniverse.hh,v 1.2 2013/06/04 19:20:18 murat Exp $
// $Author: murat $
// $Date: 2013/06/04 19:20:18 $
//
// Original author K. Genser to generalize Mu2eWorld
//
// Notes
//


// C++ includes
#include <vector>


// Mu2e includes
#include "murat/pet/PetG4Helper.hh"
#include "murat/pet/PetGeometryService.hh"

// G4 includes
#include "G4Types.hh"

//G4 forward reference
class G4VPhysicalVolume;

namespace mu2e {

  // Forward references within mu2e namespace.
  class SimpleConfig;

  class PetUniverse {
  public:

    PetUniverse();
    virtual ~PetUniverse();

    // Construct everything.
    // The non-const return type is eventually required 
    // by G4VUserDetectorConstruction::Construct();
    virtual G4VPhysicalVolume * construct() = 0;

  protected:

    // Utility functions.
    static void setUnits( std::vector<double>& V, G4double unit );

    // A helper function for debugging.  Print a subset of the physical volume store
    static void printPhys();

    // geometry service
    PetGeometryService const & _geom;

    // Stash a pointer to the config object so that all methods can get at it easily.
    SimpleConfig const & _config; // make it ref?? (some functions need to change before it...

    // Access to the G4HelperService.
    PetG4Helper * _helper;

    int  _verbosityLevel;

  };

} // end namespace mu2e
#endif /* Mu2eG4_PetUniverse_hh */
