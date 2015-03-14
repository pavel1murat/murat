//
// Umbrela for the the Mu2e G4 world classes
//
// $Id: PetUniverse.cc,v 1.2 2013/06/04 19:20:18 murat Exp $
// $Author: murat $
// $Date: 2013/06/04 19:20:18 $
//
// Original author Rob Kutschke
//
//

// C++ includes
#include <iostream>
#include <vector>
#include <iomanip>

// Framework includes
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "cetlib/exception.h"

// Mu2e includes
#include "murat/pet/PetG4Helper.hh"
#include "murat/pet/PetUniverse.hh"

// G4 includes
#include "G4PhysicalVolumeStore.hh"

using namespace std;

namespace mu2e {

  PetUniverse::PetUniverse():
    _geom(*(art::ServiceHandle<PetGeometryService>())),
    _config(_geom.config()), 
    _helper(&(*(art::ServiceHandle<PetG4Helper>())))
  {} // beware of the order of initialization/declarations

  PetUniverse::~PetUniverse(){
  }

  // Convert to base units for all of the items in the vector.
  void PetUniverse::setUnits( vector<double>& V, G4double unit ){
    for ( vector<double>::iterator b=V.begin(), e=V.end();
          b!=e; ++b){
      *b *= unit;
    }
  }

  // A helper function for debugging.  Print a subset of the physical volume store.
  void PetUniverse::printPhys() {
    G4PhysicalVolumeStore* pstore = G4PhysicalVolumeStore::GetInstance();
    int n(0);
    for ( std::vector<G4VPhysicalVolume*>::const_iterator i=pstore->begin(); i!=pstore->end(); i++){
      cout << "Physical Volume: "
           << setw(5) << n++
           << (*i)->GetName()
           << endl;
      if ( n > 25 ) break;
    }

  }

} // end namespace mu2e
