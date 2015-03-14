#ifndef PetG4Helper_PetG4Helper_hh
#define PetG4Helper_PetG4Helper_hh
//
// The design of G4 requires that users new many objects and then delete
// them at the appropriate time, usually the end of the G4 run.  This
// Service exists to manage the delete automatically.  It is also available
// as a place to create any other required singleton-like behaviour for
// support of G4.  For technical reasons, this cannot be done by making
// Mu2eG4RunManager a singleton.
//
// $Id: PetG4Helper.hh,v 1.1 2013/06/04 19:20:18 murat Exp $
// $Author: murat $
// $Date: 2013/06/04 19:20:18 $
//
// Original author Rob Kutschke
//

// C++ includes
#include <map>
#include <string>

// Framework include files
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

// Mu2e includes
#include "G4Helper/inc/AntiLeakRegistry.hh"
#include "murat/pet/PetVolumeInfo.hh"

namespace mu2e {

  class PetG4Helper {
  public:
    PetG4Helper(const fhicl::ParameterSet&, art::ActivityRegistry&);
    ~PetG4Helper();

    AntiLeakRegistry& antiLeakRegistry(){ return _antiLeakRegistry; }

    // Versions of the map [] operator that check for errors.
    PetVolumeInfo& locateVolInfo( const std::string key);
    void addVolInfo( const PetVolumeInfo& info );

  private:

    AntiLeakRegistry _antiLeakRegistry;

    std::map<std::string,PetVolumeInfo> _volumeInfoList;

  };

}

DECLARE_ART_SERVICE(mu2e::PetG4Helper, LEGACY)
#endif /* PetG4Helper_PetG4Helper_hh */
