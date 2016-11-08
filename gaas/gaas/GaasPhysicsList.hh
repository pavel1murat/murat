#ifndef murat_GaasPhysicsList_hh
#define murat_GaasPhysicsList_hh
//
// Define a minimal physics list for G4, just transportation.
// Used for debugging geometry.
//
// $Id: GaasPhysicsList.hh,v 1.3 2011/05/18 02:27:17 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:17 $
//
// Original author Rob Kutschke
//

#include "G4VUserPhysicsList.hh"

namespace mu2e {
  class GaasPhysicsList: public G4VUserPhysicsList {
  public:
    GaasPhysicsList();
    ~GaasPhysicsList();

  protected:

    // These methods are called by G4 not by users.
    void ConstructParticle();
    void ConstructProcess();
    void SetCuts();

  };

}  // end namespace mu2e

#endif /* Mu2eG4_GaasPhysicsList_hh */


