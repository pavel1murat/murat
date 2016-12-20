#ifndef murat_gaas_GaasStepLimiterPhysConstructor_hh
#define murat_gaas_GaasStepLimiterPhysConstructor_hh
//
// A Physics constructor that adds step limiters to some particles.
//
// $Id: GaasStepLimiterPhysConstructor.hh,v 1.3 2011/05/18 02:27:17 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:17 $
//
// Original author Rob Kutschke
//
#include "G4VPhysicsConstructor.hh"

namespace mu2e {

  class  GaasStepLimiterPhysConstructor: public G4VPhysicsConstructor {

  public:
    GaasStepLimiterPhysConstructor();
    ~GaasStepLimiterPhysConstructor();

    void ConstructParticle();
    void ConstructProcess();

  };

} // end namespace mu2e

#endif /* murat_gaas_GaasStepLimiterPhysConstructor_hh */
