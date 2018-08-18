//
// A PhysicsConstructor that adds step limiters to some particles.
//
// $Id: GaasStepLimiterPhysConstructor.cc,v 1.2 2011/05/18 02:27:18 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:18 $
//
// Original author Rob Kutschke
//
// Notes
// 1) These methods are simple enough that one might be tempted to put
//    them in the .hh file.  But that will not work since they are never
//    called by code that we compile; except for the constructor, are
//    all methods callbacks that are called be G4.

// Mu2e includes
#include "murat/gaas/GaasStepLimiterPhysConstructor.hh"
#include "murat/gaas/addGaasStepLimiter.hh"

using namespace std;

namespace mu2e {

  GaasStepLimiterPhysConstructor::GaasStepLimiterPhysConstructor():
    G4VPhysicsConstructor("GaasStepLimiterPhysConstructor"){
  }

  GaasStepLimiterPhysConstructor::~GaasStepLimiterPhysConstructor(){
  }

  void GaasStepLimiterPhysConstructor::ConstructParticle(){
  }

  void GaasStepLimiterPhysConstructor::ConstructProcess(){
    addGaasStepLimiter();
  }


}  // end namespace mu2e
