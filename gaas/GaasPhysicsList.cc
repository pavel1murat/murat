//
// Define a minimal physics list.
// Just transportation; used for debugging geometry.
//
// $Id: GaasPhysicsList.cc,v 1.6 2013/10/25 18:47:09 genser Exp $
// $Author: genser $
// $Date: 2013/10/25 18:47:09 $
//
// Original author Rob Kutschke
//

// CLHEP includes
#include "CLHEP/Units/SystemOfUnits.h"

// Geant4 includes
#include "murat/gaas/GaasPhysicsList.hh"
#include "murat/gaas/GaasStepLimiterPhysConstructor.hh"

#include "G4ParticleTypes.hh"
#include "G4ProcessManager.hh"
#include "G4StepLimiter.hh"
#include "G4ParticleTable.hh"
#include "G4Alpha.hh"
#include "G4String.hh"

#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4IonQMDPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4HadronElasticPhysicsLEND.hh"
#include "G4DataQuestionaire.hh"

#include <vector>

using namespace std;

namespace mu2e {
  GaasPhysicsList::GaasPhysicsList():  G4VModularPhysicsList(){
    defaultCutValue = 1.0*CLHEP::cm;
    SetVerboseLevel(1);

    int verbose = 1;
    
    G4DataQuestionaire it(photon, neutron, radioactive);

    // EM Physics
    this->RegisterPhysics( new G4EmStandardPhysics(verbose));

    // Synchroton Radiation & GN Physics
    this->RegisterPhysics( new G4EmExtraPhysics(verbose) );

    // Decays 
    this->RegisterPhysics( new G4DecayPhysics(verbose) );
    //if ( rad == true ) this->RegisterPhysics( new G4RadioactiveDecayPhysics(verbose) );
    this->RegisterPhysics( new G4RadioactiveDecayPhysics(verbose) );

    // The modular physics list takes ownership of the StepLimiterPhysConstructor.
    this->RegisterPhysics( new GaasStepLimiterPhysConstructor() );
  }

  GaasPhysicsList::~GaasPhysicsList(){
  }

    
  void GaasPhysicsList::SetCuts() {
  }

}  // end namespace mu2e
