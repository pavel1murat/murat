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
#include "Mu2eG4/inc/addStepLimiter.hh"

#include "G4ParticleTypes.hh"
#include "G4ProcessManager.hh"
#include "G4StepLimiter.hh"
#include "G4ParticleTable.hh"
#include "G4Alpha.hh"

namespace mu2e {
  GaasPhysicsList::GaasPhysicsList():  G4VUserPhysicsList(){
    defaultCutValue = 1.0*CLHEP::cm;
    SetVerboseLevel(1);
  }

  GaasPhysicsList::~GaasPhysicsList(){
  }

  void GaasPhysicsList::ConstructParticle(){

    G4ChargedGeantino::ChargedGeantinoDefinition();
    G4Electron::ElectronDefinition();
    G4Positron::PositronDefinition();
    G4MuonPlus::MuonPlusDefinition();
    G4MuonMinus::MuonMinusDefinition();
    G4Gamma::GammaDefinition();
    G4Proton::Definition();
    G4AntiProton::Definition();
    G4Alpha::Definition();
    G4GenericIon::GenericIonDefinition();

  }

  void GaasPhysicsList::ConstructProcess(){
    AddTransportation();

    // Add step limiters to a standard list of particles.

    vector<G4String> list;
    list.push_back( "e+"  );
    list.push_back( "e-"  );
    list.push_back( "mu+" );
    list.push_back( "mu-" );
    list.push_back( "pi+" );
    list.push_back( "pi-" );
    list.push_back( "kaon+"  );
    list.push_back( "kaon-"  );
    list.push_back( "proton" );
    list.push_back( "anti_proton"     );
    list.push_back( "chargedgeantino" );
    list.push_back( "alpha" );

    G4ParticleTable* ptable = G4ParticleTable::GetParticleTable();
    G4ParticleTable::G4PTblDicIterator* iter = ptable->GetIterator();
    
    // See note 1.
    iter->reset();
    
    // Check each existing particle to see if it is in the list.  See note 1.
    while( (*iter)() ){
      G4ParticleDefinition* particle = iter->value();
      G4ProcessManager* pmanager     = particle->GetProcessManager();
      G4String particleName          = particle->GetParticleName();
      
      // Is this particle in the list?
      if ( find( list.begin(), list.end(), particleName ) != list.end() ){

        // The process manager takes ownership of the G4StepLimiter object.
        pmanager->AddDiscreteProcess(new G4StepLimiter);
      }

    }
  }

  void GaasPhysicsList::SetCuts(){
  }

}  // end namespace mu2e
