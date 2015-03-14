#ifndef Mu2eG4_PrimaryGeneratorAction_hh
#define Mu2eG4_PrimaryGeneratorAction_hh
//
// Give generated tracks to G4 by copying information from a GenParticleCollection.
//
// $Id: PetPrimaryGeneratorAction.hh,v 1.2 2013/10/15 23:41:14 murat Exp $
// $Author: murat $
// $Date: 2013/10/15 23:41:14 $
//
// Original author Rob Kutschke
//

// C++ includes
#include <string>

// Framework includes
#include "art/Framework/Principal/Handle.h"

// G4 includes
#include "globals.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"

class G4ParticleDefinition;
class G4ParticleGun;
class G4Event;
class TH1D;

namespace mu2e {

  class SteppingAction;
  class SimParticlePrimaryHelper;

  typedef std::vector<art::ValidHandle<StepPointMCCollection> > HitHandles;

  class PetPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction{
  public:

    PetPrimaryGeneratorAction();

    // This is the interface specified by G4.
    void GeneratePrimaries(G4Event*);

    // Should change the interface for Primary
    void setEventData(const GenParticleCollection *gens, // may be NULL. No ownership passing.
                      const HitHandles& hitInputs,
                      SimParticlePrimaryHelper *parentMapping);

  private:

    void fromEvent( G4Event* );

    void addG4Particle(G4Event *event,
                       PDGCode::type pdgId,
                       const G4ThreeVector& pos,
                       double time,
                       double properTime,
                       const G4ThreeVector& mom);

    // Input event kinematics
    // Must be set before the call to GeneratePrimaries.
    const GenParticleCollection *genParticles_;
    const HitHandles *hitInputs_;
    SimParticlePrimaryHelper *parentMapping_;

    TH1D* _totalMultiplicity;

  };

}  // end namespace mu2e
#endif /* Mu2eG4_PrimaryGeneratorAction_hh */
