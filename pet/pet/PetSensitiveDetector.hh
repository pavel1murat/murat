#ifndef murat_pet_PetSensitiveDetector_hh
#define murat_pet_PetSensitiveDetector_hh
//
// Defines a generic Pet sensitive detector
//
// $Id: PetSensitiveDetector.hh,v 1.2 2013/10/15 23:41:14 murat Exp $
// $Author: murat $
// $Date: 2013/10/15 23:41:14 $
//
// Original author KLG
//

// Mu2e includes
#include "Mu2eG4/inc/EventNumberList.hh"
#include "Mu2eG4/inc/SimParticleHelper.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"

// G4 includes
#include "G4VSensitiveDetector.hh"

// Art includes
#include "art/Persistency/Provenance/ProductID.h"
#include "art/Framework/Principal/Event.h"

class G4Step;
class G4HCofThisEvent;

namespace mu2e {

  // Forward declarations in mu2e namespace
  class SimpleConfig;
  class PhysicsProcessInfo;

  class PetSensitiveDetector : public G4VSensitiveDetector{

  public:

    PetSensitiveDetector(G4String const name, SimpleConfig const & config);

    virtual void Initialize(G4HCofThisEvent*);

    virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);

    virtual void EndOfEvent(G4HCofThisEvent*);

    void beforeG4Event(StepPointMCCollection& outputHits,
                       PhysicsProcessInfo & processInfo,
                       const SimParticleHelper& spHelper);

  protected:

    // Non-owning pointer to the  collection into which hits will be added.
    StepPointMCCollection* _collection;

    // Non-ownning pointer and object that returns code describing physics processes.
    PhysicsProcessInfo* _processInfo;

    // Mu2e point of origin
    G4ThreeVector _mu2eOrigin;

    // List of events for which to enable debug printout.
    EventNumberList _debugList;

    // Limit maximum size of the steps collection
    int _sizeLimit;
    int _currentSize;

    // A helper to create pointers to SimParticles
    const SimParticleHelper *_spHelper;
  };

} // namespace mu2e

#endif /* murat_pet_PetSensitiveDetector_hh */
