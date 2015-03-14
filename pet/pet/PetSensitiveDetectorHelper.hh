#ifndef murat_pet_PetSensitiveDetectorHelper_hh
#define murat_pet_PetSensitiveDetectorHelper_hh
//
// Some helper functions to manage repeated tasks related to sensitive detectors.
//
// $Id: PetSensitiveDetectorHelper.hh,v 1.3 2013/10/15 23:41:14 murat Exp $
// $Author: murat $
// $Date: 2013/10/15 23:41:14 $
//
// Original author Rob Kutschke
//

// From Mu2e
#include "murat/pet/PetSensitiveDetector.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StepInstanceName.hh"

// From the art tool chain
#include "fhiclcpp/ParameterSet.h"
#include "cetlib/maybe_ref.h"

// From C++ and STL
#include <map>
#include <memory>
#include <string>
#include <vector>

// Forward references.
namespace art   { class Event; }
namespace fhicl { class ParameterSet; }

namespace mu2e {

  class SimpleConfig;
  class SimParticleHelper;

  class PetSensitiveDetectorHelper{

  public:

    PetSensitiveDetectorHelper( fhicl::ParameterSet const& pset );
    // Accept compiler generator d'tor, copy c'tor and assignment operator.

    // Register the sensitive detector with this class; to be called after G4 Initialize.
    void registerSensitiveDetectors();

    // Create data  products and pre-fill with input hits if any; to be called at the start of each event.
    void createProducts(const art::Event& evt, const SimParticleHelper& spHelper);

    // Attach the new per-event data products to the corresponding sensitive detector objects.
    void updateSensitiveDetectors( PhysicsProcessInfo&   info,
                                   const SimParticleHelper& spHelper);

    // Put the data products into the event.
    void put( art::Event& event);

    // Return all of the instances names of the data products to be produced.
    std::vector<std::string> stepInstanceNamesToBeProduced() const;

    // Query the same info
    bool enabled(StepInstanceName::enum_type instance) const;

    // Return one of the StepPointMCCollections.
    cet::maybe_ref<StepPointMCCollection> steps( StepInstanceName::enum_type id );

    // create SDs for arbitrary logical volumes as requiested
    void instantiateLVSDs(const SimpleConfig& config);

  private:

    // A helper class to hold information about each sensitive detector object.
    struct StepInstance {
      explicit StepInstance(StepInstanceName::enum_type astepName ):
        p(),
        stepName(StepInstanceName(astepName).name()),
        sensitiveDetector(0){
      }

      explicit StepInstance(const std::string &name)  : stepName(name) {}

      explicit StepInstance():
        p(),
        stepName(),
        sensitiveDetector(nullptr){
      }

      // Accept compiler written d'tor, copy c'tor and assignment operator.

      // See note 2 in the .cc file for why the collection is not held by pointer.
      // For historical reasons the two names are different; maybe some day we will synchronize them.
      StepPointMCCollection    p;
      std::string              stepName;
      PetSensitiveDetector*    sensitiveDetector;

    };

    // Enabled pre-defined StepPointMC collections, except the timevd.
    typedef std::map<StepInstanceName::enum_type,StepInstance> InstanceMap;
    InstanceMap stepInstances_;

    // Logical volumes requested to be sensitive
    typedef std::map<std::string,StepInstance> LVSDMap;
    LVSDMap lvsd_;

    // Existing hit collections that should be merged into the output
    typedef std::vector<art::InputTag> InputTags;
    InputTags preSimulatedHits_;
  };


} // end namespace mu2e
#endif /* murat_pet_PetSensitiveDetectorHelper_hh */
