#ifndef Mu2eG4_PetSensitiveDetectorName_hh
#define Mu2eG4_PetSensitiveDetectorName_hh
// Define names of Sensitive Detectors; revised to forward the names of the
// StepInstanceName names.
//
// $Id: PetSensitiveDetectorName.hh,v 1.1 2013/06/06 23:02:21 murat Exp $
// $Author: murat $
// $Date: 2013/06/06 23:02:21 $
//
// Original author KLG

#include "MCDataProducts/inc/StepInstanceName.hh"

namespace mu2e {

  class PetSensitiveDetectorName {

  public:

    // we define the functins here to avoid a run time undefined symbol error
    static char const * TrackerGas(){
      return StepInstanceName::name(StepInstanceName::tracker).c_str();
    }

    static char const * VirtualDetector(){
      return StepInstanceName::name(StepInstanceName::virtualdetector).c_str();
    }

    static char const * CaloCrystal(){
      return StepInstanceName::name(StepInstanceName::calorimeter).c_str();
    }

    static char const * CaloReadout(){
      return StepInstanceName::name(StepInstanceName::calorimeterRO).c_str();
    }

    static char const * ExtMonFNAL(){
      return "ExtMonFNAL";
    }

    static char const * ExtMonUCITof(){
      return StepInstanceName::name(StepInstanceName::ExtMonUCITof).c_str();
    }

    static char const * StoppingTarget(){
      return StepInstanceName::name(StepInstanceName::stoppingtarget).c_str();
    }

    static char const * CRSScintillatorBar(){
      return StepInstanceName::name(StepInstanceName::CRV).c_str();
    }

    static char const * TTrackerDeviceSupport(){
      return StepInstanceName::name(StepInstanceName::ttrackerDS).c_str();
    }

    static char const * ProtonAbsorber() {
      return StepInstanceName::name(StepInstanceName::protonabsorber).c_str();
    }

    static char const * PSVacuum() {
      return StepInstanceName::name(StepInstanceName::PSVacuum).c_str();
    }

    static char const * TrackerSWires(){
      return StepInstanceName::name(StepInstanceName::trackerSWires).c_str();
    }

    static char const * ITrackerFWires(){
      return StepInstanceName::name(StepInstanceName::itrackerFWires).c_str();
    }

    static char const * TrackerWalls(){
      return StepInstanceName::name(StepInstanceName::trackerWalls).c_str();
    }

  };

} // namespace mu2e

#endif /* Mu2eG4_PetSensitiveDetectorName_hh */
