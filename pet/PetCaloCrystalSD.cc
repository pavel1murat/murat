//
// Define a sensitive detector for CaloCrystal Detectors
//
// $Id: PetCaloCrystalSD.cc,v 1.2 2013/10/16 18:37:26 murat Exp $
// $Author: murat $
// $Date: 2013/10/16 18:37:26 $
//
// Original author Ivan Logashenko
//

#include <cstdio>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

// Mu2e includes
#include "murat/pet/PetCaloCrystalSD.hh"
#include "Mu2eG4/inc/Mu2eG4UserHelpers.hh"
#include "Mu2eG4/inc/PhysicsProcessInfo.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "murat/pet/PetWorldG4.hh"

// G4 includes
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4ios.hh"

using namespace std;

namespace mu2e {

  PetCaloCrystalSD::PetCaloCrystalSD(G4String name, SimpleConfig const & config ): 
    PetSensitiveDetector(name,config)
  { }
  
  G4bool PetCaloCrystalSD::ProcessHits(G4Step* aStep,G4TouchableHistory*){

    // Calculate energy deposition in the crystal
    G4double edep = aStep->GetTotalEnergyDeposit();
    if( edep<=0 ) return false;

    _currentSize += 1;

    if( _sizeLimit>0 && _currentSize>_sizeLimit ) {
      if( (_currentSize - _sizeLimit)==1 ) {
        mf::LogWarning("G4") << "Maximum number of particles reached in " 
                             << SensitiveDetectorName
                             << ": "
                             << _currentSize << endl;
      }
      return false;
    }

    const G4TouchableHandle & touchableHandle = aStep->GetPreStepPoint()->GetTouchableHandle();

    // Get crystal ID
    G4int copyNo = touchableHandle->GetCopyNumber(0);

    ProcessCode endCode(_processInfo->findAndCount(Mu2eG4UserHelpers::findStepStoppingProcessName(aStep)));

    // Originally the hit position was saved in local crystal frame.
    // Not it is saved in Mu2e frame, hence the following code is
    // commented out.
    // Calculate enerdy deposition position along the crystal
     G4AffineTransform const& toLocal = touchableHandle->GetHistory()->GetTopTransform();
     //G4AffineTransform        toWorld = toLocal.Inverse();
     G4ThreeVector posWorld = aStep->GetPreStepPoint()->GetPosition();
     G4ThreeVector posLocal = toLocal.TransformPoint(posWorld);

    // Add the hit to the framework collection.
    // The point's coordinates are saved in the mu2e coordinate system.

    _collection->
      push_back(StepPointMC(_spHelper->particlePtr(aStep->GetTrack()),
                            copyNo,
                            edep,
                            aStep->GetNonIonizingEnergyDeposit(),
                            aStep->GetPreStepPoint()->GetGlobalTime(),
                            aStep->GetPreStepPoint()->GetProperTime(),
                            aStep->GetPreStepPoint()->GetPosition() - _mu2eOrigin,
                            aStep->GetPreStepPoint()->GetMomentum(),
                            aStep->GetStepLength(),
                            endCode
                            ));

    return true;

  }

} //namespace mu2e
