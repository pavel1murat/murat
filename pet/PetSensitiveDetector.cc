//
// Defines sensitive detector for a typicaly numbered volume using mu2e reference frame
//
// $Id: PetSensitiveDetector.cc,v 1.2 2013/10/15 23:41:14 murat Exp $
// $Author: murat $
// $Date: 2013/10/15 23:41:14 $
//
// Original author KLG
//

#include <cstdio>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

// Mu2e includes
#include "murat/pet/PetSensitiveDetector.hh"
#include "Mu2eG4/inc/Mu2eG4UserHelpers.hh"
#include "Mu2eG4/inc/PhysicsProcessInfo.hh"
#include "Mu2eG4/inc/SimParticleHelper.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "murat/pet/PetGeomHandle.hh"
#include "murat/pet/PetWorldG4.hh"

// G4 includes
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4ios.hh"

using namespace std;

namespace mu2e {

  PetSensitiveDetector::PetSensitiveDetector(G4String const name, const SimpleConfig& config):
    G4VSensitiveDetector(name),
    _collection(0),
    _processInfo(0),
    _mu2eOrigin(PetGeomHandle<PetWorldG4>()->mu2eOriginInWorld()),
    _debugList(0),
    _sizeLimit(config.getInt("g4.stepsSizeLimit",0)),
    _currentSize(0),
    _spHelper()
  {

   // Get list of events for which to make debug printout.
    // we generate the key string based on the detector name
    // consult murat/pet/SensitiveDetectorName.hh for the names

    std::ostringstream sdKeyName;
    sdKeyName<<"g4."<< SensitiveDetectorName << "SDEventList";
    //G4cout << __func__ << " sdKeyName: " << sdKeyName.str() << G4endl;
 
    string key(sdKeyName.str());
    if ( config.hasName(key) ){
      vector<int> list;
      config.getVectorInt(key,list);
      _debugList.add(list);
    }

  }

  void PetSensitiveDetector::Initialize(G4HCofThisEvent* HCE){

    _currentSize=0;

  }


  G4bool PetSensitiveDetector::ProcessHits(G4Step* aStep,G4TouchableHistory*){

    _currentSize += 1;

    if ( _sizeLimit>0 && _currentSize>_sizeLimit ) {
      if( (_currentSize - _sizeLimit)==1 ) {
        mf::LogWarning("G4") << "Maximum number of particles reached in " 
                             << SensitiveDetectorName
                             << ": "
                             << _currentSize << endl;
      }
      return false;
    }

    if ( _debugList.inList() )  {
            G4cout<<"edep "<<aStep->GetTotalEnergyDeposit()
                  <<" nidep "<<aStep->GetNonIonizingEnergyDeposit()
                  <<" step "<<aStep->GetStepLength()<<G4endl;
            G4cout<<"Step vol name "<<aStep->GetTrack()->GetVolume()->GetName()<<G4endl;
    }


    // Which process caused this step to end?
    ProcessCode endCode(_processInfo->
                        findAndCount(Mu2eG4UserHelpers::findStepStoppingProcessName(aStep)));

    // Add the hit to the framework collection.
    // The point's coordinates are saved in the mu2e coordinate system.
    _collection->
      push_back(StepPointMC(_spHelper->particlePtr(aStep->GetTrack()),
                            aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(),
                            aStep->GetTotalEnergyDeposit(),
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


  void PetSensitiveDetector::EndOfEvent(G4HCofThisEvent*){

    if( _sizeLimit>0 && _currentSize>=_sizeLimit ) {
      mf::LogWarning("G4") << "Total of " << _currentSize << " "
                           << SensitiveDetectorName
                           << " hits were generated in the event."
                           << endl
                           << "Only " << _sizeLimit << " are saved in output collection."
                           << endl;

      G4cout << "Total of " << _currentSize << " "
           << SensitiveDetectorName
           << " hits were generated in the event."
           << G4endl
           << "Only " << _sizeLimit << " are saved in output collection."
           << G4endl;

    }

    if (verboseLevel>0) {
      G4int NbHits = _collection->size();
      G4cout << "\n-------->Hits Collection: in this event they are " << NbHits
             << " hits in " << SensitiveDetectorName << ": " << G4endl;
      for (G4int i=0;i<NbHits;i++) (*_collection)[i].print(G4cout, true, false);
    }

  }


  void PetSensitiveDetector::beforeG4Event(StepPointMCCollection& outputHits,
                                            PhysicsProcessInfo& processInfo,
                                            const SimParticleHelper& spHelper){
    _collection  = &outputHits;
    _processInfo = &processInfo;
    _spHelper    = &spHelper;

    return;

  }

} //namespace mu2e
