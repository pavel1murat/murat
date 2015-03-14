//
// Construct the Mu2e G4 world and serve information about that world.
//
// $Id: PetWorld.cc,v 1.5 2013/11/20 18:15:52 murat Exp $
// $Author: murat $
// $Date: 2013/11/20 18:15:52 $
//
// Original author Rob Kutschke
//
//  Heirarchy is:
//  0      World (air)
//  1      Earthen Overburden
//  2      Concrete walls of the hall
//  3      Air inside the hall
//  4      Effective volume representing the DS coils+cryostats.
//  4      DS Vaccum
//
//  4      Effective volume representing the PS coils+cryostats.
//  4      PS Vacuum
//
//  The Earth overburden is modeled in two parts: a box that extends
//  to the surface of the earth plus a cap above grade.  The cap is shaped
//  as a G4Paraboloid.
//

// C++ includes
#include <iostream>
#include <vector>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

// Mu2e includes
#include "murat/pet/PetG4Helper.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"

#include "murat/pet/PetWorld.hh"
#include "murat/pet/constructPetWorldVolume.hh"

// #include "Mu2eG4/inc/constructVirtualDetectors.hh"
#include "murat/pet/PetconstructVisualizationRegions.hh"

#include "Mu2eG4/inc/MaterialFinder.hh"
#include "murat/pet/PetCaloCrystalSD.hh"
#include "murat/pet/PetCaloReadoutSD.hh"

#include "murat/pet/PetSensitiveDetector.hh"

#include "Mu2eG4/inc/findMaterialOrThrow.hh"

//#include "Mu2eG4/inc/nestCons.hh"

#include "Mu2eG4/inc/nestTubs.hh"

// #include "Mu2eG4/inc/nestTorus.hh"

#include "murat/pet/PetnestBox.hh"
#include "murat/pet/PetfinishNesting.hh"

#include "murat/pet/PetGeometryService.hh"
#include "murat/pet/PetWorldG4.hh"
#include "murat/pet/PetGeomHandle.hh"

#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "murat/pet/BrainImager.hh"

#include "murat/pet/constructBrainImager.hh"
#include "murat/pet/constructBrainPhantom.hh"

#include "murat/pet/PetSensitiveDetectorHelper.hh"

// G4 includes
#include "G4SDManager.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Paraboloid.hh"
#include "G4Colour.hh"
#include "G4Tubs.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4TwoVector.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ExactHelixStepper.hh"
#include "G4ChordFinder.hh"
#include "G4TransportationManager.hh"
#include "G4UserLimits.hh"
#include "G4ClassicalRK4.hh"
#include "G4ImplicitEuler.hh"
#include "G4ExplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4GDMLParser.hh"

#include "Mu2eG4/inc/Mu2eGlobalField.hh"
#include "Mu2eG4/inc/FieldMgr.hh"

using namespace std;

namespace mu2e {

  // This is the callback called by PetG4 via G4VPhysicalVolume* WorldMaker::Construct()
  G4VPhysicalVolume * PetWorld::construct(){
    // Construct all of the Mu2e world, hall, detectors, beamline ...
    return constructWorld();
  }

//-----------------------------------------------------------------------------
// Construct all of the Mu2e world, hall, detectors, beamline ...
//-----------------------------------------------------------------------------
  G4VPhysicalVolume* PetWorld::constructWorld(){

    _verbosityLevel = _config.getInt("world.verbosityLevel", 0);

    // If you play with the order of these calls, you may break things.
    PetGeomHandle<PetWorldG4> worldGeom;
    // G4ThreeVector tmpTrackercenter = GeomHandle<Mu2eBuilding>()->relicMECOOriginInMu2e() 
    //   + worldGeom->mu2eOriginInWorld()
    //   - G4ThreeVector(0.0,0.0,12000.-_config.getDouble("itracker.z0",0.0));

    instantiateSensitiveDetectors();

    PetVolumeInfo worldVInfo = constructPetWorldVolume(_config);

    if ( _verbosityLevel > 0) {
      cout << __func__ << " worldVInfo.centerInParent : " <<  worldVInfo.centerInParent << endl;
      cout << __func__ << " worldVInfo.centerInWorld  : " <<  worldVInfo.centerInWorld  << endl;
    }
//-----------------------------------------------------------------------------
// brain imager is always present, phantom - not necessarily
//-----------------------------------------------------------------------------
    PetVolumeInfo brainImagerInfo;
    PetVolumeInfo phantom;
    PetVolumeInfo const & world = _helper->locateVolInfo("World");

    if ( _config.getBool("hasBrainImager",false) ) {
      brainImagerInfo = constructBrainImager(world,_config);
    }

    if (_config.getBool("hasBrainPhantom",false)) {
//-----------------------------------------------------------------------------
// kludge: returns pointer to world...
//-----------------------------------------------------------------------------
      phantom = constructBrainPhantom(world,_config);
    }

    //    constructVirtualDetectors(_config); // beware of the placement order of this function
    
    PetconstructVisualizationRegions(worldVInfo, _config);

    sdHelper_->instantiateLVSDs(_config);

    if ( _verbosityLevel > 0) {
      mf::LogInfo log("GEOM");
      log << "Mu2e Origin:          " << worldGeom->mu2eOriginInWorld() << "\n";
    }

    constructStepLimiters();

    // Write out mu2e geometry into a gdml file.
    if ( _config.getBool("writeGDML",false) ) {
      string gdmlFileName = _config.getString("GDMLFileName","mu2e.gdml");
      G4GDMLParser parser;
      parser.Write(gdmlFileName, worldVInfo.logical);
    }

    return worldVInfo.physical;
  }


  // Adding a step limiter is a two step process.
  // 1) In the physics list constructor add a G4StepLimiter to the list of discrete
  //    physics processes attached to each particle species of interest.
  //
  // 2) In this code, create a G4UserLimits object and attach it to the logical
  //    volumes of interest.
  // The net result is specifying a step limiter for pairs of (logical volume, particle species).
  //
  void PetWorld::constructStepLimiters(){

    // Maximum step length, in mm.
    //    double maxStep = _config.getDouble("bfield.maxStep", 20.);

    //    AntiLeakRegistry& reg = art::ServiceHandle<PetG4Helper>()->antiLeakRegistry();
  }

//-----------------------------------------------------------------------------
  void PetWorld::instantiateSensitiveDetectors(){

    art::ServiceHandle<PetGeometryService> geom;

    G4SDManager* SDman      = G4SDManager::GetSDMpointer();

    // G4 takes ownership and will delete the detectors at the job end

    if(sdHelper_->enabled(StepInstanceName::virtualdetector)) {
      PetSensitiveDetector* vdSD =
        new PetSensitiveDetector(    SensitiveDetectorName::VirtualDetector(), _config);
      SDman->AddNewDetector(vdSD);
    }
//-----------------------------------------------------------------------------
// this is an illustration how using templates creates a problem - one has to 
// instantiate a class with a specific class name...
//-----------------------------------------------------------------------------
    if (   const_cast<PetGeometryService&>(_geom).hasElement<BrainImager>() ) {
      if(sdHelper_->enabled(StepInstanceName::calorimeter)) {
        PetCaloCrystalSD* ccSD     =
          new PetCaloCrystalSD(          SensitiveDetectorName::CaloCrystal(),     _config);
        SDman->AddNewDetector(ccSD);
      }

      if(sdHelper_->enabled(StepInstanceName::calorimeterRO)) {
        PetCaloReadoutSD* crSD     =
          new PetCaloReadoutSD(          SensitiveDetectorName::CaloReadout(),     _config);
        SDman->AddNewDetector(crSD);
      }
    }
  }

} // end namespace mu2e
