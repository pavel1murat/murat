//
// Free function to create the Vane calorimeter.
//
// $Id: constructBrainPhantom.cc,v 1.2 2013/11/20 18:15:52 murat Exp $
// $Author: murat $
// $Date: 2013/11/20 18:15:52 $
//
// Original author Ivan Logashenko
// Modified by Bertrand Echenard
//
// Notes
//
//  1. a crystal has readouts at the back, both are surrounded by the wrapping, and the wrapping by a shell
//  2. by default, the wrapping surrounds the front/back face of the crystal+ro, the shell does not (shell is a casing)
//  3. The vanes are placed directly into DS3.  We did not make a mother volume for them.
//  4. The argument zOff is the zlocation of the center of the mother volume, as mesaured in the mu2e coordinate system.
//
//  5) Modified version  builds the calorimeter by making a physical mother volume.

#include <iostream>

// Mu2e includes.
#include "murat/pet/PetG4Helper.hh"
#include "G4Helper/inc/AntiLeakRegistry.hh"

#include "murat/pet/constructBrainPhantom.hh"
#include "murat/pet/PetGeometryService.hh"
#include "murat/pet/BrainPhantom.hh"
#include "murat/pet/PetGeomHandle.hh"

#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "murat/pet/PetSensitiveDetectorName.hh"

#include "murat/pet/PetnestBox.hh"
#include "murat/pet/PetnestTubs.hh"

#include "ConfigTools/inc/SimpleConfig.hh"

#include "CalorimeterGeom/inc/Vane.hh"
#include "CalorimeterGeom/inc/Crystal.hh"

#include "Mu2eG4/inc/CaloCrystalSD.hh"
#include "Mu2eG4/inc/CaloReadoutSD.hh"
#include "GeomPrimitives/inc/TubsParams.hh"


// G4 includes
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4Tubs.hh"

using namespace std;

namespace mu2e {

  PetVolumeInfo constructBrainPhantom(PetVolumeInfo const&  Mother, SimpleConfig const& config) {

    //-- Read parameters from config file
    //   Control of graphics for debugging the geometry.
    //   Only instantiate sectors to be drawn.

    //    int const verbosityLevel             = config.getInt("brainPhantom.verbosityLevel",0);

    //    bool const isVisible      = config.getBool("brainPhantom.calorimeterVisible",false);
    bool const isSolid        = config.getBool("brainPhantom.calorimeterSolid",false);

    bool const forceAuxEdgeVisible       = config.getBool("g4.forceAuxEdgeVisible",false);
    G4bool const doSurfaceCheck          = config.getBool("g4.doSurfaceCheck",false);
    bool const placePV                   = true;

    //calorimeter mother enveloppe
    double mother_radius   = config.getDouble("brainPhantom.motherRadius",110); 
    double mother_z0       = config.getDouble("brainPhantom.motherZ0",-100); 
    double mother_z1       = config.getDouble("brainPhantom.motherZ1",100); 

    //-- A helper class for parsing the config file.
    MaterialFinder materialFinder(config);
    G4Material*    material    = materialFinder.get("brainPhantom.material");

    //-- Get brain phantom handle
    const BrainPhantom* bf = PetGeomHandle<BrainPhantom>().get();

    //-- Construct brain phantom mother volume

    double mother_zlength  = mother_z1-mother_z0;
    double mother_zCenter  = (mother_z1+mother_z0)/2.0;
//-----------------------------------------------------------------------------
//  Make the mother volume for the brain phantom.
//-----------------------------------------------------------------------------
    // CLHEP::Hep3Vector const& posDS3         = Mother.centerInMu2e();
    // G4ThreeVector posBrainPhantomMother     = G4ThreeVector(posDS3.x(), 0, mother_zCenter);
    // G4ThreeVector posBrainPhantomMotherInDS = posBrainPhantomMother - posDS3;

    // TubsParams mother_params(0,mother_radius,mother_zlength/2.0, 0., CLHEP::twopi);

    // bool is_mother_solid = false;

    // PetVolumeInfo brainPhantomMotherInfo = PetnestTubs( "BrainPhantomMother",
    // 							mother_params,
    // 							findMaterialOrThrow("DSVacuum"),
    // 							0,
    // 							posBrainPhantomMotherInDS,
    // 							Mother,
    // 							0,
    // 							true,
    // 							G4Colour::Blue(),
    // 							is_mother_solid,
    // 							forceAuxEdgeVisible,
    // 							placePV,
    // 							doSurfaceCheck);

    G4ThreeVector  phantom_position(0.,0.,0);
    double         phantom_radius  = bf->RMax();
    double         phantom_zlength = bf->Dz()*2;

    TubsParams phantom_params(0,phantom_radius,phantom_zlength/2.0,0.,CLHEP::twopi);

    PetVolumeInfo brainPhantomInfo = PetnestTubs("BrainPhantom",
						 phantom_params,
						 material,
						 0,
						 phantom_position,
						 Mother,
						 0,
						 true,
						 G4Colour::Red(),
						 isSolid,
						 forceAuxEdgeVisible,
						 placePV,
						 doSurfaceCheck);

    G4ThreeVector  phantom_body_position(0.,0.,-500);
    double         phantom_body_radius  = bf->BodyRMax();
    double         phantom_body_zlength = bf->BodyDz()*2;

    TubsParams phantom_body_params(0,bf->BodyRMax(),bf->BodyDz(),0.,CLHEP::twopi);

    PetVolumeInfo bodyPhantomInfo = PetnestTubs("BodyPhantom",
						 phantom_body_params,
						 material,
						 0,
						 phantom_body_position,
						 Mother,
						 0,
						 true,
						 G4Colour::Red(),
						 isSolid,
						 forceAuxEdgeVisible,
						 placePV,
						 doSurfaceCheck);
    return Mother;

  }


} // end namespace mu2e
