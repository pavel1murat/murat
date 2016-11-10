//
// Free function to create a geant4 "calorimetric" test environment geometry
//
// $Id: constructStudyEnv_v003.cc,v 1.2 2013/04/09 23:18:45 genser Exp $
// $Author: genser $
// $Date: 2013/04/09 23:18:45 $
//
// Original author KLG 
//
// Notes:
//
// one can nest volume inside other volumes if needed
// see other construct... functions for examples
//

// Mu2e includes.

#include "murat/gaas/constructGaasPbSandwich.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "ConfigTools/inc/SimpleConfig.hh"

// G4 includes
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4Colour.hh"
#include "G4Box.hh"
#include "G4NistManager.hh"

// c++ includes
#include <sstream>
#include <iostream>

using namespace std;

namespace mu2e {

  void constructGaasPbSandwich(VolumeInfo const & parentVInfo,SimpleConfig const & config) {

    const bool forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible");
    const bool doSurfaceCheck      = config.getBool("g4.doSurfaceCheck");
    const bool placePV             = true;

    // Extract calorimeter information from the config file and construct it

    // 2016-11-06 P.Murat: add GaAs and InGaAs

    G4NistManager* nistMan = G4NistManager::Instance();
    G4Material* mat;

    mat = G4Material::GetMaterial("InGaAs",false);
    if (mat == NULL) {
      G4Material* GaAs = new G4Material("InGaAs", 5.5*CLHEP::g/CLHEP::cm3, 3 );

      G4Element* In = nistMan->FindOrBuildElement("In",true);
      GaAs->AddElement( In, 1);

      G4Element* ga = nistMan->FindOrBuildElement("Ga",true);
      GaAs->AddElement( ga, 1);

      G4Element* as = nistMan->FindOrBuildElement("As",true);
      GaAs->AddElement( as, 1);
    }

    mat = G4Material::GetMaterial("GaAs",false);
    if (mat == NULL) {
      G4Material* GaAs = new G4Material("GaAs", 5.3*CLHEP::g/CLHEP::cm3, 2 );

      G4Element* ga = nistMan->FindOrBuildElement("Ga",true);
      GaAs->AddElement( ga, 1);

      G4Element* as = nistMan->FindOrBuildElement("As",true);
      GaAs->AddElement( as, 1);
    }

    mat = G4Material::GetMaterial("LYSO",false);
    if (mat == NULL) {
      G4Material* lyso = new G4Material("LYSO", 7.4*CLHEP::g/CLHEP::cm3, 4 );

      G4Element* Lu  = nistMan->FindOrBuildElement("Lu",true);
      G4Element* Si  = nistMan->FindOrBuildElement("Si",true);
      G4Element* O   = nistMan->FindOrBuildElement("O" ,true);
      G4Element* Y   = nistMan->FindOrBuildElement("Y" ,true);
      //      G4Element* Ce  = nistMan->FindOrBuildElement("Ce",true);
      
      lyso->AddElement( Lu, 71.0*CLHEP::perCent );
      lyso->AddElement( Si,  7.0*CLHEP::perCent );      
      lyso->AddElement( O , 18.0*CLHEP::perCent );      
      lyso->AddElement( Y ,  4.0*CLHEP::perCent );
    }

    // it will be a set of thin plates (G4Box'es)

    G4int  verbosityLevel  = config.getInt("calo.verbosityLevel",-1);
    G4bool isVisible       = config.getBool("calo.visible",true);
    G4bool forceSolid      = config.getBool("calo.solid",true);

    vector<double> tHL;
    config.getVectorDouble( "calo.transverseHalfLengths", tHL);

    vector<double> tCInParent;
    config.getVectorDouble( "calo.transverseCenterInWorld", tCInParent);

    vector<double> rMHL;
    config.getVectorDouble( "calo.moduleLayerHalfLengths",rMHL);

    vector<string> rMMat;
    config.getVectorString( "calo.moduleLayerMaterials", rMMat);

    int activeVolumeCopyNumber  = config.getInt("calo.activeVolumeStartingCopyNumber");
    int passiveVolumeCopyNumber = config.getInt("calo.passiveVolumeStartingCopyNumber");

    double    rMLCenterInParent = config.getDouble("calo.moduleStartingLongitPosition");
    int       rMNumberOfLayers  = config.getInt   ("calo.moduleNumberOfLayers");

    constructDoubleLayerdModule(parentVInfo,
                                tHL,
                                tCInParent,
                                rMHL,
                                rMMat,
                                rMNumberOfLayers,
                                verbosityLevel,
                                rMLCenterInParent,
                                passiveVolumeCopyNumber,
                                activeVolumeCopyNumber,
                                isVisible,
                                2, // G4Colour::Brown(),
                                3, // G4Colour::Yellow(),
                                forceSolid,
                                forceAuxEdgeVisible,
                                placePV,
                                doSurfaceCheck
                                );

  } // constructGaAsPbSandwich;

  void constructDoubleLayerdModule(VolumeInfo const & parentVInfo,
				   std::vector<double> const & tHL,
				   std::vector<double> const & tCInParent,
				   vector<double> const      & mHL,
				   vector<string> const      & mMat,
				   G4int                       mNumberOfLayers,
				   G4int                       verbosityLevel,
				   G4double                    mLCenterInParent,
				   G4int                       passiveVolumeStartingCopyNumber,
				   G4int                       activeVolumeStartingCopyNumber,
				   bool const                  isVisible,
				   G4Colour const & passiveVolumeColour,
				   G4Colour const &  activeVolumeColour,
				   bool const forceSolid,
				   bool const forceAuxEdgeVisible,
				   bool const placePV,
				   bool const doSurfaceCheck
				   ) { 
    //
    //  construct a calorimetric module (could be made to a n-layerd if needed)
    //  ( P - Passive and A - Active layers)

    G4int passiveVolumeCopyNumber = passiveVolumeStartingCopyNumber;
    G4int activeVolumeCopyNumber  = activeVolumeStartingCopyNumber;

    G4Material* mPMaterial = findMaterialOrThrow(mMat[0]); // passive
    G4Material* mAMaterial = findMaterialOrThrow(mMat[1]); // active

    G4double const mPParams[] = {tHL[0], tHL[1], mHL[0]}; // passive
    G4double const mAParams[] = {tHL[0], tHL[1], mHL[1]}; // active

    // we place all the volumes directly "in" the parent

    ostringstream vsPNumber("");
    vsPNumber.width(3);
    vsPNumber.fill('0');
    string vPName;

    ostringstream vsANumber("");
    vsANumber.width(3);
    vsANumber.fill('0');
    string vAName;

    G4double      mLayerStep  = 2.*(mHL[0]+mHL[1]);

    // we initialize the Z position to be one "step earlier" where it
    // should eventually be and then add the step to position of the
    // individual layers

    G4double      mPZCenterInParent = mLCenterInParent + mHL[0] - mLayerStep;
    G4double      mAZCenterInParent = mLCenterInParent - mHL[1];

    G4ThreeVector mPCenterInParent(tCInParent[0],tCInParent[1],mPZCenterInParent);
    G4ThreeVector mACenterInParent(tCInParent[0],tCInParent[1],mAZCenterInParent);

    for( G4int nl = 0; nl<mNumberOfLayers; ++nl ) {

      vsPNumber.str("");
      vsPNumber << passiveVolumeCopyNumber;
      vPName = "lead_" + vsPNumber.str();

      mPZCenterInParent += mLayerStep;
      mPCenterInParent.setZ(mPZCenterInParent);

      vsANumber.str("");
      vsANumber << activeVolumeCopyNumber;
      vAName = "gaas_" + vsANumber.str();

      mAZCenterInParent += mLayerStep;
      mACenterInParent.setZ(mAZCenterInParent);

      if (verbosityLevel > 0 ) {
        cout << __func__ << " constructing: " << vPName 
             << " " << "at " << mPCenterInParent 
             << endl;
        cout << __func__ << " constructing: " << vAName 
             << " " << "at " << mACenterInParent 
             << endl;
      }

      // passive medium
      VolumeInfo mPVInfo(nestBox( vPName,
                                   mPParams,
                                   mPMaterial,
                                   0x0, // no rotation
                                   mPCenterInParent,
                                   parentVInfo,
                                   passiveVolumeCopyNumber++,
                                   // non 0 for volume tracking purposes
                                   isVisible,
                                   passiveVolumeColour,
                                   forceSolid,
                                   forceAuxEdgeVisible,
                                   placePV,
                                   doSurfaceCheck
                                   ));

      // active medium
      VolumeInfo mAVInfo(nestBox( vAName,
                                   mAParams,
                                   mAMaterial,
                                   0x0, // no rotation
                                   mACenterInParent,
                                   parentVInfo,
                                   activeVolumeCopyNumber++,
                                   // non 0 for volume tracking purposes
                                   isVisible,
                                   activeVolumeColour,
                                   forceSolid,
                                   forceAuxEdgeVisible,
                                   placePV,
                                   doSurfaceCheck
                                   ));

    }
  } // constructDoubleLayerdModule

}
