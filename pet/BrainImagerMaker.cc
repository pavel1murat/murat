//
// Make a Vane Calorimeter.
//
// $Id: BrainImagerMaker.cc,v 1.5 2013/11/04 21:09:32 murat Exp $
// $Author: murat $
// $Date: 2013/11/04 21:09:32 $

// original authors Julie Managan and Robert Bernstein

//
// C++ includes
#include <math.h>
#include <boost/shared_ptr.hpp>


// Mu2e includes
#include "murat/pet/BrainImagerMaker.hh"
#include "murat/pet/BrainImager.hh"
#include "CalorimeterGeom/inc/Vane.hh"

// Framework include files
#include "cetlib/exception.h"
#include "Mu2eUtilities/inc/hep3VectorFromStdVector.hh"

// other includes
#include "CLHEP/Vector/RotationX.h"
#include "CLHEP/Vector/RotationY.h"
#include "CLHEP/Vector/RotationZ.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Matrix/DiagMatrix.h"


BrainImagerMaker::BrainImagerMaker( mu2e::SimpleConfig const& config) {

  _calo = std::unique_ptr<BrainImager>(new BrainImager());

  _calo->_nSections             = config.getInt   ("calorimeter.numberOfVanes");
  _calo->_caseThickness         = config.getDouble("calorimeter.caseThickness");
  _calo->_crystalHW             = config.getDouble("calorimeter.crystalHalfTrans");
  _calo->_crystalHL             = config.getDouble("calorimeter.crystalHalfLong");
  _calo->_nCrystalR             = config.getInt   ("calorimeter.nCrystalRSlices");
  _calo->_nCrystalZ             = config.getInt   ("calorimeter.nCrystalZSlices");
  _calo->_rMin                  = config.getDouble("calorimeter.rInscribed");
  _calo->_rMax                  = _calo->_rMin + 2.0*(_calo->_crystalHW + _calo->_wrapperThickness + _calo->_shellThickness) * _calo->_nCrystalR;

  _calo->_enveloppeRadius       = config.getDouble("calorimeter.caloMotherRadius",850); 
  _calo->_enveloppeZ0           = config.getDouble("calorimeter.caloMotherZ0",11740); 
  _calo->_enveloppeZ1           = config.getDouble("calorimeter.caloMotherZ1",13910); 

  _calo->_shieldHalfThickness   = config.getDouble("calorimeter.shieldHalfThickness");
  _calo->_absorberHalfThickness = config.getDouble("calorimeter.neutronAbsorberHalfThickness");

  _calo->_wrapperThickness      = config.getDouble("calorimeter.crystalWrapperThickness",0.0);
  _calo->_shellThickness        = config.getDouble("calorimeter.crystalShellThickness",0.0);

  _calo->_nROPerCrystal         = config.getInt("calorimeter.crystalReadoutChannelCount");
  _calo->_roHalfTrans           = config.getDouble("calorimeter.crystalReadoutHalfTrans");
  _calo->_roHalfThickness       = config.getDouble("calorimeter.crystalReadoutHalfThickness");
  
  _calo->_nonUniformity         = config.getDouble("calorimeter.crystalNonUniformity",0.0);
  _calo->_timeGap               = config.getDouble("calorimeter.timeGap",100.0);
  _calo->_electronEdep          = config.getDouble("calorimeter.electronDepositionAPD",1000.0);
  _calo->_electronEmin          = config.getDouble("calorimeter.electronMinEnergyAPD",0.1);

  _calo->_apdMeanNoise          = config.getDouble("calorimeter.meanNoiseAPD", 0.0);
  _calo->_apdSigmaNoise         = config.getDouble("calorimeter.sigmaNoiseAPD", 0.03);
  _calo->_lysoLightYield        = config.getDouble("calorimeter.lysoLightYield", 2000.0);
  _calo->_apdQuantumEff         = config.getDouble("calorimeter.quantumEffAPD", 0.68);
  _calo->_lightCollectEffAPD    = config.getDouble("calorimeter.lightCollectEffAPD", 0.11);


  //COORDINATE OF VOLUME CONTAINING THE CRYSTALS ONLY (NO READOUT,..) W.R.T OUTSIDE VANE VOLUME
  double absorberHalfLength =  _calo->_shieldHalfThickness + _calo->_absorberHalfThickness;
  _calo->_crystalShift = CLHEP::Hep3Vector(_calo->_crystalHL, 0 ,absorberHalfLength);
  

  //THE CALORIMETER ORIGIN IS TAKEN AS THE POINT CLOSEST TO THE TRACKER IN MU2E COORDINATES
  double xOrigin               = 0; // -config.getDouble("mu2e.solenoidOffset");
  double zOrigin               = config.getDouble("calorimeter.calorimeterZOrigin",11740);
  _calo->_origin               = CLHEP::Hep3Vector(xOrigin,0,zOrigin);
	
	
  //make sure above information is consistent
  CheckIt();
	
  // Create vanes
  MakeVanes();
}


BrainImagerMaker::~BrainImagerMaker() {}




void BrainImagerMaker::MakeVanes() {

  // Local reference frame of the vane:
  // Z is along Mu2e z
  // Y is radial towards large radius
  // X is transverse to vane plane

  double absorberHalfLength =  _calo->_shieldHalfThickness + _calo->_absorberHalfThickness;
  //  double caloHalfLength     =  (_calo->_crystalHW + _calo->_wrapperThickness + _calo->_shellThickness)*_calo->_nCrystalZ;
          
  double crystalCellRadius = _calo->_crystalHW + _calo->_wrapperThickness + _calo->_shellThickness;
  CLHEP::Hep3Vector crystalShift(_calo->_roHalfThickness,0,0);

  double dX     =  _calo->_crystalHL + _calo->_wrapperThickness + _calo->_roHalfThickness + _calo->_caseThickness;
  double dR     = (_calo->_crystalHW + _calo->_wrapperThickness + _calo->_shellThickness) * _calo->_nCrystalR + _calo->_caseThickness;
  double dZ     = (_calo->_crystalHW + _calo->_wrapperThickness + _calo->_shellThickness) * _calo->_nCrystalZ + absorberHalfLength + _calo->_caseThickness;
  double radius = _calo->_rMin + dR - _calo->_caseThickness;
  double dphi   = 2*CLHEP::pi/_calo->_nSections;


  for (unsigned int i=0; i<_calo->_nSections; ++i ) {
    double phi = i*dphi;

    boost::shared_ptr<mu2e::Vane> thisVane( new mu2e::Vane(i,_calo->_rMin,_calo->_nCrystalR,_calo->_nCrystalZ, crystalCellRadius, crystalShift) );	 
    _calo->_sections.push_back(thisVane);

    CLHEP::Hep3Vector localOrigin(radius*cos(phi),radius*sin(phi),0.); // absorberHalfLength+caloHalfLength+_calo->_caseThickness);
    
    thisVane->setSize(        CLHEP::Hep3Vector(dX,dR,dZ) );
    //             thisVane->setRotation(    (CLHEP::HepRotation::IDENTITY)*CLHEP::HepRotationZ(-CLHEP::pi/2 - i*dphi) );
    thisVane->setRotation(    (CLHEP::HepRotation::IDENTITY)*CLHEP::HepRotationZ(CLHEP::pi-i*dphi));
    thisVane->setOriginLocal( localOrigin );
    thisVane->setOrigin(      localOrigin + _calo->origin() );
  }

}



void BrainImagerMaker::CheckIt(void) {

  //check that calorimeter fits in the mother volume
  double absorberHalfLength =  _calo->_shieldHalfThickness + _calo->_absorberHalfThickness;
  double dR  = (_calo->_crystalHW + _calo->_wrapperThickness + _calo->_shellThickness) * _calo->_nCrystalR + _calo->_caseThickness;
  double dZ  = (_calo->_crystalHW + _calo->_wrapperThickness + _calo->_shellThickness) * _calo->_nCrystalZ + absorberHalfLength + _calo->_caseThickness;
  double calozBegin = _calo->_origin.z();
  double calozEnd   = _calo->_origin.z() + 2*dZ;
  
  std::cout<<dZ<<" "<<calozEnd<<" "<<(_calo->_crystalHW + _calo->_wrapperThickness + _calo->_shellThickness)<<std::endl;
  
  if ( (_calo->_rMin + 2*dR) > _calo->_enveloppeRadius) 
    {throw cet::exception("VaneCaloGeom") << "calorimeter radius larger than calorimeter mother \n";} 
  if (calozBegin < _calo->_enveloppeZ0 || calozBegin > _calo->_enveloppeZ1) 
    {throw cet::exception("VaneCaloGeom") << "calorimeter.calorimeterZOrigin   outside calorimeter mother.\n";}  
  if (calozEnd   > _calo->_enveloppeZ1)                       
    {throw cet::exception("VaneCaloGeom") << "calorimeter z-coordinate extends outside calorimeter mother.\n";}  
  
	   
  // Check number of readouts
  int nRO    = _calo->_nROPerCrystal;
  if( ! (nRO==1 || nRO==2 || nRO==4) ) 
    {throw cet::exception("VaneCaloGeom") << "calorimeter.crystalReadoutChannelCount can only be 1,2 or 4.\n";}


  // Check size of readouts
  if (nRO == 1) {
    if( _calo->_roHalfTrans > _calo->_crystalHW ) {
      throw cet::exception("VaneCaloGeom") << 
	"calorimeter.crystalReadoutHalfTrans > calorimeter.crystalHalfTrans.\n";
    }
  } 
  else {
    if( _calo->_roHalfTrans > 0.5*_calo->_crystalHW ) {
      throw cet::exception("VaneCaloGeom") << 
	"calorimeter.crystalReadoutHalfTrans > 0.5*calorimeter.crystalHalfTrans.\n";
    }
  }
}

