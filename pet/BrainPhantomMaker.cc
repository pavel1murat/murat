///////////////////////////////////////////////////////////////////////////////
// Construct and return an Target.
//
//
// $Id: BrainPhantomMaker.cc,v 1.4 2013/11/30 03:41:58 murat Exp $
// $Author: murat $
// $Date: 2013/11/30 03:41:58 $
//
// Original author Peter Shanahan
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <cmath>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "ConditionsService/inc/GlobalConstantsHandle.hh"
#include "ConditionsService/inc/PhysicsParams.hh"

#include "murat/pet/BrainPhantomMaker.hh"
#include "murat/pet/BrainPhantom.hh"

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/RotationY.h"
#include "CLHEP/Vector/RotationZ.h"

using namespace std;

//-----------------------------------------------------------------------------
BrainPhantomMaker::BrainPhantomMaker() {
}

//-----------------------------------------------------------------------------
BrainPhantomMaker::BrainPhantomMaker(mu2e::SimpleConfig const& config) 
  : _detSysOrigin(0.,0.,0.) {

//-----------------------------------------------------------------------------
// create phantom, initialize its geometry, Maker is a C++ friend
//-----------------------------------------------------------------------------
  double    brain_z0, body_z0;

  fPhantom = unique_ptr<::BrainPhantom>(new ::BrainPhantom());

  // ####   fMaterial = config.getString("brainPhantom.material");

  fPhantom->fDz       = config.getDouble("brainPhantom.dz");
  fPhantom->fRMax     = config.getDouble("brainPhantom.radius");
  fPhantom->fBodyDz   = config.getDouble("brainPhantom.bodyDz");
  fPhantom->fBodyRMax = config.getDouble("brainPhantom.bodyRadius");

  brain_z0 = config.getDouble("brainPhantom.brainZ0");
  fPhantom->fCenter.set(0.,0.,brain_z0);
    
  body_z0 = config.getDouble("brainPhantom.bodyZ0");
  fPhantom->fBodyCenter.set(0.,0.,body_z0);
    
  // debugging print...
  //  if ( verbosity > 0 ) PrintConfig();
  
  // Do the real work.
  // BuildIt( );
}

//-----------------------------------------------------------------------------
BrainPhantomMaker::~BrainPhantomMaker() {
}

void BrainPhantomMaker::BuildIt() {
  
  // Build the Target Geometry.  This means MU2E internal geometry, not
  // Root, G4, or any other scheme.
  
  //  _targ->_centerInMu2e = CLHEP::Hep3Vector(0.,0.,0.);
  
  //  _targ->_fillMaterial=_fillMaterial;
} 

void BrainPhantomMaker::PrintConfig() {
  // printout the BrainPhantomMaker's understanding of what it needs to build.
  //  for debugging...

  std::cout<<"\n BrainPhantomMaker Input Configuration -----------------"<<std::endl;
}
