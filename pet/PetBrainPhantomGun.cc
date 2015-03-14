//
// Shoots a single particle gun and puts its output into a generated event.
//
// $Id: PetBrainPhantomGun.cc,v 1.7 2013/12/02 02:03:24 murat Exp $
// $Author: murat $
// $Date: 2013/12/02 02:03:24 $
//
// Original author Rob Kutschke
// Modified by MyeongJae Lee. See docdb-2049
//

#include "murat/pet/PetGeomHandle.hh"
#include "murat/pet/BrainPhantom.hh"
#include "murat/pet/PetBrainPhantomGun.hh"

#include <iostream>

// Mu2e includes
#include "MCDataProducts/inc/PDGCode.hh"
#include "ConfigTools/inc/SimpleConfig.hh"

// Other external includes.
#include "CLHEP/Units/PhysicalConstants.h"


// Conversion energy for Al.  Should come from conditions.
// static const double pEndPoint = 104.96 * CLHEP::MeV;

//-----------------------------------------------------------------------------
PetBrainPhantomGun::PetBrainPhantomGun(art::Run const&, fhicl::ParameterSet& Parameters)
  : mu2e::GeneratorBase() 
{
//-----------------------------------------------------------------------------
// in general, mean depends on the dose, decide how many events to generate
// getEngine() comes from base class.
// this event generator is created from PetEventGenerator::beginRun, 
// at this point geometry should be already initialized and BrainPhantom - exist
// dose - in mCi (1 Ci = 3.7*10^10 decays per second)
//-----------------------------------------------------------------------------
  double   phantom_dz, phantom_r, phantom_body_r, phantom_body_dz;

  fRandFlat       = new CLHEP::RandFlat(getEngine());

  mu2e::PetGeomHandle<BrainPhantom> phantom;
  phantom_dz = phantom->Dz();
  phantom_r  = phantom->RMax();

  fDose         = Parameters.get<double> ("dose"        ,-1.);
  fNMean        = Parameters.get<double> ("nMean"       ,-1.);

  fHeadFraction = Parameters.get<double> ("headFraction",0.1);
  fBodyFraction = Parameters.get<double> ("bodyFraction",0.9);

  fDzMax        = Parameters.get<double> ("brainDz"     ,phantom_dz);
  fRMax         = Parameters.get<double> ("brainRadius" ,phantom_r );

  fBodyDzMax    = Parameters.get<double> ("bodyDz"      ,phantom_body_dz);
  fBodyRMax     = Parameters.get<double> ("bodyRadius"  ,phantom_body_r );

					// default: 1024 cells at 5 GHz
  fTimeWindow   = 1024*0.2;

  if (fDose > 0) {
//-----------------------------------------------------------------------------
// if fDose < 0, expect NMean also to be negative and defined explicitly
// fNMean: average number of the radioactive decays within the readout window,
//         currently set to 1024 cells , 200psec each
//-----------------------------------------------------------------------------
    fNMean       = fDose*3.7e10*1e-3*fTimeWindow*1e-9;
    fRandPoisson = new CLHEP::RandPoissonQ(getEngine(),fNMean);
  }
  else {
    fRandPoisson = 0;
  }

  fRandUnitSphere = new mu2e::RandomUnitSphere(getEngine());
}

//-----------------------------------------------------------------------------
// generate a random position within the phantom
//-----------------------------------------------------------------------------
void PetBrainPhantomGun::generatePos(CLHEP::Hep3Vector* Pos, int IBrain) {

  mu2e::PetGeomHandle<BrainPhantom> phantom;
  const CLHEP::Hep3Vector*          center;
  double                            r, dz, phi;

  if (IBrain == 1) {
//-----------------------------------------------------------------------------
// generate random point inside the brain
//-----------------------------------------------------------------------------
    center = phantom->Center();
    r      = fRMax*sqrt(fRandFlat->fire());
    dz     = (-1.+2.*fRandFlat->fire())*fDzMax;
    phi    = CLHEP::twopi*fRandFlat->fire();
    
    Pos->set(center->x()+r*cos(phi),center->y()+r*sin(phi),center->z()+dz );
  }
  else if (IBrain == 0) {
//-----------------------------------------------------------------------------
// generate random point inside the body
//-----------------------------------------------------------------------------
    center = phantom->BodyCenter();
    r      = phantom->BodyRMax()*sqrt(fRandFlat->fire());
    dz     = (-1.+2.*fRandFlat->fire())*phantom->BodyDz();
    phi    = CLHEP::twopi*fRandFlat->fire();
    
    Pos->set(center->x()+r*cos(phi),center->y()+r*sin(phi),center->z()+dz );
  }
}

//-----------------------------------------------------------------------------
// generate a set of back-to-back photons, 2  per annihilation
// for now, do not simulate the radioactive decay and the positron emission
//-----------------------------------------------------------------------------
void PetBrainPhantomGun::generate(mu2e::GenParticleCollection& GenParts) {

  CLHEP::Hep3Vector       pos;
  CLHEP::HepLorentzVector p2;

  int    n, ibrain; 
  double p, e, time, r1;

  p     = 0.511;
  e     = p;
  
  if (fNMean < 0) n = -fNMean;
  else            n = fRandPoisson->fire();

  for (int i=0; i<n; i++) {
//-----------------------------------------------------------------------------
// first figure whether it is the brain or the body
//-----------------------------------------------------------------------------
    r1 = fRandFlat->fire();

    if (r1 < fHeadFraction) ibrain = 1;
    else                    ibrain = 0;
//-----------------------------------------------------------------------------
// step 2: generate position and time
//-----------------------------------------------------------------------------
    generatePos(&pos,ibrain);
    time = fTimeWindow*fRandFlat->fire();

    // generate direction of the first photon, assume that the second is back-to-back

    CLHEP::HepLorentzVector p1(fRandUnitSphere->fire(p), e);

    p2.setX(-p1.x());
    p2.setY(-p1.y());
    p2.setZ(-p1.z());
    p2.setT( p1.t());
    
    GenParts.push_back( mu2e::GenParticle(mu2e::PDGCode::gamma, mu2e::GenId::particleGun, pos, p1, time));
    GenParts.push_back( mu2e::GenParticle(mu2e::PDGCode::gamma, mu2e::GenId::particleGun, pos, p2, time));
  }
}
