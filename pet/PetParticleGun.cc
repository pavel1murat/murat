//
// Shoots a single particle gun and puts its output into a generated event.
//
// $Id: PetParticleGun.cc,v 1.5 2013/11/04 21:09:32 murat Exp $
// $Author: murat $
// $Date: 2013/11/04 21:09:32 $
//
// Original author Rob Kutschke
// Modified by MyeongJae Lee. See docdb-2049
//

#include "murat/pet/PetParticleGun.hh"

#include <iostream>

// Mu2e includes
#include "MCDataProducts/inc/PDGCode.hh"
#include "ConfigTools/inc/SimpleConfig.hh"

// Other external includes.
#include "CLHEP/Units/PhysicalConstants.h"


// Conversion energy for Al.  Should come from conditions.
static const double pEndPoint = 104.96 * CLHEP::MeV;

//-----------------------------------------------------------------------------
PetParticleGun::PetParticleGun(art::Run const&, const mu2e::SimpleConfig& config)
  : mu2e::GeneratorBase()

    // Random number distributions; getEngine() comes from base class.
  , _randomUnitSphere( getEngine())

{
}

//-----------------------------------------------------------------------------
// generate 2 back-to-back photons
//-----------------------------------------------------------------------------
void PetParticleGun::generate( mu2e::GenParticleCollection& genParts) {

  //  m_gun.generate(genParts);

  CLHEP::Hep3Vector pos(0.,0.,0.);

  double p = 0.511;
  double e = p;

  CLHEP::HepLorentzVector p1( _randomUnitSphere.fire(p), e);

  // hack - for debugging only

  // double phi   = 2.5/120; // want to hit center of the crystal
  // double theta = 2.5/120; // want to hit center of the crystal

  // double px, py, pz;

  // px = p*cos(theta)*cos(phi);
  // py = p*cos(theta)*sin(phi);
  // pz = p*sin(theta);

  // CLHEP::HepLorentzVector p1( px, py, pz, e);

  CLHEP::HepLorentzVector p2;

  p2.setX(-p1.x());
  p2.setY(-p1.y());
  p2.setZ(-p1.z());
  p2.setT( p1.t());

  //  int pdg = 22; // photon

  double time = 0;

  genParts.push_back( mu2e::GenParticle(mu2e::PDGCode::gamma, mu2e::GenId::particleGun, pos, p1, time));
  genParts.push_back( mu2e::GenParticle(mu2e::PDGCode::gamma, mu2e::GenId::particleGun, pos, p2, time));
}
