#ifndef EventGenerator_PetParticleGun_hh
#define EventGenerator_PetParticleGun_hh
//
// Shoots a single particle gun and puts its output into a generated event.
//
// $Id: PetParticleGun.hh,v 1.2 2013/06/05 00:47:07 murat Exp $
// $Author: murat $
// $Date: 2013/06/05 00:47:07 $
//
// Original author Rob Kutschke
// Modified  by MyeongJae Lee. See docdb-2049
//
// The position is given in the Mu2e coordinate system.
//
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandExponential.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "Mu2eUtilities/inc/RandomUnitSphere.hh"

#include "EventGenerator/inc/GeneratorBase.hh"
// #include "EventGenerator/inc/ParticleGunImpl.hh"

// Forward references.
namespace art{ class Run; }

  // Forward reference.

namespace mu2e {
  class SimpleConfig;
}

class PetParticleGun: public mu2e::GeneratorBase {
  
public:
  PetParticleGun( art::Run const& run, const mu2e::SimpleConfig& config );

  // adds generated particles to the collection
  virtual void generate(mu2e::GenParticleCollection& out);
  
private:
  //  mu2e::ParticleGunImpl m_gun;
  
  mu2e::RandomUnitSphere _randomUnitSphere;
};

#endif /* murat_pet_PetParticleGun_hh */
