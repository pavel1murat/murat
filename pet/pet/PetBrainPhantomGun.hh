#ifndef EventGenerator_PetBrainPhantomGun_hh
#define EventGenerator_PetBrainPhantomGun_hh
//
// Shoots a single particle gun and puts its output into a generated event.
//
// $Id: PetBrainPhantomGun.hh,v 1.4 2013/12/02 02:03:24 murat Exp $
// $Author: murat $
// $Date: 2013/12/02 02:03:24 $
//
// Original author Rob Kutschke
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

class PetBrainPhantomGun: public mu2e::GeneratorBase {
private:
  //  mu2e::ParticleGunImpl m_gun;
  
  mu2e::RandomUnitSphere* fRandUnitSphere;

  CLHEP::RandFlat*        fRandFlat;
  CLHEP::RandPoissonQ*    fRandPoisson;

  double fNMean;
  double fDose;		                // in mCi, if negative: fixed number of events
  double fHeadFraction;
  double fBodyFraction;	                // so far, 1-fHeadFraction, by construction
  double fDzMax;
  double fRMax;
  double fBodyDzMax;
  double fBodyRMax;
  double fTimeWindow; 			// digitization time window
 
public:
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
  PetBrainPhantomGun( art::Run const& run, fhicl::ParameterSet&);
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  double DzMax     () { return fDzMax; }
  double RMax      () { return fRMax;  }
  double BodyDzMax () { return fBodyDzMax; }
  double BodyRMax  () { return fBodyRMax;  }
  double TimeWindow() { return fTimeWindow;}
//-----------------------------------------------------------------------------
// setters
//-----------------------------------------------------------------------------
  void SetDzMax(double DzMax) { fDzMax = DzMax; }
  void SetRMax (double RMax ) { fRMax  = RMax;  }

  void SetBodyDzMax(double DzMax) { fBodyDzMax = DzMax; }
  void SetBodyRMax (double RMax ) { fBodyRMax  = RMax;  }

  void SetTimeWindow(double TMax) { fTimeWindow  = TMax;  }
//-----------------------------------------------------------------------------
// other methods
//-----------------------------------------------------------------------------
  void generatePos(CLHEP::Hep3Vector* Pos, int Brain);

					// adds generated particles to the collection

  virtual void generate(mu2e::GenParticleCollection& out);

};

#endif /* murat_pet_PetBrainPhantomGun_hh */
