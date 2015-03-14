#ifndef BrainPhantomGeom_BrainPhantomMaker_hh
#define BrainPhantomGeom_BrainPhantomMaker_hh
//
// Construct and return an Target.
//
//
// $Id: BrainPhantomMaker.hh,v 1.1 2013/11/04 21:09:32 murat Exp $
// $Author: murat $
// $Date: 2013/11/04 21:09:32 $
//
// Original author Peter Shanahan
//

#include <memory>
#include <string>
#include <vector>

#include "CLHEP/Vector/ThreeVector.h"
#include "ConfigTools/inc/SimpleConfig.hh"

class BrainPhantom;

class BrainPhantomMaker {

private:

  void BuildIt();
  void PrintConfig();// prints the target configuration as TargetMaker
                     // understands it...

  //  SimpleConfig const& _config;

  // pointer to the head phantom 

  std::unique_ptr<BrainPhantom>  fPhantom;

  // variables needed to build the Brain Phantom.  Read in from config file,
  // data base, etc.  These (and the TargetMaker object itself) only need
  // to persist long enough to make the Target.  After that, the definition
  // resides entirely in the Target object.

  CLHEP::Hep3Vector _detSysOrigin;

  std::string  fMaterial; // material of enclosing cylinder

public:

  BrainPhantomMaker();
  BrainPhantomMaker(mu2e::SimpleConfig const& config);
  ~BrainPhantomMaker();

  // Use compiler-generated copy c'tor, copy assignment

  // This is the accessor that will remain.
  //  std::unique_ptr<BrainPhantom> getBrainPhantomPtr() { return std::move(fPhantom); }
  std::unique_ptr<BrainPhantom> detectorPtr() { return std::move(fPhantom); }

};

#endif /* BrainPhantomGeom_BrainPhantomMaker_hh */
