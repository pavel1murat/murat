#ifndef murat_pet_BrainPhantom_hh
#define murat_pet_BrainPhantom_hh
//
// Hold all geometry and identifier information about
// a Calorimeter.  In order to insulate this class from
// knowledge of databases etc, this class must not know
// how to make itself.
//
// $Id: BrainPhantom.hh,v 1.2 2013/11/20 18:15:52 murat Exp $
// $Author: murat $
// $Date: 2013/11/20 18:15:52 $
//
// Original author R. Bernstein and Rob Kutschke
//

//C++ includes
#include <vector>
#include <boost/shared_ptr.hpp>

// Mu2e includes
#include "Mu2eInterfaces/inc/Detector.hh"

//CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"


using boost::static_pointer_cast;
using boost::shared_ptr;

class BrainPhantom: public mu2e::Detector {

  friend class BrainPhantomMaker;

protected:
					// brain
  CLHEP::Hep3Vector fCenter;
  double            fRMax;
  double            fDz;
					// body
  CLHEP::Hep3Vector fBodyCenter;
  double            fBodyRMax;
  double            fBodyDz;

public:
  
  BrainPhantom ();
  ~BrainPhantom();

  const CLHEP::Hep3Vector* Center    () const { return (const CLHEP::Hep3Vector*)  &fCenter;     }
  const CLHEP::Hep3Vector* BodyCenter() const { return (const CLHEP::Hep3Vector*)  &fBodyCenter; }

  double   Dz  () const { return fDz;   }
  double   RMax() const { return fRMax; }

  double   BodyDz  () const { return fBodyDz;   }
  double   BodyRMax() const { return fBodyRMax; }
//-----------------------------------------------------------------------------
// overloaded methods
//-----------------------------------------------------------------------------
};

#endif /* murat_pet_BrainPhantom_hh */
