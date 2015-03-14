//
// Geometry and identifier info about the BrainPhantom.
//
//
// $Id: BrainPhantom.cc,v 1.2 2013/11/20 18:15:52 murat Exp $
// $Author: murat $
// $Date: 2013/11/20 18:15:52 $
//
// Original author R. Bernstein and Rob Kutschke
//
//C++ includes
#include <algorithm>

//mu2e includes
#include "murat/pet/BrainPhantom.hh"

//-----------------------------------------------------------------------------
// by default, the brain phantom is small - need the geometry file!
//-----------------------------------------------------------------------------
BrainPhantom::BrainPhantom(): 
  fCenter    (0.,0.,0),
  fBodyCenter(0.,0.,0)
{
  fDz       = 0.0;
  fRMax     = 0.0;
  fBodyDz   = 0.0;
  fBodyRMax = 0.0;
}

//-----------------------------------------------------------------------------
BrainPhantom::~BrainPhantom() {
}
