// This define a box, in the Mu2e coordinate system, which contains
// all the pieces of "real" mu2e geometry, such as detectors,
// building, and shielding dirt.
//
// cloned from Andrei Gaponenko, 2012

#ifndef PETENVELOPE_HH
#define PETENVELOPE_HH

#include <ostream>

#include "Mu2eInterfaces/inc/Detector.hh"
#include "art/Persistency/Common/Wrapper.h"

namespace mu2e {
  class PetGeometryService;
}

//   class Mu2eBuilding;
//   class ProtonBeamDump;
//   class ExtMonFNALBuilding;

class PetEnvelope : virtual public mu2e::Detector {
public:

  double xmin() const { return xmin_; }
  double xmax() const { return xmax_; }
  double ymin() const { return ymin_; }
  double ymax() const { return ymax_; }
  double zmin() const { return zmin_; }
  double zmax() const { return zmax_; }

  //----------------------------------------------------------------
  private:
  // Private ctr: the class should be only obtained via GeometryService
  friend class mu2e::PetGeometryService;

  PetEnvelope(double XMin, double XMax, double YMin, double YMax, double ZMin, double ZMax);

  template<class T> friend class art::Wrapper; // Needed for persistency
  PetEnvelope(); // Needed for persistency

  double xmin_;
  double xmax_;
  double ymin_;
  double ymax_;
  double zmin_;
  double zmax_;
  
};

std::ostream& operator<<(std::ostream& os, const PetEnvelope& env);

#endif/*PETENVELOPE_HH*/
