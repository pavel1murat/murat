// 2013-05-28 P.Murat

#include "murat/pet/PetEnvelope.hh"

#include <limits>
#include <algorithm>
#include <iterator>
#include <utility>
#include <cmath>

#include "cetlib/exception.h"

// #include "Mu2eBuildingGeom/inc/Mu2eBuilding.hh"
// #include "ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"
// #include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALBuilding.hh"

  //----------------------------------------------------------------
  PetEnvelope::PetEnvelope()
    : xmin_(0), xmax_(0), ymin_(0), ymax_(0), zmin_(0), zmax_(0)
  {}

  //----------------------------------------------------------------
  PetEnvelope::PetEnvelope(double XMin, double XMax, double YMin, double YMax, double ZMin, double ZMax) {

    xmin_ = XMin;
    xmax_ = XMax;
    ymin_ = YMin;
    ymax_ = YMax;
    zmin_ = ZMin;
    zmax_ = ZMax;
  }

  //================================================================
  std::ostream& operator<<(std::ostream& os, const PetEnvelope& env) {
    return os<<"PetEnvelope(xmin="<<env.xmin()
	     <<",xmax="<<env.xmax()
	     <<",ymin="<<env.ymin()
	     <<",ymax="<<env.ymax()
	     <<",zmin="<<env.zmin()
	     <<",zmax="<<env.zmax()
	     <<" )";
  }

  //================================================================
