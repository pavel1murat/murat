//
// Print the information about the TTracker
//
// $Id: DumpTrackerGeometry_module.cc,v 1.2 2015/03/09 00:56:23 murat Exp $
// $Author: murat $
// $Date: 2015/03/09 00:56:23 $
//
// Original author Rob Kutschke
//


#include "GeometryService/inc/GeomHandle.hh"
#include "TTrackerGeom/inc/TTracker.hh"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"

// C++ includes.
#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

namespace mu2e {

  class DumpTrackerGeometry : public art::EDAnalyzer {
  public:

    explicit DumpTrackerGeometry(fhicl::ParameterSet const& pset);

    void analyze(const art::Event& e) override;

    void beginRun ( const art::Run& r) override;

  private:

  };

  DumpTrackerGeometry::DumpTrackerGeometry(fhicl::ParameterSet const& pset ):
    EDAnalyzer(pset){
  }

  void DumpTrackerGeometry::analyze(const art::Event& ){}

  void DumpTrackerGeometry::beginRun(const art::Run& run){

    TTracker const& tracker(*GeomHandle<TTracker>());

    int nstations = tracker.nDevices();

    printf("Tracker N(stations): %i\n", nstations);

    for ( auto const& dev : tracker.getDevices() ) {

    printf("----------------------------------------------------------------------------------------------------------------------\n");
    printf(" Chamber Sector Layer   StrawID  R(straw)  X(straw)   Y(straw)  Z(straw)   Rho(straw)    L/2       phi    wireNx   wireNy\n"); 
    printf("----------------------------------------------------------------------------------------------------------------------\n");
      //      printf("Station ID = %2i Z = %10.3f nsectors = %3i\n",dev.id(), dev.origin().z(),dev.nSectors());

      for ( auto const& sec : dev.getSectors() ) {
	int nlayers = sec.nLayers();

	for (int il=0; il<nlayers; il++) {

	  StrawId sid( sec.id(), il, 0 );

	  Straw const& straw = sec.getStraw(sid);

	  double x    = straw.getMidPoint().x();  
	  double y    = straw.getMidPoint().y();  
	  double rho  = sqrt(x*x+y*y);

	  double phi  = straw.direction().phi();
	  double nx   = straw.direction().x();
	  double ny   = straw.direction().y();

	  double z    = straw.getMidPoint().z();  
	  double hl   = straw.getHalfLength();
	  double r    = straw.getRadius();
	  
	  printf(" %5i %5i %5i",dev.id(), sec.id().getSector(), il);

	  double phi1 = phi/M_PI*180.;
	  printf("  %10i %8.3f %10.3f %10.3f %10.3f %10.3f %10.3f %8.2f %8.4f %8.4f\n",
		 straw.index().asInt(),r,x,y,z,rho, hl,phi1,nx,ny);
	}
      }
    }
  }


}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::DumpTrackerGeometry);
