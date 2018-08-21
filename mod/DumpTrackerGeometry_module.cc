//
// Print the TTracker geometry
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

    int nplanes = tracker.nPlanes();

    printf("Tracker N(planes): %i\n", nplanes);

    for ( auto const& plane : tracker.getPlanes() ) {

    printf("----------------------------------------------------------------------------------------------------------------------\n");
    printf(" Chamber Sector Layer   StrawID  R(straw)  X(straw)   Y(straw)  Z(straw)   Rho(straw)    L/2       phi    wireNx   wireNy\n"); 
    printf("----------------------------------------------------------------------------------------------------------------------\n");

      for ( auto const& panel : plane.getPanels() ) {
	int nlayers = panel.nLayers();

	for (int il=0; il<nlayers; il++) {
//-----------------------------------------------------------------------------
// assume straw 'il' (0,1) corresponds to the layer 'il'
//-----------------------------------------------------------------------------
//	  StrawId sid( panel.id(), il);

	  Straw const& straw = panel.getStraw(il);

	  double x    = straw.getMidPoint().x();  
	  double y    = straw.getMidPoint().y();  
	  double rho  = sqrt(x*x+y*y);

	  double phi  = straw.direction().phi();
	  double nx   = straw.direction().x();
	  double ny   = straw.direction().y();

	  double z    = straw.getMidPoint().z();  
	  double hl   = straw.getHalfLength();
	  double r    = straw.getRadius();
	  
	  printf(" %5i %5i %5i",plane.id().asUint16(), panel.id().getPanel(), il);

	  double phi1 = phi/M_PI*180.;
	  printf("  %10i %8.3f %10.3f %10.3f %10.3f %10.3f %10.3f %8.2f %8.4f %8.4f\n",
		 straw.id().asUint16(),r,x,y,z,rho, hl,phi1,nx,ny);
	}
      }
    }
  }


}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::DumpTrackerGeometry);
