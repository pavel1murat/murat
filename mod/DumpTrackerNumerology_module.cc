//
// Print the Tracker geometry
//

#include "GeometryService/inc/GeomHandle.hh"
#include "TrackerGeom/inc/Tracker.hh"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"

// C++ includes.
#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

namespace mu2e {

  class DumpTrackerNumerology : public art::EDAnalyzer {
  public:

    explicit DumpTrackerNumerology(fhicl::ParameterSet const& pset);

    void analyze(const art::Event& e) override;

    void beginRun ( const art::Run& r) override;

  private:

  };

  DumpTrackerNumerology::DumpTrackerNumerology(fhicl::ParameterSet const& pset ):
    EDAnalyzer(pset){
  }

  void DumpTrackerNumerology::analyze(const art::Event& ){}

  void DumpTrackerNumerology::beginRun(const art::Run& run){

    GeomHandle<Tracker> handle;

    const Tracker*  tracker = handle.get();

    int             nstations, nplanes, iface, npanels, nlayers;
    //    int             nstraws;

    //    const Station*  station;

    const Straw* straw;
    StrawId      sid;

    nstations = StrawId::_nstations; // tracker->nStations();

    printf("Tracker N(stations): %i\n", nstations);

    for (int ist=0; ist<nstations; ist++) {

      //      station = &tracker->getStations().at(ist);

      printf("---------------------------------------------------------------------------------");
      printf("--------------------------------------------------------------------\n");
      printf(" Station Plane Face Panel Layer Straw  StrawID   R(straw)  X(straw)");
      printf("  Y(straw)   Z(straw)   Rho(straw)     L/2       phi    wireNx   wireNy\n"); 
      printf("---------------------------------------------------------------------------------");
      printf("--------------------------------------------------------------------\n");
      //      printf("Station ID = %2i Z = %10.3f nsectors = %3i\n",dev.id(), dev.origin().z(),dev.nSectors());

      nplanes = 2; // station->nPlanes();

      for (int iplane=0; iplane<nplanes; iplane++) {
	int ipl = nplanes*ist+iplane;

	const Plane* plane = &tracker->getPlane(ipl);

	npanels = plane->nPanels();

	for (int ipanel=0; ipanel<npanels; ipanel++) {

	  const Panel* panel = &plane->getPanel(ipanel);

	  iface   = ipanel%2;
	  nlayers = panel->nLayers();

	  for (int il=0; il<nlayers; il++) {
	    //	    const Layer* lay = &panel->getLayer(il);
	      
	    // nstraws = zl->nStraws();

	    // 	      for (int is=0; is<nstraws; is++) {
	    int is = 0;
	    straw       = &panel->getStraw(il);
	    sid         = straw->id();
	    
	    double x    = straw->getMidPoint().x();  
	    double y    = straw->getMidPoint().y();  
	    double rho  = sqrt(x*x+y*y);

	    double phi  = straw->direction().phi();
	    double nx   = straw->direction().x();
	    double ny   = straw->direction().y();
	    
	    double z    = straw->getMidPoint().z();  
	    double hl   = straw->halfLength();
	    double r    = straw->getRadius();
	    
	    double phi1 = phi/M_PI*180.;
	    printf("  %3i %6i %5i %4i %5i %5i %10i %8.3f %10.3f %10.3f %10.3f %10.3f %10.3f %8.2f %8.4f %8.4f\n",
		   ist,iplane,iface, ipanel,
		   il,is, straw->id().asUint16(),r,x,y,z,rho, hl,phi1,nx,ny);
	  }
	}
      }
    }
  }


}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::DumpTrackerNumerology);
