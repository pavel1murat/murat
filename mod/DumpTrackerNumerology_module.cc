//
// Print the Tracker geometry
//

#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"

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

    const Straw* straw;
    StrawId      sid;

    int nstations = StrawId::_nstations ; 
    int nplanes   = StrawId::_nplanes   ; 

    printf("Tracker N(stations): %i\n", nstations);
//-----------------------------------------------------------------------------
// station as a concept doesn't exist in the Mu2e offline software
// so loop over planels
//-----------------------------------------------------------------------------
    for (int ipl=0; ipl<nplanes; ipl++) {

      const Plane* plane = &tracker->getPlane(ipl);
      int iplane = plane->getPanel(0).getStraw(0).id().plane();
      int ist = iplane/2;

      int npanels = plane->nPanels();

      for (int ipanel=0; ipanel<npanels; ipanel++) {

        printf("---------------------------------------------------------------");
        printf("-----------------------------------------------------------------------\n");
        printf("Station Plane Face Panel Layer Straw StrawID R(straw)  X(straw)");
        printf("   Y(straw)   Z(straw)  Rho(straw)     L/2       phi    wireNx   wireNy\n"); 
        printf("---------------------------------------------------------------");
        printf("-----------------------------------------------------------------------\n");

        const Panel* panel = &plane->getPanel(ipanel);

        int iface   = panel->getStraw(0).id().face();
        int nstraws = panel->nStraws();

        for (int is=0; is<nstraws; is++) {
          int il = is % 2;
          straw       = &panel->getStraw(is);
          sid         = straw->id();
	    
          double x    = straw->getMidPoint().x();  
          double y    = straw->getMidPoint().y();  
          double rho  = sqrt(x*x+y*y);

          double phi  = straw->direction().phi();
          double nx   = straw->direction().x();
          double ny   = straw->direction().y();
	    
          double z    = straw->getMidPoint().z();  
          double hl   = straw->halfLength();
          double r    = tracker->strawProperties()._strawOuterRadius;;
	    
          double phi1 = phi/M_PI*180.;
          printf(" %3i %6i %5i %4i %5i %5i %8i %8.3f %10.3f %10.3f %10.3f %10.3f %10.3f %8.2f %8.4f %8.4f\n",
                 ist,iplane,iface, ipanel,
                 il,is, straw->id().asUint16(),r,x,y,z,rho, hl,phi1,nx,ny);
        }
      }
    }
  }


}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::DumpTrackerNumerology);
