//
// Print the Tracker geometry
//

#include "art/Framework/Principal/Event.h"
//#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"
#include "TrackerGeom/inc/Tracker.hh"

// C++ includes.
#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

namespace mu2e {

  class DumpGeometry : public art::EDAnalyzer {
  public:

    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<int> dumpVirtualDetectors { Name("dumpVirtualDetectors" ),  Comment("dump virtual detectors"   ) }; 
      fhicl::Atom<int> dumpTrackerNumerology{ Name("dumpTrackerNumerology"),  Comment("dump tracker numerology"  ) }; 
      fhicl::Atom<int> dumpCRVNumerology    { Name("dumpCRVNumerology"    ),  Comment("dump CRV     numerology"  ) }; 
      fhicl::Atom<int> dumpCaloGeometry     { Name("dumpCaloGeometry"     ),  Comment("dump calorimeter geometry") }; 
    };

    class VirtualDetectorA : public VirtualDetector {
    public:
      const std::map<int,CLHEP::Hep3Vector>*  map_local() { return &_local; }
    };

    explicit DumpGeometry(const art::EDAnalyzer::Table<Config>& config);

    void  analyze (const art::Event& e) override;
    void  beginRun(const art::Run&   r) override;

    void  dumpVirtualDetectors ();
    void  dumpTrackerNumerology();
    void  dumpCRVNumerology    ();
    void  dumpCaloGeometry     ();

  private:
    int  _dumpVirtualDetectors;
    int  _dumpTrackerNumerology;
    int  _dumpCRVNumerology;
    int  _dumpCaloGeometry;
  };

//-----------------------------------------------------------------------------
  DumpGeometry::DumpGeometry(const art::EDAnalyzer::Table<Config>& config):
    EDAnalyzer(config),
    _dumpVirtualDetectors (config().dumpVirtualDetectors ()),
    _dumpTrackerNumerology(config().dumpTrackerNumerology()),
    _dumpCRVNumerology    (config().dumpCRVNumerology    ()),
    _dumpCaloGeometry     (config().dumpCaloGeometry     ())
  {
  }

//-----------------------------------------------------------------------------
  void DumpGeometry::analyze(const art::Event& ){}

//-----------------------------------------------------------------------------
  void DumpGeometry::dumpVirtualDetectors() {
    GeomHandle<VirtualDetector> handle;

    int ndet = handle->nDet();
    
    printf(" [DumpGeometry::beginRun] N(virtual detectors) = %i\n",ndet);

    VirtualDetectorA*  vd = (VirtualDetectorA*) handle.get();
    printf("  id name                                        x(loc)   y(loc)   z(loc)");
    printf("  x(glob)  y(glob)  z(glob)\n");
    printf(" ------------------------------------------------------------------------");
    printf(" --------------------------\n");

    for (auto const &pair: *vd->map_local()) {
      int i                             = pair.first;
      std::string name                  = vd->volumeName(i);
      const CLHEP::Hep3Vector& xyz_loc  = vd->getLocal(i);
      const CLHEP::Hep3Vector& xyz_glob = vd->getGlobal(i);

      printf("%4i %-40s  %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n",
	     i,name.data(),
	     xyz_loc.x(),xyz_loc.y(),xyz_loc.z(),
	     xyz_glob.x(),xyz_glob.y(),xyz_glob.z());

    }
  }

//-----------------------------------------------------------------------------
  void DumpGeometry::dumpTrackerNumerology() {
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

//-----------------------------------------------------------------------------
  void DumpGeometry::dumpCRVNumerology() {
  }

//-----------------------------------------------------------------------------
  void DumpGeometry::dumpCaloGeometry() {
  }

//-----------------------------------------------------------------------------  
  void DumpGeometry::beginRun(const art::Run& run) {

    if (_dumpVirtualDetectors  != 0) dumpVirtualDetectors ();
    if (_dumpTrackerNumerology != 0) dumpTrackerNumerology();
    if (_dumpCRVNumerology     != 0) dumpTrackerNumerology();
    if (_dumpCaloGeometry      != 0) dumpCaloGeometry     ();

  }
}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::DumpGeometry);
