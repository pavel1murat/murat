//
// Print the Tracker geometry
//

#include "art/Framework/Principal/Event.h"
//#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"

#include "Offline/GeometryService/inc/GeomHandle.hh"

#include "Offline/GeometryService/inc/VirtualDetector.hh"
#include "Offline/DataProducts/inc/VirtualDetectorId.hh"

#include "Offline/TrackerGeom/inc/Tracker.hh"

#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/CRVResponse/inc/CrvHelper.hh"

#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Offline/CalorimeterGeom/inc/Disk.hh"
#include "Offline/CalorimeterGeom/inc/Crystal.hh"

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
	    double r    = tracker->strawOuterRadius();
	    
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
    GeomHandle<CosmicRayShield> crvh; 

    const CosmicRayShield*  crv = crvh.get();

    const std::vector<CRSScintillatorShield>* shields = &crv->getCRSScintillatorShields();

    int nsectors = shields->size(); 
    printf("CRV N(sectors): %i\n", nsectors);

    int bar_index = 0;
      
    for (int is=0; is<nsectors; is++) { 
      const CRSScintillatorShield* shield = &shields->at(is);

      const std::vector<CRSScintillatorModule>* modules = &shield->getCRSScintillatorModules();
      const CRSScintillatorModule*    m0    = &modules->at(0);
      const CRSScintillatorLayer*     l0    = &m0->getLayers().at(0);
      const std::shared_ptr<CRSScintillatorBar> b0 = l0->getBars().at(0);
      const CRSScintillatorBarDetail* bd0   = &b0->getBarDetail();

      int nmodules = modules->size();
      
      printf( "sector %2i name : %-10s nmodules: %3i nlayers: %3lu nbars: %3lu ",
	      is,shield->getName().data(),
	      nmodules,
	      m0->getLayers().size(),
	      l0->getBars().size());
      
      printf(" ----------------------- bar_index: %4i %4i %4i %4i  %10.3f %10.3f %10.3f\n", 
	     bar_index,
	     bd0->getThicknessDirection(),bd0->getWidthDirection(),bd0->getLengthDirection(),
	     b0->getHalfThickness(),b0->getHalfWidth(),b0->getHalfLength());

      for(int im = 0; im<nmodules; im++) {

	const CRSScintillatorModule* module = &modules->at(im);

        const std::vector<CRSScintillatorLayer>* layers = &module->getLayers();

	int nlayers = layers->size();

	printf(" module: %3i nlayers : %3i\n", im, nlayers);

        for (int il=0; il<nlayers; il++) {

	  const CRSScintillatorLayer* layer = &layers->at(il);
	  const CLHEP::Hep3Vector* lpos =  &layer->getPosition();
	  
          // const std::vector<double> & hl = layer->getHalfLengths();

          // const CLHEP::Hep3Vector &layerCenterInMu2e=layer->getPosition();

          const std::vector<std::shared_ptr<CRSScintillatorBar> >* bars = &layer->getBars();
          std::vector<std::shared_ptr<CRSScintillatorBar> >::const_iterator ibar;

	  int nbars = bars->size();

          for (int ib=0; ib<nbars; ib++) {
	    const std::shared_ptr<CRSScintillatorBar> bar = bars->at(ib);
	    int ibb, iss, imm, ill;
	    CrvHelper::GetCrvCounterInfo(crvh, mu2e::CRSScintillatorBarIndex(bar_index), iss, imm, ill, ibb);
	    
	    if (ib == 0) {
	      printf(" ----------------------- layer : %3i nbars: %4i bar_index: %4i %4i %4i %4i %4i %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n", 
		     il, nbars, bar_index, iss, imm, ill, ibb,
		     lpos->x(),lpos->y(),lpos->z(),
		     layer->getHalfThickness(),
		     layer->getHalfWidth(),
		     layer->getHalfLength()
		     );
	    }

	    const CRSScintillatorBarDetail* bd = &bar->getBarDetail();
	    
	    printf("%5i %3i %3i %2i %3i %2i %2i %2i %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n",
		   bar_index,is,im,il,ib,
		   bd->getThicknessDirection(),bd->getWidthDirection(),bd->getLengthDirection(),
		   bar->getPosition().x(),bar->getPosition().y(),bar->getPosition().z(),
		   bar->getHalfThickness(),bar->getHalfWidth(),bar->getHalfLength());

	    bar_index++;
	  }
	}
      }
    }
  }

//-----------------------------------------------------------------------------
  void DumpGeometry::dumpCaloGeometry() {
    GeomHandle<DiskCalorimeter> ch;
    const DiskCalorimeter* cal = ch.get();

    const CaloInfo&     ci = cal->caloInfo();

    int ndisks = cal->nDisk();

    printf("Calorimeter N(disks): %i\n", ndisks);

    printf("crystal halfLength     : %10.3f\n",ci.getDouble("crystalZLength")/2.);
    printf("crystal halfTrans      : %10.3f\n",ci.getDouble("crystalXYLength")/2.);
    printf("crystal wrap thickness : %10.3f\n",ci.getDouble("wrapperThickness"));
    printf("crystal case thickness : %10.3f\n",ci.getDouble("crystalFrameThickness"));

    for ( int i=0; i<ndisks; i++) {
      const Disk& disk = cal->disk(i);

      int ncrystals = disk.nCrystals();

      printf(" -- id, ncrystals, Rin, Rout: %i, %3i, %10.4f, %10.4f",
	     disk.id(),ncrystals,disk.innerRadius(),disk.outerRadius());

      const DiskGeomInfo& gi = disk.geomInfo();

      printf(" center: %12.3f %12.3f %12.3f\n",gi.origin().x(),gi.origin().y(),gi.origin().z());

      const Crystal cr0 = disk.crystal(0);

      printf ("crystal Z : %12.4f\n",cr0.position().z());
    }
  }

//-----------------------------------------------------------------------------  
  void DumpGeometry::beginRun(const art::Run& run) {

    if (_dumpVirtualDetectors  != 0) dumpVirtualDetectors ();
    if (_dumpTrackerNumerology != 0) dumpTrackerNumerology();
    if (_dumpCRVNumerology     != 0) dumpCRVNumerology    ();
    if (_dumpCaloGeometry      != 0) dumpCaloGeometry     ();

  }
}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::DumpGeometry);
