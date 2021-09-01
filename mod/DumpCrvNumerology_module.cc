//
// Print the CRV geometry
//

#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/CRVResponse/inc/CrvHelper.hh"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"

// C++ includes.
#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

namespace mu2e {

  class DumpCrvNumerology : public art::EDAnalyzer {

  public:
    explicit DumpCrvNumerology(fhicl::ParameterSet const& pset);
    void     analyze(const art::Event& e) override;
    void     beginRun ( const art::Run& r) override;

  private:

  };

  DumpCrvNumerology::DumpCrvNumerology(fhicl::ParameterSet const& pset ):
    EDAnalyzer(pset){
  }

  void DumpCrvNumerology::analyze(const art::Event& ){}
//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
  void DumpCrvNumerology::beginRun(const art::Run& run) {

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
}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::DumpCrvNumerology);
