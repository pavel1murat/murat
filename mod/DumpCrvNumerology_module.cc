//
// Print the CRV geometry
//

#include "GeometryService/inc/GeomHandle.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "CRVResponse/inc/CrvHelper.hh"

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

      int nmodules = modules->size();
      printf( "sector %2i name : %-10s nmodules: %3i\n", is , shield->getName().data(), nmodules);
      
      for(int im = 0; im<nmodules; im++) {

	const CRSScintillatorModule* module = &modules->at(im);

        const std::vector<CRSScintillatorLayer>* layers = &module->getLayers();

	int nlayers = layers->size();

	printf(" module: %3i nlayers : %3i\n", im, nlayers);

        for (int il=0; il<nlayers; il++) {

	  const CRSScintillatorLayer* layer = &layers->at(il);

          // const std::vector<double> & hl = layer->getHalfLengths();

          // const CLHEP::Hep3Vector &layerCenterInMu2e=layer->getPosition();

          const std::vector<std::shared_ptr<CRSScintillatorBar> >* bars = &layer->getBars();
          std::vector<std::shared_ptr<CRSScintillatorBar> >::const_iterator ibar;

	  int nbars = bars->size();

          for (int ib=0; ib<nbars; ib++) {
	    int ibb, iss, imm, ill;
	    CrvHelper::GetCrvCounterInfo(crvh, mu2e::CRSScintillatorBarIndex(bar_index), iss, imm, ill, ibb);
	    
	    if (ib == 0) {
	      printf(" ------------- layer : %3i nbars: %4i bar_index: %4i %4i %4i %4i %4i\n", 
		     il, nbars, bar_index, iss, imm, ill, ibb);
	    }
	    
	    bar_index++;
	  }
	}
      }
    }
  }
}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::DumpCrvNumerology);
