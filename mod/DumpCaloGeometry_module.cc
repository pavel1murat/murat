//
// Print the information about the TTracker
//
// $Id: DumpCaloGeometry_module.cc,v 1.2 2015/03/09 00:56:23 murat Exp $
// $Author: murat $
// $Date: 2015/03/09 00:56:23 $
//
// Original author Rob Kutschke
//


#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Offline/CalorimeterGeom/inc/Disk.hh"
#include "Offline/CalorimeterGeom/inc/Crystal.hh"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"

// C++ includes.
#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

namespace mu2e {

  class DumpCaloGeometry : public art::EDAnalyzer {
  public:

    explicit DumpCaloGeometry(fhicl::ParameterSet const& pset);

    void analyze(const art::Event& e) override;

    void beginRun ( const art::Run& r) override;
  };

  DumpCaloGeometry::DumpCaloGeometry(fhicl::ParameterSet const& pset ):
    EDAnalyzer(pset){
  }

  void DumpCaloGeometry::analyze(const art::Event& ){}

  void DumpCaloGeometry::beginRun(const art::Run& run){

    GeomHandle<DiskCalorimeter> ch;
    const DiskCalorimeter* cal = ch.get();

    const CaloInfo&     ci = cal->caloInfo();

    int ndisks = cal->nDisks();

    printf("Calorimeter N(disks): %i\n", ndisks);

    printf("crystal halfLength     : %10.3f\n",ci.getDouble("crystalZLength")/2.);
    printf("crystal halfTrans      : %10.3f\n",ci.getDouble("crystalXYLength")/2.);
    printf("crystal wrap thickness : %10.3f\n",ci.getDouble("wrapperThickness"));
    printf("crystal case thickness : %10.3f\n",ci.getDouble("crystalFrameThickness"));

    for ( int i=0; i<ndisks; i++) {
      const Disk& disk = cal->disk(i);

      int ncrystals = disk.nCrystals();

      const DiskGeomInfo& gi = disk.geomInfo();

      printf(" -- id, ncrystals, Rin, Rout: %i, %3i, %10.4f, %10.4f",
	     disk.id(),ncrystals,gi.innerEnvelopeR(),gi.outerEnvelopeR());

      printf(" center: %12.3f %12.3f %12.3f\n",gi.origin().x(),gi.origin().y(),gi.origin().z());

      const Crystal cr0 = disk.crystal(0);

      printf ("crystal Z : %12.4f\n",cr0.position().z());
    }
  }


}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::DumpCaloGeometry)
