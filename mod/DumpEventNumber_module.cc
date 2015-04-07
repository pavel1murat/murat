//
//  A module to look at the provenance of a product.
//
//  $Id: DumpEventNumber_module.cc,v 1.4 2013/10/21 21:15:46 kutschke Exp $
//  $Author: kutschke $
//  $Date: 2013/10/21 21:15:46 $
//
//  Original author Rob Kutschke
//

#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Provenance.h"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/ParameterSetRegistry.h"

#include <iostream>

using namespace std;

namespace mu2e {

  class DumpEventNumber : public art::EDAnalyzer {
  protected: 
    int       _ientry;

  public:
    explicit DumpEventNumber(fhicl::ParameterSet const& pset);

    void analyze(const art::Event& event) override;

  private:

    void printProvenance( art::Provenance const& );

  };

  DumpEventNumber::DumpEventNumber(fhicl::ParameterSet const& pset)
    : art::EDAnalyzer(pset){

    _ientry = 0;
  }


  void DumpEventNumber::analyze(const art::Event&  Event) {
    printf(" run, event  %6i %10i  ientry: %10i\n",
	   Event.run(),Event.event(),_ientry);
    _ientry++;

  } // end analyze

} // end namespace mu2e

using mu2e::DumpEventNumber;
DEFINE_ART_MODULE(DumpEventNumber);
