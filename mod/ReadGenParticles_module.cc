// Reads out all GenParticles in a collection.
//
// $Id: ReadGenParticles_module.cc,v 1.2 2013/10/21 20:44:04 genser Exp $
// $Author: genser $
// $Date: 2013/10/21 20:44:04 $
//
// Original author Andrei Gaponenko
//

#include <iostream>
#include <string>

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"

#include "MCDataProducts/inc/GenParticle.hh"
<<<<<<< HEAD
#include "MCDataProducts/inc/GenParticle.hh"
=======
>>>>>>> 34081a1 (work with the last version of OFFLINE)

#include "TH1F.h"

namespace mu2e {

  class ReadGenParticles : public art::EDAnalyzer {
  public:

    explicit         ReadGenParticles(fhicl::ParameterSet const& pset);

    void             analyze(const art::Event& e);
    virtual void     beginJob();

  private:
					// The two strings specify what collection to print
    std::string  _moduleLabel;
    std::string  _instanceName;

    struct Hist_t {
      TH1F*   _momentum;                     // 
      TH1F*   _energy;                       // 
      TH1F*   _ekin;                         // 
    } ;

    Hist_t _hist;
  };

  ReadGenParticles::ReadGenParticles(const fhicl::ParameterSet& pset) : 
    art::EDAnalyzer(pset),
    _moduleLabel (pset.get<std::string>("inputModuleLabel" )),
    _instanceName(pset.get<std::string>("inputInstanceName"))
  {}


  void ReadGenParticles::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;

    _hist._momentum  = tfs->make<TH1F>("mom" ,"Momentum"      ,  500,  0,  500);
    _hist._energy    = tfs->make<TH1F>("e"   ,"Energy"        , 1000,900, 1100);
    _hist._ekin      = tfs->make<TH1F>("ekin","Kinetic Energy", 1000,  0,  200);
  }


  void ReadGenParticles::analyze(const art::Event& event) {

    art::Handle<GenParticleCollection> ih;
    event.getByLabel(_moduleLabel, _instanceName, ih);
    const GenParticleCollection& gp(*ih);

    for (GenParticleCollection::const_iterator i = gp.begin(); i != gp.end(); ++i) {
      double p    = (*i).momentum().vect().mag();
      double e    = (*i).momentum().e();
      double ekin = e-(*i).momentum().m();
      _hist._momentum->Fill(p);
      _hist._energy->Fill(e);
      _hist._ekin->Fill(ekin);
    }

  } // analyze()

}  // end namespace mu2e

// Register the module with the framework
//using mu2e::ReadGenParticles;
DEFINE_ART_MODULE(mu2e::ReadGenParticles);
