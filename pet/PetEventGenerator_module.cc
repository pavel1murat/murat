/*

  A plug_in for running a variety of event generators.

  $Id: PetEventGenerator_module.cc,v 1.8 2013/11/30 03:41:58 murat Exp $
  $Author: murat $
  $Date: 2013/11/30 03:41:58 $

  Original author Rob Kutschke

  Eventually this will support a variety of event generators, controllable
  from the run time configuration.  A given call might invoke one or more
  of these generators.

  1) A full featured single particle gun.
  2) Single conversion track, uniformly from the targets.
  3) (Emax-E)**5 DIO model.
  4) Other DIO models.
  5) protons, neutrons, gammas and nuclear fragments from muon capture.
  6) Mockups of CLHEP::pion capture on nuclei and of CLHEP::pion and muon decay in flight.
  I say mock-ups because I see this starting from an known CLHEP::pion and muon
  flux distributions, not by starting from a CLHEP::pion or a muon entering
  the DS.
  7) Simplified models of cosmics.

  At present I expect that the highest fidelity generation of cosmics will be
  done by running an external generator and then reading "events" from the output
  of that generator.  Perhaps the merge will be done in this module, perhaps
  it will be done in a separate module?

*/

// Mu2e includes.

#include "ConfigTools/inc/SimpleConfig.hh"
#include "ConfigTools/inc/requireUniqueKey.hh"
#include "MCDataProducts/inc/GenId.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"

// Particular generators that this code knows about.
#include "murat/pet/PetConversionGun.hh"
#include "murat/pet/PetParticleGun.hh"
#include "murat/pet/PetBrainPhantomGun.hh"
#include "murat/pet/PetGeomHandle.hh"
#include "murat/pet/BrainPhantom.hh"

#include "SeedService/inc/SeedService.hh"

// Includes from art and its toolchain.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Other external includes.
#include <boost/shared_ptr.hpp>

// C++ includes.
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

namespace mu2e {

  class PetEventGenerator : public art::EDProducer {

  public:

    explicit PetEventGenerator(fhicl::ParameterSet const& pSet);
    // Accept compiler written d'tor.  Modules are never moved or copied.

    virtual void produce (art::Event& e);
    virtual void beginRun(art::Run&   r);

  private:

    // Name of the run-time configuration file.
    string _configfile;

    // Control the behaviour of messages from the SimpleConfig object holding
    // the geometry parameters.
    bool _allowReplacement;
    bool _messageOnReplacement;
    bool _messageOnDefault;

    fhicl::ParameterSet _particleGun;
    fhicl::ParameterSet _conversionGun;
    fhicl::ParameterSet _brainPhantomGun;

    // Print final config file after all replacements.

    bool _printConfig;

    int  _configStatsVerbosity;
    int  _verbose;
    int  _printFrequency;


    // A collection of all of the generators that we will run.
    typedef  boost::shared_ptr<GeneratorBase> GeneratorBasePtr;
    std::vector<GeneratorBasePtr> _generators;

    void checkConfig( const SimpleConfig&  config);

  };

  PetEventGenerator::PetEventGenerator(fhicl::ParameterSet const& pSet):
    _configfile          (pSet.get<std::string>        ("inputfile"           ,"none")),
    _allowReplacement    (pSet.get<bool>               ("allowReplacement"    ,true )),
    _messageOnReplacement(pSet.get<bool>               ("messageOnReplacement",false)),
    _messageOnDefault    (pSet.get<bool>               ("messageOnDefault"    ,false)),
    _particleGun         (pSet.get<fhicl::ParameterSet>("particleGun"         ,fhicl::ParameterSet())),
    _conversionGun       (pSet.get<fhicl::ParameterSet>("conversionGun"       ,fhicl::ParameterSet())),
    _brainPhantomGun     (pSet.get<fhicl::ParameterSet>("brainPhantomGun"     ,fhicl::ParameterSet())),
    _printConfig         (pSet.get<bool>               ("printConfig"         ,false)),
    _configStatsVerbosity(pSet.get<int>                ("configStatsVerbosity",0    )),
    _verbose             (pSet.get<int>                ("verbose"             ,0    )),
    _printFrequency      (pSet.get<int>                ("printFrequency"      ,10   ))
 {

    produces<GenParticleCollection>();

    // A common random engine for the generators to use.
    createEngine( art::ServiceHandle<SeedService>()->getSeed() );
  }

//-----------------------------------------------------------------------------
// At beginRun time, update any derived geometry information.
//-----------------------------------------------------------------------------
  void PetEventGenerator::beginRun(art::Run &run) {
    //    double dzmax, rmax;

    static int ncalls(0);
    if ( ++ncalls > 1){
      mf::LogInfo("PetEventGenerator")
        << "PetEventGenerator does not change state at beginRun.  Hope that's OK.";
      return;
    }

    cout << "Event generator configuration file: " << _configfile << endl << endl;

    SimpleConfig config(_configfile, _allowReplacement, _messageOnReplacement, _messageOnDefault );

    if (_configfile != "none") {
      checkConfig(config);
      if (_printConfig) config.print(cout,"EvtGen: ");
    }

    // Change this to modify rather than delete and make an new one??

    bool use_particle_gun      = _particleGun.get    <bool> ("use",false);
    bool use_conversion_gun    = _conversionGun.get  <bool> ("use",false);
    bool use_brain_phantom_gun = _brainPhantomGun.get<bool> ("use",false);

    // Delete generators from the previous run.
    _generators.clear();

    // Instantiate generators for this run.
    if (use_particle_gun  ) _generators.push_back(GeneratorBasePtr(new PetParticleGun    (run,config)));
    if (use_conversion_gun) _generators.push_back(GeneratorBasePtr(new PetConversionGun  (run,config)));

    if (use_brain_phantom_gun) {
//-----------------------------------------------------------------------------
// by default, generate decays uniformly over the phantom volume 
// (defined by the geometry file)
// can redefine the dimensions talking to the PetEventGenerator
//-----------------------------------------------------------------------------
      PetBrainPhantomGun* gun = new PetBrainPhantomGun(run,_brainPhantomGun);
      _generators.push_back(GeneratorBasePtr(gun));
    }

    if ( _generators.size() == 0 ){
      mf::LogWarning("CONTROL")
        << "PetEventGenerator has no generators enabled. Hope that's OK.";
    }

    config.printAllSummaries( cout, _configStatsVerbosity, "EvtGen: ");
  }

  void PetEventGenerator::produce(art::Event& evt) {
    const char* oname = "PetEventGenerator::produce";

    if ((_verbose > 0) && (evt.event() % _printFrequency == 0)) {
      printf(" >>>>>>> [%s] RUN:EVENT : %10i:%10i\n",oname,evt.run(),evt.event());
    }
    
    // Make the collection to hold the output.
    unique_ptr<GenParticleCollection> genParticles(new GenParticleCollection);

    // Run all of the registered generators.
    for (std::vector<GeneratorBasePtr>::const_iterator i = _generators.begin(); i != _generators.end(); ++i) {
      (*i)->generate(*genParticles);
    }

    // Put the generated particles into the event.
    evt.put(std::move(genParticles));

  }

  // Look for inconsistencies in the config file.
  void PetEventGenerator::checkConfig( const SimpleConfig&  config){
  }


}


using mu2e::PetEventGenerator;
DEFINE_ART_MODULE(PetEventGenerator);
