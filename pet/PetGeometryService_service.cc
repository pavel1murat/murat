//
// Maintain up to date geometry information and serve it to
// other services and to the modules.
//
// $Id: PetGeometryService_service.cc,v 1.7 2013/11/20 18:15:52 murat Exp $
// $Author: murat $
// $Date: 2013/11/20 18:15:52 $
//
// Original author Rob Kutschke
//

// C++ include files
#include <iostream>
#include <typeinfo>

// Framework include files
#include "art/Persistency/Provenance/ModuleDescription.h"
#include "art/Persistency/Provenance/EventID.h"
#include "art/Persistency/Provenance/Timestamp.h"
#include "art/Persistency/Provenance/SubRunID.h"
#include "art/Persistency/Provenance/RunID.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


// Mu2e include files

#include "murat/pet/PetGeometryService.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/src/DetectorSystemMaker.hh"
#include "murat/pet/PetWorldG4.hh"
#include "murat/pet/PetWorldG4Maker.hh"

#include "murat/pet/BrainPhantomMaker.hh"
#include "murat/pet/BrainPhantom.hh"

#include "murat/pet/BrainImagerMaker.hh"
#include "murat/pet/BrainImager.hh"

#include "murat/pet/PetEnvelope.hh"

using namespace std;

namespace mu2e {

  PetGeometryService::PetGeometryService(fhicl::ParameterSet const& pset,
                                   art::ActivityRegistry&iRegistry) :
    _inputfile(            pset.get<std::string> ("inputFile",            "geom000.txt")),
    _allowReplacement(     pset.get<bool>        ("allowReplacement",     true)),
    _messageOnReplacement( pset.get<bool>        ("messageOnReplacement", false)),
    _messageOnDefault(     pset.get<bool>        ("messageOnDefault",     false)),
    _configStatsVerbosity( pset.get<int>         ("configStatsVerbosity", 0)),
    _printConfig(          pset.get<bool>        ("printConfig",          false)),
    _config(nullptr),
    _detectors(),
    _run_count()
  {
    iRegistry.sPreBeginRun.watch(this, &PetGeometryService::preBeginRun);
    iRegistry.sPostEndJob.watch (this, &PetGeometryService::postEndJob );
  }

  PetGeometryService::~PetGeometryService(){
  }

  // This template can be defined here because this is a private method which is only
  // used by the code below in the same file.
  template <typename DET>
  void PetGeometryService::addDetector(std::unique_ptr<DET> d)
  {
    if(_detectors.find(typeid(DET).name())!=_detectors.end())
      throw cet::exception("GEOM") << "failed to install detector with type name "
                                   << typeid(DET).name() << "\n";

      DetectorPtr ptr(d.release());
      _detectors[typeid(DET).name()] = ptr;
  }

  template <typename DETALIAS, typename DET> 
  void PetGeometryService::addDetectorAliasToBaseClass(std::unique_ptr<DET> d)
  {

	std::string OriginalName = typeid(DET).name();
	DetMap::iterator it(_detectors.find(OriginalName));

	if(it==_detectors.end())
          throw cet::exception("GEOM")
            << "Can not alias an inexistant detector, detector " << OriginalName << "\n";

	std::string detectorName= typeid(DETALIAS).name() ;
	_detectors[detectorName] = it->second;
  }

//-----------------------------------------------------------------------------
  void PetGeometryService::preBeginRun(art::Run const &) {

    if(++_run_count > 1) {
      mf::LogWarning("GEOM") << "This test version does not change geometry on run boundaries.";
      return;
    }

    cout  << "Geometry input file is: " << _inputfile << "\n";

    _config = unique_ptr<SimpleConfig>(new SimpleConfig(_inputfile,
                                                      _allowReplacement,
                                                      _messageOnReplacement,
                                                      _messageOnDefault ));

    // Print final state of file after all substitutions.
    if ( _printConfig      ){ _config->print(cout, "PetGeom: ");       }

    // decide if this is standard Mu2e detector or something else ...

    if (!_config->getBool("mu2e.standardDetector",true)) {
      cout  << "Non standard mu2e configuration, assuming it is intentional" << endl;
      return;
    }

    // Throw if the configuration is not self consistent.
    //    checkConfig();

    // This must be the first detector added since other makers may wish to use it.
    addDetector(DetectorSystemMaker::make( *_config));

    // Make a detector for every component present in the configuration.

    if (_config->getBool("hasBrainPhantom",false)) {
      BrainPhantomMaker brain_phantom_maker(*_config);
      addDetector(brain_phantom_maker.detectorPtr());
    }

    if (_config->getBool("hasBrainImager",false)) {
      BrainImagerMaker brain_imager_maker(*_config);
      addDetector(brain_imager_maker.detectorPtr());
    }

    addDetector(std::unique_ptr<PetEnvelope>(new PetEnvelope(-400.,4000.,-400.,400.,-900.,100)));

  } // preBeginRun()

  // WorldG4 could be added along with all the other detectors in preBeginRun().
  // However we don't want to make WorldG4 available in non-Geant jobs.
  // Therefore it is added by G4_module via this dedicated call.
  void PetGeometryService::addWorldG4() {
    addDetector(PetWorldG4Maker::make(*_config));
  }

  // Called after all modules have completed their end of job.
  void   PetGeometryService::postEndJob(){
    _config->printAllSummaries( cout, _configStatsVerbosity, "Geom: " );
  }



} // end namespace mu2e

DEFINE_ART_SERVICE(mu2e::PetGeometryService);
