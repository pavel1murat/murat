//
// A Producer Module that runs Geant4 and adds its output to the event.
// Still under development.
//
// $Id: PetG4_module.cc,v 1.7 2014/02/02 02:15:53 murat Exp $
// $Author: murat $
// $Date: 2014/02/02 02:15:53 $
//
// Original author Rob Kutschke
//
//
// Notes:
// 1) According to Sunanda Banerjee, the various SetUserAction methods
//    take ownership of the object that is passed to it.  So we must
//    not delete them.
//

// Mu2e includes
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "murat/pet/PetG4RunManager.hh"
#include "murat/pet/PetWorldMaker.hh"
#include "murat/pet/PetWorld.hh"

#include "Mu2eG4/inc/addPointTrajectories.hh"
#include "Mu2eG4/inc/exportG4PDT.hh"

#include "murat/pet/PetGeometryService.hh"

#include "murat/pet/PetGeomHandle.hh"
#include "murat/pet/PetWorldG4.hh"

#include "Mu2eBuildingGeom/inc/Mu2eBuilding.hh"
#include "Mu2eG4/inc/DetectorConstruction.hh"

#include "murat/pet/PetPrimaryGeneratorAction.hh"
#include "murat/pet/PetStackingAction.hh"

#include "Mu2eG4/inc/EventAction.hh"
#include "Mu2eG4/inc/SteppingAction.hh"
#include "Mu2eG4/inc/SteppingVerbose.hh"

#include "Mu2eG4/inc/TrackingAction.hh"
#include "Mu2eG4/inc/PhysicalVolumeHelper.hh"
#include "Mu2eG4/inc/PhysicsProcessInfo.hh"
#include "Mu2eG4/inc/physicsListDecider.hh"
#include "Mu2eG4/inc/postG4InitializeTasks.hh"

#include "murat/pet/PetSensitiveDetector.hh"
#include "murat/pet/PetSensitiveDetectorName.hh"
#include "murat/pet/PetSensitiveDetectorHelper.hh"

#include "murat/pet/DiagnosticsPetG4.hh"

#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Mu2eG4/inc/SimParticleHelper.hh"
#include "Mu2eG4/inc/SimParticlePrimaryHelper.hh"
#include "SeedService/inc/SeedService.hh"
#include "Mu2eUtilities/inc/SimParticleCollectionPrinter.hh"

// Data products that will be produced by this module.
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoMultiCollection.hh"
#include "MCDataProducts/inc/PointTrajectoryCollection.hh"
#include "MCDataProducts/inc/StatusG4.hh"
#include "MCDataProducts/inc/StepInstanceName.hh"
#include "MCDataProducts/inc/MCTrajectoryCollection.hh"

// From art and its tool chain.
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Utilities/InputTag.h"

// Geant4 includes
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4Run.hh"
#include "G4Timer.hh"
#include "G4VUserPhysicsList.hh"
#include "G4RunManagerKernel.hh"
#include "G4SDManager.hh"

// ROOT includes
#include "TNtuple.h"

// C++ includes.
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <memory>
#include <iomanip>
#include <utility>

using namespace std;

namespace mu2e {

  class PetG4 : public art::EDProducer {

  public:
    PetG4(fhicl::ParameterSet const& pSet);
    // Accept compiler supplied d'tor

    virtual void produce(art::Event& e);

    virtual void beginJob();
    virtual void endJob();

    virtual void beginRun(art::Run &r);
    virtual void endRun(art::Run &);

    virtual void beginSubRun(art::SubRun &sr);

  private:
    typedef std::vector<art::InputTag> InputTags;
    typedef std::vector<std::string> Strings;

    unique_ptr<PetG4RunManager> _runManager;

    // Do we issue warnings about multiple runs?
    bool _warnEveryNewRun;

    // Do we want to export the G4 particle data table.
    bool  _exportPDTStart;
    bool  _exportPDTEnd;

    PetPrimaryGeneratorAction* _genAction;
    TrackingAction*         _trackingAction;
    SteppingAction*         _steppingAction;
    PetStackingAction*         _stackingAction;

    G4UIsession  *_session;
    G4UImanager  *_UI;
    std::unique_ptr<G4VisManager> _visManager;
    int _rmvlevel;
    int _tmvlevel;
    int _checkFieldMap;

    // Name of a macro file for visualization.
    string _visMacro;

    // Name of a macro file to be used for controling G4 parameters after
    // the initialization phase.
    string _g4Macro;

    art::InputTag _generatorModuleLabel;
    InputTags _genInputHitLabels;

    string _inputPhysVolumeMultiInfoLabel;
    bool   _doWriteLegacyPhysVolumeInfo;

    // Helps with indexology related to persisting info about G4 volumes.
    PhysicalVolumeHelper _physVolHelper;

    // Helps with recording information about physics processes.
    PhysicsProcessInfo _processInfo;
    bool _printPhysicsProcessSummary;

    SimParticleCollectionPrinter _simParticlePrinter;

    PetSensitiveDetectorHelper _sensitiveDetectorHelper;
    //    ExtMonFNALPixelSD       *_extMonFNALPixelSD;

    // Instance name of the timeVD StepPointMC data product.
    const StepInstanceName _tvdOutputName;

    unsigned _simParticleNumberOffset;
    art::InputTag _inputSimParticles;

    // A class to make some standard histograms.
    DiagnosticsPetG4 _diagnostics;

    // A parameter extracted from the geometry file at beginRun and used in produce.
    int _pointTrajectoryMinSteps;

    // Do the PetG4 initialization that must be done only once per job, not once per run
    void initializeG4( PetGeometryService& geom, art::Run const& run );

  }; // end PetG4 header

  PetG4::PetG4(fhicl::ParameterSet const& pSet):
    _runManager(nullptr),
    _warnEveryNewRun(pSet.get<bool>("warnEveryNewRun",false)),
    _exportPDTStart(pSet.get<bool>("exportPDTStart",false)),
    _exportPDTEnd(pSet.get<bool>("exportPDTEnd",false)),
    _genAction(nullptr),
    _trackingAction(nullptr),
    _steppingAction(nullptr),
    _stackingAction(nullptr),
    _session(nullptr),
    _UI(nullptr),
    _visManager(nullptr),
    _rmvlevel(pSet.get<int>("diagLevel",0)),
    _tmvlevel(pSet.get<int>("trackingVerbosityLevel",0)),
    _checkFieldMap(pSet.get<int>("checkFieldMap",0)),
    _visMacro(pSet.get<std::string>("visMacro","")),
    _g4Macro(pSet.get<std::string>("g4Macro","")),
    _generatorModuleLabel(pSet.get<std::string>("generatorModuleLabel", "")),
    _inputPhysVolumeMultiInfoLabel(pSet.get<string>("inputPhysVolumeMultiInfoLabel", "")),
    _doWriteLegacyPhysVolumeInfo(pSet.get<bool>("doWriteLegacyPhysVolumeInfo", true)),
    _physVolHelper(),
    _processInfo(),
    _printPhysicsProcessSummary(false),
    _simParticlePrinter(pSet.get<fhicl::ParameterSet>("SimParticlePrinter", SimParticleCollectionPrinter::defaultPSet())),
    _sensitiveDetectorHelper(pSet.get<fhicl::ParameterSet>("SDConfig", fhicl::ParameterSet())),
  //    _extMonFNALPixelSD(),
    _tvdOutputName(StepInstanceName::timeVD),
    _simParticleNumberOffset(pSet.get<unsigned>("simParticleNumberOffset", 0)),
    _inputSimParticles(pSet.get<std::string>("inputSimParticles", "")),
    _diagnostics(),
    _pointTrajectoryMinSteps(-1){

    Strings genHitsStr(pSet.get<Strings>("genInputHits", Strings()));
    for(const auto& s : genHitsStr) {
      _genInputHitLabels.emplace_back(s);
    }

    if((_generatorModuleLabel == art::InputTag()) && _genInputHitLabels.empty()) {
      throw cet::exception("CONFIG")
        << "Error: both generatorModuleLabel and genInputHits are empty - nothing to do!\n";
    }

    produces<StatusG4>();
    produces<SimParticleCollection>();

    // The main group of StepPointMCCollections.
    vector<string> const& instanceNames = _sensitiveDetectorHelper.stepInstanceNamesToBeProduced();
    for ( vector<string>::const_iterator i=instanceNames.begin();
          i != instanceNames.end(); ++i){
      produces<StepPointMCCollection>(*i);
    }

    // The timevd collection is special.
    produces<StepPointMCCollection>(_tvdOutputName.name());

    produces<PointTrajectoryCollection>();
    produces<PhysicalVolumeInfoCollection,art::InRun>();
    produces<PhysicalVolumeInfoMultiCollection,art::InSubRun>();

    // The string "G4Engine" is magic; see the docs for RandomNumberGenerator.
    createEngine( art::ServiceHandle<SeedService>()->getSeed(), "G4Engine");

  } // end PetG4:PetG4(fhicl::ParameterSet const& pSet);

  // Create an instance of the run manager.
  void PetG4::beginJob(){
    _runManager = unique_ptr<PetG4RunManager>(new PetG4RunManager);
  }

//-----------------------------------------------------------------------------
  void PetG4::beginRun( art::Run &run){

    static int ncalls(0);
    ++ncalls;

    art::ServiceHandle<PetGeometryService> geom;

    // Do the main initialization of PetG4; only once per job.
    if ( ncalls == 1 ) {
      initializeG4(*geom,run);
    } 
    else {
      if ( ncalls ==2 || _warnEveryNewRun ){
        mf::LogWarning log("G4");
        log << "G4 does not change state when we cross run boundaries - hope this is OK .... ";
        if ( ncalls == 2 && !_warnEveryNewRun ){
          log << "\nThis message will not be repeated on subsequent new runs.";
        }
      }
    }

    // Tell G4 that we are starting a new run.
    _runManager->BeamOnBeginRun( run.id().run() );

    // Helps with indexology related to persisting G4 volume information.
    _physVolHelper.beginRun();
    _processInfo.beginRun();

    if(_doWriteLegacyPhysVolumeInfo) {
      // Add info about the G4 volumes to the run-data.
      // The framework rules requires we make a copy and add the copy.
      const PhysicalVolumeInfoCollection& vinfo = _physVolHelper.persistentInfo();
      unique_ptr<PhysicalVolumeInfoCollection> volumes(new PhysicalVolumeInfoCollection(vinfo));
      run.put(std::move(volumes));
    }

    // Some of the user actions have beginRun methods.
    PetGeomHandle<PetWorldG4>  worldGeom;
    _trackingAction->beginRun( _physVolHelper, _processInfo, worldGeom->mu2eOriginInWorld() );
    _steppingAction->beginRun( _processInfo, worldGeom->mu2eOriginInWorld() );
    _stackingAction->beginRun( -1.e6, 1.e6);

    // A few more things that only need to be done only once per job,
    // not once per run, but which need to be done after the call to
    // BeamOnBeginRun.

    if ( ncalls == 1 ) {

      _steppingAction->finishConstruction();

      //      if( _checkFieldMap>0 ) generateFieldMap(worldGeom->mu2eOriginInWorld(),_checkFieldMap);

      if ( _exportPDTStart ) exportG4PDT( "Start:" );
    }

    // Get some run-time configuration information that is stored in the geometry file.
    SimpleConfig const& config  = geom->config();
    _printPhysicsProcessSummary = config.getBool("g4.printPhysicsProcessSummary",false);

    _pointTrajectoryMinSteps = config.getInt("g4.pointTrajectoryMinSteps",5);
  }

  void PetG4::initializeG4( PetGeometryService& geom, art::Run const& run ){

    SimpleConfig const& config = geom.config();

    geom.addWorldG4();

    if (_rmvlevel > 0) {
      mf::LogInfo logInfo("GEOM");
      logInfo << "Initializing Geant 4 for " << run.id()
              << " with verbosity " << _rmvlevel << endl;
      logInfo << "Configured simParticleNumberOffset = "<< _simParticleNumberOffset << endl;
    }

    // Create user actions and register them with G4.

    PetWorldMaker<PetWorld>* allMu2e;

    allMu2e = new PetWorldMaker<PetWorld>(std::unique_ptr<PetWorld>(new PetWorld(&_sensitiveDetectorHelper)));

    _runManager->SetVerboseLevel(_rmvlevel);

    _runManager->SetUserInitialization(allMu2e);

    G4VUserPhysicsList* pL = physicsListDecider(config);
    pL->SetVerboseLevel(_rmvlevel);

    _runManager->SetUserInitialization(pL);

    _genAction = new PetPrimaryGeneratorAction();
    _runManager->SetUserAction(_genAction);

    _steppingAction = new SteppingAction(config);
    _runManager->SetUserAction(_steppingAction);

    G4UserEventAction* event_action = new EventAction(_steppingAction);
    _runManager->SetUserAction(event_action);

    _stackingAction = new PetStackingAction(config);
    _runManager->SetUserAction(_stackingAction);

    _trackingAction = new TrackingAction(config,_steppingAction);
    _runManager->SetUserAction(_trackingAction);

    // setting tracking/stepping verbosity level; tracking manager
    // sets stepping verbosity level as well; 

    G4RunManagerKernel const * rmk = G4RunManagerKernel::GetRunManagerKernel();
    G4TrackingManager* tm  = rmk->GetTrackingManager();
    tm->SetVerboseLevel(_tmvlevel);

    _UI = G4UImanager::GetUIpointer();

    // Any final G4 interactive commands ...
    if ( !_g4Macro.empty() ) {
      G4String command("/control/execute ");
      ConfigFileLookupPolicy path;
      command += path(_g4Macro);
      _UI->ApplyCommand(command);

    }

    // Initialize G4 for this run.
    _runManager->Initialize();

    // At this point G4 geometry and physics processes have been initialized.
    // So it is safe to modify physics processes and to compute information
    // that is derived from the G4 geometry or physics processes.

    // Mu2e specific customizations that must be done after the call to Initialize.
    postG4InitializeTasks(config);
    _sensitiveDetectorHelper.registerSensitiveDetectors();

    // Setup the graphics if requested.
    if ( !_visMacro.empty() ) {

      _visManager = std::unique_ptr<G4VisManager>(new G4VisExecutive);
      _visManager->Initialize();

      ConfigFileLookupPolicy visPath;

      G4String command("/control/execute ");
      command += visPath(_visMacro);

      _UI->ApplyCommand( command );
    }

    // Book some diagnostic histograms.
    art::ServiceHandle<art::TFileService> tfs;
    _diagnostics.book("Outputs");

  } // end PetG4::initializeG4


//-----------------------------------------------------------------------------
  void PetG4::beginSubRun(art::SubRun& sr) {
    unique_ptr<PhysicalVolumeInfoMultiCollection> mvi(new PhysicalVolumeInfoMultiCollection());

    if(!_inputPhysVolumeMultiInfoLabel.empty()) {
      // Copy over data from the previous simulation stages
      art::Handle<PhysicalVolumeInfoMultiCollection> ih;
      sr.getByLabel(_inputPhysVolumeMultiInfoLabel, ih);
      mvi->reserve(1 + ih->size());
      mvi->insert(mvi->begin(), ih->cbegin(), ih->cend());
    }

    // Append info for the current stage
    mvi->emplace_back(std::make_pair(_simParticleNumberOffset, _physVolHelper.persistentSingleStageInfo()));

    sr.put(std::move(mvi));
  }

//-----------------------------------------------------------------------------
// Create one G4 event and copy its output to the art::event.
//-----------------------------------------------------------------------------
  void PetG4::produce(art::Event& event) {

    // Handle to the generated particles; need when building art::Ptr to a GenParticle.
    art::Handle<GenParticleCollection> gensHandle;
    if(!(_generatorModuleLabel == art::InputTag())) {
      event.getByLabel(_generatorModuleLabel, gensHandle);
    }

    // input hits from the previous simulation stage
    HitHandles genInputHits;
    for(const auto& i : _genInputHitLabels) {
      genInputHits.emplace_back(event.getValidHandle<StepPointMCCollection>(i));
    }

    art::Handle<SimParticleCollection> inputSimHandle;
    if(!(art::InputTag() == _inputSimParticles)) {
      event.getByLabel(_inputSimParticles, inputSimHandle);
    }

    // ProductID for the SimParticleCollection.
    art::ProductID simPartId(getProductID<SimParticleCollection>(event));
    SimParticleHelper spHelper(_simParticleNumberOffset, simPartId, event);
    SimParticlePrimaryHelper parentHelper(event, simPartId, gensHandle);

    // Create empty data products.
    unique_ptr<SimParticleCollection>     simParticles(      new SimParticleCollection);
    unique_ptr<StepPointMCCollection>     tvdHits(           new StepPointMCCollection);
    unique_ptr<PointTrajectoryCollection> pointTrajectories( new PointTrajectoryCollection);
    unique_ptr<MCTrajectoryCollection>     mcTrajectories(    new MCTrajectoryCollection);
    _sensitiveDetectorHelper.createProducts(event, spHelper);

    // Some of the user actions have begin event methods. These are not G4 standards.
    _trackingAction->beginEvent(inputSimHandle, spHelper, parentHelper, *mcTrajectories );
    _genAction->setEventData(gensHandle.isValid() ? &*gensHandle : 0, genInputHits, &parentHelper);
    _steppingAction->BeginOfEvent(*tvdHits,  spHelper);

    // Connect the newly created StepPointMCCollections to their sensitive detector objects.
    _sensitiveDetectorHelper.updateSensitiveDetectors( _processInfo, spHelper);

    // Run G4 for this event and access the completed event.
    _runManager->BeamOnDoOneEvent( event.id().event() );
    G4Event const* g4event = _runManager->getCurrentEvent();

    // Populate the output data products.
    PetGeomHandle<PetWorldG4>  world;
    //    GeomHandle<Mu2eBuilding>  building;
    addPointTrajectories( g4event, *pointTrajectories, spHelper, world->mu2eOriginInWorld(), _pointTrajectoryMinSteps);

    // Run self consistency checks if enabled.
    _trackingAction->endEvent(*simParticles);

    // Fill the status object.
    G4Timer const* timer = _runManager->getG4Timer();
    float cpuTime  = timer->GetSystemElapsed()+timer->GetUserElapsed();

    int status(0);
    if (  _steppingAction->nKilledStepLimit() > 0 ) status =  1;
    if (  _trackingAction->overflowSimParticles() ) status = 10;

    unique_ptr<StatusG4> g4stat(new StatusG4( status,
                                            _trackingAction->nG4Tracks(),
                                            _trackingAction->overflowSimParticles(),
                                            _steppingAction->nKilledStepLimit(),
                                            cpuTime,
                                            timer->GetRealElapsed() )
                              );

    _diagnostics.fill( &*g4stat,
                       &*simParticles,

                       _sensitiveDetectorHelper.steps(StepInstanceName::calorimeter) ? 
                       &_sensitiveDetectorHelper.steps(StepInstanceName::calorimeter).ref() : nullptr,

                       _sensitiveDetectorHelper.steps(StepInstanceName::calorimeterRO) ?
                       &_sensitiveDetectorHelper.steps(StepInstanceName::calorimeterRO).ref() : nullptr,

                       _sensitiveDetectorHelper.steps(StepInstanceName::virtualdetector) ?
                       &_sensitiveDetectorHelper.steps(StepInstanceName::virtualdetector).ref() : nullptr,

                       &*pointTrajectories,
                       &_physVolHelper.persistentInfo() );

    _simParticlePrinter.print(std::cout, *simParticles);

    // Add data products to the event.
    event.put(std::move(g4stat));
    event.put(std::move(simParticles));
    event.put(std::move(tvdHits),          _tvdOutputName.name()          );
    event.put(std::move(pointTrajectories));
    _sensitiveDetectorHelper.put(event);

    // Pause to see graphics.
    if ( !_visMacro.empty() ){

      // Prompt to continue and wait for reply.
      cout << "Enter a character to go to the next event (q quits, v enters G4 interactive session)" <<
        endl;
      cout << "(Once in G4 interactive session to quit it type exit): ";
      string userinput;
      cin >> userinput;
      G4cout << userinput << G4endl;

      // Check if user is requesting an early termination of the event loop.
      if ( !userinput.empty() ){
        // Checks only the first character; we should check first non-blank.
        char c = tolower( userinput[0] );
        if ( c == 'q' ){
          throw cet::exception("CONTROL")
            << "Early end of event loop requested inside G4, \n";
        } else if ( c == 'v' ){
          G4int argc=1;
          // Cast away const-ness; required by the G4 interface ...
          char* dummy = (char *)"dummy";
          char** argv = &dummy;
          G4UIExecutive* UIE = new G4UIExecutive(argc, argv);
          UIE->SessionStart();
          delete UIE;
        }
      } // end !userinput.empty()

    }   // end !_visMacro.empty()

    // This deletes the object pointed to by currentEvent.
    _runManager->BeamOnEndEvent();

  }

  // Tell G4 that this run is over.
  void PetG4::endRun(art::Run & run){
    _runManager->BeamOnEndRun();
  }

  void PetG4::endJob(){

    if ( _exportPDTEnd ) exportG4PDT( "End:" );

    // Yes, these are named endRun, but they are really endJob actions.
    _physVolHelper.endRun();
    _trackingAction->endRun();

    if ( _printPhysicsProcessSummary ){
      _processInfo.endRun();
    }

  }

} // End of namespace mu2e

using mu2e::PetG4;
DEFINE_ART_MODULE(PetG4);
