//
// A helper class manage repeated tasks related to sensitive detectors.
//
// $Id: PetSensitiveDetectorHelper.cc,v 1.3 2013/10/15 23:41:14 murat Exp $
// $Author: murat $
// $Date: 2013/10/15 23:41:14 $
//
// Original author Rob Kutschke
//
//
// Notes.
//
// 1) The timeVD virtual detector does not use the G4 sensitve detector methodology.
//    It is, instead, implemented inside SteppingAction.  So it is not managed
//    by this class.
//
// 2) The embedded class StepIntance needs to be copyable so that it can be put into
//    the std::vector.  Therefore it cannot contain an unique_ptr to its StepPointMC
//    collection.  The work around is to hold the collection by value and to use swap
//    to transfer it into the unique_ptr that will be given to the event.  This is
//    a very small CPU time penalty but it saves us from doing any explicit memory management.
//

// From Mu2e.
#include "murat/pet/PetSensitiveDetectorHelper.hh"
#include "MCDataProducts/inc/StatusG4.hh"

#include "murat/pet/PetSensitiveDetectorName.hh"
#include "murat/pet/PetG4Helper.hh"

#include "Mu2eG4/inc/SimParticleHelper.hh"

// From art and its tool chain.
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

// From G4.
#include "G4VSensitiveDetector.hh"
#include "G4SDManager.hh"

#include <map>

using namespace std;

namespace mu2e {

  // pset argument reserved for future use; see note 3.
  PetSensitiveDetectorHelper::PetSensitiveDetectorHelper( fhicl::ParameterSet const& pset ):
    stepInstances_(){

    // Backward compatible defaults
    bool enableAllSDs = true;
    std::vector<std::string> enabledSDs;

    if(!pset.is_empty()) {
      enableAllSDs = pset.get<bool>("enableAllSDs", false);
      const bool explicitList = pset.get_if_present<std::vector<std::string> >("enableSD", enabledSDs);
      if(enableAllSDs && explicitList) {
        throw cet::exception("CONFIG")<<"SensitiveDetectorHelper: enableAllSDs=true conflicts with explicit list enableSD\n";
      }
      if(!enableAllSDs && !explicitList) {
        throw cet::exception("CONFIG")<<"SensitiveDetectorHelper: no SDs are defined.  Either say enableAllSDs:true, or enableSD:[list]\n";
      }
    }

    typedef std::set<std::string> TODO;
    TODO todo(enabledSDs.begin(), enabledSDs.end());

    // Build list of StepInstances to look after.  See note 1.
    std::vector<StepInstanceName> const& preDefinedSD(StepInstanceName::allValues());
    for ( std::vector<StepInstanceName>::const_iterator i=preDefinedSD.begin(); i != preDefinedSD.end(); ++i ){
      if(enableAllSDs || (todo.find(i->name()) != todo.end())) {
        todo.erase(i->name());

        // stepper does not have a sensitive detector associated with it...
        if (*i == StepInstanceName::stepper) continue;
        StepInstanceName::enum_type id(*i);
        if ( id != StepInstanceName::unknown && id !=StepInstanceName::timeVD ){
          stepInstances_[id] = StepInstance(id);
        }
      }
    }

    if(!todo.empty()) {
      std::ostringstream os;
      os<<"SensitiveDetectorHelper: the following detectors in the enableSD list are not recognized: ";
      std::copy(todo.begin(), todo.end(), std::ostream_iterator<std::string>(os, " "));
      throw cet::exception("CONFIG")<<os.str();
    }

    //----------------
    std::vector<string> lvlist(pset.get<vector<string>>("sensitiveVolumes", {}));
    for(const auto& name : lvlist) {
      lvsd_[name] = StepInstance(name);
    }

    //----------------
    // Careful here: we are in the constructor, so class methods may not be used freely.
    // However calling the following at OK at this point (but not earlier!)

    const vector<string> outputs(stepInstanceNamesToBeProduced());

    //----------------
    // New hits will be added to hits from these existing collections
    for(const auto& s : pset.get<vector<string> >("preSimulatedHits", vector<string>())) {
      preSimulatedHits_.emplace_back(s);
      // check that the destination is known and enabled
      if(std::find(outputs.begin(), outputs.end(), preSimulatedHits_.back().instance()) == outputs.end()) {
        throw cet::exception("CONFIG")<<"PetSensitiveDetectorHelper: no matching destination for preSimulatedHits = "<<s<<"\n";
      }
    }

    //----------------
  }

  void PetSensitiveDetectorHelper::instantiateLVSDs(const SimpleConfig& config){
    G4SDManager* SDman      = G4SDManager::GetSDMpointer();
    art::ServiceHandle<PetG4Helper> helper;

    for(auto& iter : lvsd_) {
      iter.second.sensitiveDetector = new PetSensitiveDetector(iter.first, config);
      SDman->AddNewDetector(iter.second.sensitiveDetector);
      helper->locateVolInfo(iter.first).logical->SetSensitiveDetector(iter.second.sensitiveDetector);
    }
  }

  // Find the sensitive detector objects and attach them to each StepInstance object.
  // Must not be called until G4 has been initialized.
  void PetSensitiveDetectorHelper::registerSensitiveDetectors(){
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();
    for ( InstanceMap::iterator i=stepInstances_.begin();
            i != stepInstances_.end(); ++i ){
      StepInstance& step(i->second);
      step.sensitiveDetector = 
        dynamic_cast<PetSensitiveDetector*>(sdManager->FindSensitiveDetector(step.stepName.c_str()));
    }
  }

  // Create new data products.  To be called at start of each event.
  void PetSensitiveDetectorHelper::createProducts(const art::Event& event, const SimParticleHelper& spHelper){

    //----------------
    // Read in pre-simulated hits that we want to merge into the outputs

    // instance => multiple hit collections
    typedef std::multimap<string, art::ValidHandle<StepPointMCCollection> > InHits;
    InHits inputHits;
    for(const auto& tag : preSimulatedHits_) {
      inputHits.insert(make_pair(tag.instance(), event.getValidHandle<StepPointMCCollection>(tag)));
    }

    //----------------
    // Clean and pre-fill pre-defined SD collections in stepInstances_

    for (auto& i : stepInstances_) {
      auto& out = i.second.p;
      out.clear();

      // Copy all input collection with the current instance name
      const auto rr = inputHits.equal_range(i.second.stepName);
      for(auto in = rr.first; in != rr.second; ++in) {
        out.insert(out.end(), in->second->cbegin(), in->second->cend());
      }

      // Update SimParticle ptr to point to the new collection
      for(auto& hit : out) {
        hit.simParticle() =
          art::Ptr<SimParticle>(spHelper.productID(),
                                hit.simParticle()->id().asUint(),
                                spHelper.productGetter());
      }
    }

    //----------------
    // Clean and pre-fill the logical volume collections

    for(auto& i : lvsd_) {
      auto& out = i.second.p;
      out.clear();

      // Copy all input collection with the current instance name
      const auto rr = inputHits.equal_range(i.first);
      for(auto in = rr.first; in != rr.second; ++in) {
        out.insert(out.end(), in->second->cbegin(), in->second->cend());
      }

      // Update SimParticle ptr to point to the new collection
      for(auto& hit : out) {
        hit.simParticle() =
          art::Ptr<SimParticle>(spHelper.productID(),
                                hit.simParticle()->id().asUint(),
                                spHelper.productGetter());
      }
    }
  }

  // Hand each sensitive detector object the StepPointMC collection that it will fill.
  // Also hand it the information needed to create art::Ptr objects.
  void PetSensitiveDetectorHelper::updateSensitiveDetectors( PhysicsProcessInfo&   info,
                                                          const SimParticleHelper& spHelper){

    for ( InstanceMap::iterator i=stepInstances_.begin();
          i != stepInstances_.end(); ++i ){
      StepInstance& instance(i->second);
      if ( instance.sensitiveDetector ){
        instance.sensitiveDetector->beforeG4Event(instance.p, info, spHelper);
      }
    }

    for(auto& i : lvsd_) {
      i.second.sensitiveDetector->beforeG4Event(i.second.p, info, spHelper);
    }

  }

  // Put all of the data products into the event.
  // To be called at the end of each event.
  void PetSensitiveDetectorHelper::put( art::Event& event ){

    for ( InstanceMap::iterator i=stepInstances_.begin();
            i != stepInstances_.end(); ++i ){
      unique_ptr<StepPointMCCollection> p(new StepPointMCCollection);
      StepInstance& instance(i->second);
      std::swap( instance.p, *p);
      event.put(std::move(p), instance.stepName );
    }

    for (auto& i: lvsd_) {
      unique_ptr<StepPointMCCollection> p(new StepPointMCCollection);
      std::swap( i.second.p, *p);
      event.put(std::move(p), i.second.stepName );
    }

  }

  // Return all of the instances names of the data products to be produced.
  vector<string> PetSensitiveDetectorHelper::stepInstanceNamesToBeProduced() const{
    vector<string> names;
    for ( InstanceMap::const_iterator i=stepInstances_.begin();
            i != stepInstances_.end(); ++i ){
      names.push_back( i->second.stepName );
    }

    for(const auto& i : lvsd_) {
      names.emplace_back(i.second.stepName);
    }

    return names;
  }

  bool PetSensitiveDetectorHelper::enabled(StepInstanceName::enum_type instance) const{
    return stepInstances_.find(instance) != stepInstances_.end();
  }

  // Return one of the StepPointMCCollections.
  cet::maybe_ref<StepPointMCCollection> PetSensitiveDetectorHelper::steps( StepInstanceName::enum_type id ){
    InstanceMap::iterator i = stepInstances_.find(id);
    if ( i != stepInstances_.end() ) {
      return cet::maybe_ref<StepPointMCCollection>(i->second.p);
    }

    // Cannot find the requested id; return an empty maybe_ref.
    return cet::maybe_ref<StepPointMCCollection>();
  }

}  // end namespace mu2e
