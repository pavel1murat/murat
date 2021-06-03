///////////////////////////////////////////////////////////////////////////////
// initial version cloned from Rob's Analyses/src/ReadMCTrajectories_module.cc
///////////////////////////////////////////////////////////////////////////////
// C++ includes.
#include <iostream>
#include <string>

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"

// Mu2e includes.
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/MCTrajectoryCollection.hh"

using namespace std;

namespace mu2e {

//-----------------------------------------------------------------------------
class DumpMCTrajectories : public art::EDAnalyzer {
public:

  explicit DumpMCTrajectories(fhicl::ParameterSet const& pset);
  
  virtual void beginJob();
  void analyze(const art::Event& e);

private:

  art::InputTag _simpCollTag;
  art::InputTag _mctrCollTag;
};

//-----------------------------------------------------------------------------
DumpMCTrajectories::DumpMCTrajectories(fhicl::ParameterSet const& pset) :
  art::EDAnalyzer(pset),
  _simpCollTag(pset.get<std::string>("simpCollTag")),
  _mctrCollTag(pset.get<std::string>("mctrCollTag"))
{
}

//-----------------------------------------------------------------------------
void DumpMCTrajectories::beginJob(){
}

//-----------------------------------------------------------------------------
void DumpMCTrajectories::analyze(const art::Event& event) {

  auto simpColl     = event.getValidHandle<SimParticleCollection>(_simpCollTag);

  auto trajectories = event.getValidHandle<MCTrajectoryCollection>(_mctrCollTag);

  char fn[10000];

  const mu2e::SimParticle* sim = &simpColl->begin()->second;
  

  for (auto const& itr : *trajectories ) {
    MCTrajectory const& traj        = itr.second;

    sprintf(fn,"%4i_%06i_%06i.txt",event.run(),event.subRun(),event.event());
//-----------------------------------------------------------------------------
// print parameters of the first simParticle , so far, handle only one trajectory
//-----------------------------------------------------------------------------
    FILE* f = fopen(fn,"w");
    int index = 0;
      
    fprintf(f,"    0 %10.3f %10.3f %10.3f %10.3f %10.3f\n",
	    sim->startPosition().x(),
	    sim->startPosition().y(),
	    sim->startPosition().z(),
	    sim->startGlobalTime(),
	    sim->startMomentum().e()-sim->startMomentum().m());
      
    for ( auto const& pt : traj.points() ){
      fprintf(f,"%5i %10.3f %10.3f %10.3f %10.3f %10.3f\n",index,pt.x(),pt.y(),pt.z(),pt.t(),pt.kineticEnergy());
      index++;
    }
      
    fclose(f);
    break;
  }
}
  
}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::DumpMCTrajectories);
