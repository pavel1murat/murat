//
// G4 begin and end of event actions for Mu2e test environment
//
// $Id: GaasEventAction.cc,v 1.1 2012/11/16 23:53:04 genser Exp $
// $Author: genser $
// $Date: 2012/11/16 23:53:04 $
//
// Original author KLG
//

#include "murat/gaas/GaasEventAction.hh"
#include "murat/gaas/GaasSteppingAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"

namespace mu2e {

GaasEventAction::GaasEventAction(GaasSteppingAction *stepping_action) {
  _stepping = stepping_action;
}

GaasEventAction::~GaasEventAction()
{}


void GaasEventAction::BeginOfEventAction(const G4Event*) {
  //  _stepping->BeginOfEvent();
}


void GaasEventAction::EndOfEventAction(const G4Event* evt)
{
  _stepping->EndOfEvent();

  // Change  G4_plugin so that we copy to the art::event from here.

}

} // end namespace mu2e
