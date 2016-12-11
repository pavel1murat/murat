//
// - Serves as the link between ROOT "events" (e.g. mouse-clicks) and the ART
//   event display service by providing a receiver slot for signals generated
//   by the ROOT events.  A ROOT dictionary needs to be generated for this.
//


#include "murat/ana/NavState.hh"
#include "murat/ana/TEventDisplayUtils.hh"
#include <string>

ClassImp(TEventDisplayUtils)

//-----------------------------------------------------------------------------
TEventDisplayUtils::TEventDisplayUtils(): TObject() {
  fTbRun = 0;
  fTbEvt = 0;
}

//-----------------------------------------------------------------------------
TEventDisplayUtils::~TEventDisplayUtils() {
}

void TEventDisplayUtils::PrevEvent(){
  NavState::Set(kPREV_EVENT);
}

void TEventDisplayUtils::NextEvent(){
  NavState::Set(kNEXT_EVENT);
}

void TEventDisplayUtils::GotoEvent(){
  int run   = std::stoi(fTbRun->GetString());
  int event = std::stoi(fTbEvt->GetString());
  
  NavState::SetTarget(run, event);
  NavState::Set(kGOTO_EVENT);
}
