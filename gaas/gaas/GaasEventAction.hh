#ifndef Mu2eG4_GaasEventAction_hh
#define Mu2eG4_GaasEventAction_hh
//
// G4 begin and end of event actions for Mu2e test environment
//
// $Id: GaasEventAction.hh,v 1.1 2012/11/16 23:29:03 genser Exp $
// $Author: genser $
// $Date: 2012/11/16 23:29:03 $
//
// Original author KLG
//

#include "G4UserEventAction.hh"

class G4Event;

namespace mu2e {

class GaasSteppingAction;

class GaasEventAction : public G4UserEventAction
{
  public:
    GaasEventAction(GaasSteppingAction*);
   ~GaasEventAction();

  public:
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);

  private:
    GaasSteppingAction * _stepping;
};

}  // end namespace mu2e
#endif /* Mu2eG4_GaasEventAction_hh */


