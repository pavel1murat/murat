#ifndef Mu2eG4_PetCaloCrystalSD_hh
#define Mu2eG4_PetCaloCrystalSD_hh
//
// Define a sensitive detector for calorimetric crystals
//
// $Id: PetCaloCrystalSD.hh,v 1.1 2013/06/06 23:02:21 murat Exp $
// $Author: murat $
// $Date: 2013/06/06 23:02:21 $
//
// Original author Ivan Logashenko
//

#include <map>
#include <vector>

// Mu2e includes
#include "Mu2eG4/inc/EventNumberList.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "murat/pet/PetSensitiveDetector.hh"

// Art includes
#include "art/Persistency/Provenance/ProductID.h"
#include "art/Framework/Principal/Event.h"

namespace mu2e {

  class PetCaloCrystalSD : public PetSensitiveDetector{

  public:

    PetCaloCrystalSD(G4String, const SimpleConfig& config);

    G4bool ProcessHits(G4Step*, G4TouchableHistory*);

  };

} // namespace mu2e

#endif /* Mu2eG4_PetCaloCrystalSD_hh */
