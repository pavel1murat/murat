#ifndef CalorimeterGeom_CalorimeterMaker_hh
#define CalorimeterGeom_BrainImagerMaker_hh
// $Id: BrainImagerMaker.hh,v 1.3 2013/11/04 21:09:32 murat Exp $
// $Author: murat $
// $Date: 2013/11/04 21:09:32 $

// original authors Julie Managan and Robert Bernstein

// C++ includes
#include <iomanip>
#include <memory>
#include <iostream>
#include <cmath>
#include <vector>
#include <string>

//Mu2e includes
#include "murat/pet/BrainImager.hh"
#include "CalorimeterGeom/inc/Vane.hh"
#include "ConfigTools/inc/SimpleConfig.hh"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"


//forward declarations
class BrainImager;

class BrainImagerMaker{

public:
     
  BrainImagerMaker(mu2e::SimpleConfig const& config);
  ~BrainImagerMaker();
      
  // Accessor and unique_ptr to BrainImager needed by GeometryService.

  std::unique_ptr<BrainImager> _calo;
  std::unique_ptr<BrainImager> detectorPtr() { return std::move(_calo); }
  
private:
  
  void CheckIt  (void);
  void MakeVanes(void);

};

#endif /* CalorimeterGeom_BrainImagerMaker_hh */
