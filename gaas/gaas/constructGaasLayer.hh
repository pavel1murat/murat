#ifndef Mu2eG4_constructGaasLayer_hh
#define Mu2eG4_constructGaasLayer_hh
//
// Free function to create one of the study environment geometries
//
// $Id: constructGaasEnv_v002.hh,v 1.2 2013/04/05 19:51:34 genser Exp $
// $Author: genser $
// $Date: 2013/04/05 19:51:34 $
//
// Original author KLG
//

namespace mu2e {

  class VolumeInfo;
  class SimpleConfig;

  void constructGaasLayer(VolumeInfo const & parentVInfo, SimpleConfig const & _config);

}

#endif /* Mu2eG4_constructGaasLayer_hh */
