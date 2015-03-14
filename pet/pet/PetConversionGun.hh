#ifndef EventGenerator_PetConversionGun_hh
#define EventGenerator_PetConversionGun_hh

//
// Generate an electron with the conversion energy
// Uses FoilParticleGenerator to extract a random spot
// within the target system at
// a random time during the accelerator cycle.
//
// $Id: PetConversionGun.hh,v 1.1 2013/06/04 05:42:44 murat Exp $
// $Author: murat $
// $Date: 2013/06/04 05:42:44 $
//

// C++ includes
#include <memory>

// Mu2e includes
#include "EventGenerator/inc/FoilParticleGenerator.hh"
#include "EventGenerator/inc/GeneratorBase.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"

// Forward declarations in other namespaces.
namespace art {
  class Run;
}


// More forward declarations.
class TH1F;
class TH2F;

namespace mu2e {

  // Forward reference.
  class SimpleConfig;

  class PetConversionGun: public GeneratorBase{

  public:
    PetConversionGun( art::Run& run, const SimpleConfig& config );
    virtual ~PetConversionGun();

    virtual void generate( GenParticleCollection&  );

  private:

    // Conversion momentum.
    double _p;

    // Limits on the generated direction.
    double _czmin;
    double _czmax;
    double _phimin;
    double _phimax;

    bool _PStoDSDelay;
    bool _pPulseDelay;
    double _pPulseShift;

    // Activate the folding procedure on generation time. Default is on
    bool _timeFolding;

    // Limits on the generated time.
    double _tmin;
    double _tmax;

    // Select the position, type and time type for the generation
    std::string _foilGen;
    std::string _posGen;
    std::string _timeGen;

    // Control histograms.
    bool _doHistograms;

    //Utility to generate direction of the momentum, random on the unit sphere.
    RandomUnitSphere _randomUnitSphere;
    std::string _STfname;
    int _nToSkip;


    // Class object to generate position and time of the particle
    std::unique_ptr<FoilParticleGenerator> _fGenerator;

    //Particle mass
    double _mass;

    // Histograms.
    TH1F* _hMultiplicity;
    TH1F* _hcz;
    TH1F* _hphi;
    TH1F* _hmomentum;
    TH1F* _hradius;
    TH1F* _hzPos;
    TH1F* _htime;
    TH1F* _hmudelay;
    TH1F* _hpulsedelay;
    TH2F* _hxyPos;
    TH2F* _hrzPos;

    void bookHistograms();

  };

} // end namespace mu2e,

#endif /* EventGenerator_PetConversionGun_hh */


