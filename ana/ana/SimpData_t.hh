#ifndef __murat_ana_SimpData__
#define __murat_ana_SimpData__

namespace murat {
  struct SimpData_t {
    int           fIndex;		// so far, not used
    TSimParticle* fParent;              // parent (in production vertex)
    int           fStage;               // for this particle
    TStepPointMC* fStepVD9;
    float         fPVD9;                // P(VD9), defined for stopped particles
    float         fY0;
    double        fWeight;
    double        fTau;                 // lifetime , in ns
    bool          fMuonDecay;
    bool          fEleFromMuonDecay;
    TSimParticle* fEle;
  };
}

#endif
