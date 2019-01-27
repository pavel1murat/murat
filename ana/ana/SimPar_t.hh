#ifndef __murat_ana_SimPar_t__
#define __murat_ana_SimPar_t__

  struct SimPar_t {
    TGenParticle*  fGenp;		// 
    TSimParticle*  fParticle;		// pointer to the signal particle
    TStepPointMC*  fTFront;		// VD hit at the tracker front
    TStepPointMC*  fTMid;		// VD hit at the tracker center
    TStepPointMC*  fTBack;		// VD hit at the tracker back (+1610 mm)
  };

#endif
