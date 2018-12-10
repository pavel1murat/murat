#ifndef __murat_ana_SimPar_t__
#define __murat_ana_SimPar_t__

  struct SimPar_t {
    TGenParticle*  fGenp;		// 
    TSimParticle*  fParticle;		// pointer to the signal particle
    TVDetHitData*  fTFront;		// VD hit at the tracker front
    TVDetHitData*  fTMid;		// VD hit at the tracker center
    TVDetHitData*  fTBack;		// VD hit at the tracker back (+1610 mm)
  };

#endif
