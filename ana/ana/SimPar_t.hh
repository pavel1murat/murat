#ifndef __murat_ana_SimPar_t__
#define __murat_ana_SimPar_t__

  struct SimPar_t {
    TGenParticle*  fGenp;		// 
    TSimParticle*  fParticle;		// pointer to the signal particle
    TVdetHitData*  fTFront;		// VD hit at the tracker front
    TVdetHitData*  fTMid;		// VD hit at the tracker center
  };

#endif
