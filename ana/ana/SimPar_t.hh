#ifndef __murat_ana_SimPar_t__
#define __murat_ana_SimPar_t__

  struct SimPar_t {
    TSimParticle*  fParticle;		// pointer to the signal particle
    TVdetHitData*  fTFront;		// VD hit at the tracker front
  };

#endif
