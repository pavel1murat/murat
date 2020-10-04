#ifndef __murat_ana_EventPar_t_hh__
#define __murat_ana_EventPar_t_hh__


namespace murat {

  struct EventPar_t {
    int           fNGenp;
    int           fNTracks[2];
    int           fNGoodTracks[2];
    int           fNStrawHits;
    TGenParticle* fParticle;   // generator particle
    int           fNCrvClusters;
    int           fNCrvCoincidences;
    int           fNCrvPulses;
    int           fCandidate_BOX;           // if 1, passes analysis cuts = event candidate
    int           fCandidate_MVA;           // if 1, passes analysis cuts = event candidate
  };

}
#endif
