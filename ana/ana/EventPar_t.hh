#ifndef __murat_ana_EventPar_t_hh__
#define __murat_ana_EventPar_t_hh__


namespace murat {

  struct EventPar_t {
    int           fNGenp;
    TGenParticle* fParticle;                // generator "signal" particle
    float         fPartE;

    int           fNTracks    [2];          // might need two for different reasons: DAR<->PAR, ELE<->MUO
    int           fNGoodTracks[2];
    int           fNStrawHits;

    int           fNClusters;               // calorimeter clusters

    int           fNCrvClusters;
    int           fNCrvCoincidences;
    int           fNCrvPulses;

    int           fCandidate_BOX;           // if 1, passes analysis cuts = event candidate
    int           fCandidate_MVA;           // if 1, passes analysis cuts = event candidate

    float         fInstLum;
    double        fOneBatchWeight;
    double        fTwoBatchWeight;
  };
}
#endif
