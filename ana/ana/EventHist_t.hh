#ifndef __murat_ana_EventHist_t_hh__
#define __murat_ana_EventHist_t_hh__

#include "murat/ana/HistBase_t.h"

namespace murat {
  struct EventHist_t : public HistBase_t {
    TH1F*    fEventWeight[2];           // MC truth
    TH1F*    fEventE;                   // MC truth relevant event energy
    TH1F*    fInstLumi[3];              // MC truth lumi info: 0: nominal 1: undo batch weight 2: opposite batch weight
    TH1F*    fBatchWeight[2];           // MC truth
    TH1F*    fLumWeight;                // luminosity-related event weight (PM)
    TH1F*    fRv;			// MC truth information
    TH1F*    fZv;

    TH1F*    fPdgCode;
    TH1F*    fMomTargetEnd;
    TH1F*    fMomTrackerFront;
    TH1F*    fNshCE;

    TH1F*    fMcMom;
    TH1D*    fDioMom;
    TH1F*    fMcCosTh;
    TH1F*    fNHelices;
    TH1F*    fNTracks[2];
    TH1F*    fNGoodTracks[2];
    TH1F*    fNShTot [2];
    TH1F*    fNGoodSH;
    TH1F*    fNChTot [2];
    TH1F*    fDtClT;
    TH1F*    fDtClS;
    TH1F*    fSHTime;
    TH1F*    fNHyp;
    TH1F*    fBestHyp[2];		// [0]: by chi2, [1]: by fit consistency
    TH1F*    fNGenp;                    // N(particles in GENP block)

    TH1F*    fNClusters;
    TH1F*    fEClMax;			// energy of the first (highest) reconstructed cluster
    TH1F*    fTClMax;			// time   of the first (highest) reconstructed cluster
    TH1F*    fDp;                       // P(TrkPatRec)-P(CalPatRec)
    TH1F*    fWeight;			// weight, need with statistics
    TH1F*    fGMom;			// photon momentum
    TH1F*    fGMomRMC;                  // photon momentum, RMC weighted

    TH1F*    fNCrvClusters;
    TH1F*    fNCrvCoincidences[2];
    TH1F*    fNCrvPulses[2];
  };
}
#endif
