#ifndef __murat_ana_TrackSeedHist_t_hh__
#define __murat_ana_TrackSeedHist_t_hh__

#include "murat/ana/HistBase_t.h"

namespace murat {

  struct TrackSeedHist_t: public HistBase_t {
    TH1F*    fNHits;	 
    TH1F*    fClusterTime;
    TH1F*    fClusterEnergy;
    TH1F*    fRadius;
    TH1F*    fMom;
    TH1F*    fPt;
    TH1F*    fTanDip;   
    TH1F*    fChi2;
    TH1F*    fFitCons;
    TH1F*    fD0;
  };
}
#endif
