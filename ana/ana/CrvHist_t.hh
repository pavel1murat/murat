#ifndef __murat_ana_CrvHist_t_hh__
#define __murat_ana_CrvHist_t_hh__

#include "murat/ana/HistBase_t.h"

#include "TH1.h"
#include "TH2.h"

namespace murat {

  struct CrvClusterHist_t {
    TH1F*    fSectorType;
    TH1F*    fNPulses;
    TH1F*    fNPe;
    TH1F*    fStartTime;
    TH1F*    fEndTime;
    TH1F*    fWidth;
    TH2F*    fXVsZ;
    TH2F*    fYVsZ;
  };

  struct CrvPulseHist_t {
    TH1F*    fNPe;
    TH1F*    fNPeHeight;
    TH1F*    fNDigis;
    TH1F*    fBar;
    TH1F*    fSipm;
    TH1F*    fTime;
    TH1F*    fHeight;
    TH1F*    fWidth;
    TH1F*    fChi2;
    TH1F*    fLeTime;
    TH1F*    fDt;
  };

  struct CrvCoincidenceHist_t {
    TH1F*    fSectorType;
    TH1F*    fNPulses;
  };

}
#endif
