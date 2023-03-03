#ifndef __murat_ana_HelixHist_t_hh__
#define __murat_ana_HelixHist_t_hh__

#include "murat/ana/HistBase_t.h"

#include "TH1.h"
#include "TH2.h"

namespace murat {
  struct HelixHist_t : public HistBase_t {
    TH1F*    fBestAlg;
    TH1F*    fCosTh; 
    TH1F*    fChi2XY;
    TH1F*    fChi2ZPhi;
    TH1F*    fD0;
    TH1F*    fHelicity;			// number of ss hits
    TH1F*    fLambda;
    TH1F*    fNCh;			// number of combohits
    TH1F*    fNSh;			// number of ss hits
    TH1F*    fP;			// total momentum, 2 hists with different binning
    TH1F*    fRadius;
  };
}
#endif
