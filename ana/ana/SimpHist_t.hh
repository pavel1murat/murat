#ifndef __murat_ana_SimpHist_t_hh__
#define __murat_ana_SimpHist_t_hh__

#include "murat/ana/HistBase_t.h"

namespace murat {

 struct SimpHist_t  : public HistBase_t {
   TH1F*    fPdgCode[2];
   TH1F*    fMomTargetEnd;
   TH1F*    fMomTrackerFront;
   TH1F*    fNStrawHits;
   TH1F*    fTime;
   TH1F*    fEndVolumeID;
 };

}
#endif
