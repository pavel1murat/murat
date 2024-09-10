#ifndef __murat_ana_SimpHist_t_hh__
#define __murat_ana_SimpHist_t_hh__

#include "murat/ana/HistBase_t.h"

namespace murat {

 struct SimpHist_t  : public HistBase_t {
   TH1F*    fStage;
   TH1F*    fPdgCode[2];
   TH1F*    fParentMom;
   TH1F*    fParentPDG;
   TH1F*    fGeneratorID;

   TH1F*    fStartVolumeID;
   TH1F*    fEndVolumeID;

   TH1F*    fStartTime;
   TH1F*    fEndTime;
   TH1F*    fStageDt;

   TH1F*    fStartMom[3];
   TH1F*    fEndMom  [3];
   TH1F*    fNStrawHits;

   TH1F*    fMomTargetEnd;
   TH1F*    fMomTrackerFront;
   TH1F*    fCosTh;
   TH1F*    fZStart;
   TH1F*    fZEnd;

   TH2F*    fNshVsCosTh;
   
   TH2F*    fYVsX;
   TH2F*    fYeVsP;
   TH2F*    fXEndVsZEnd;
   TH2F*    fYcVsZEnd;
   TH2F*    fPVD9VsZEnd;
   TH2F*    fYVsX_2480;
   TH2F*    fYVsX_2513;

   TH2F*    fCosThVsMom[2];
   TH2F*    fCosThVsMomPV;
 };

}
#endif
