#ifndef __murat_ana_GenpHist_t_hh__
#define __murat_ana_GenpHist_t_hh__

#include "murat/ana/HistBase_t.h"

namespace murat {

 struct GenpHist_t  : public HistBase_t {
   TH1F*    fPdgCode[2];		// same distribution in different scale
   TH1F*    fGenID;			// 
   TH1F*    fZ0;			// 
   TH1F*    fT0;			// 
   TH1F*    fR0;			// 
   TH1F*    fP;			// 
   TH1F*    fCosTh;			// 
 };

}
#endif
