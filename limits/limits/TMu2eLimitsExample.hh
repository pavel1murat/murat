#ifndef zzx_limits_TMu2eLimitsExample
#define zzx_limits_TMu2eLimitsExample

#include "TNamed.h"
#include "TH1.h"
#include "murat/mclimit/mclimit_csm.h"

class TMu2eLimitsExample: public TNamed {

public:

  mclimit_csm*    fMcLimit;

  csm_model*      fNullHyp;
  csm_model*      fNullHypPe;

  csm_model*      fTestHyp;
  csm_model*      fTestHypPe;

  int             fNPseudoExp;              // number of pseudo-experiments

  TH1*            fDataHist;
  TH1*            fBgrHist;
  TH1*            fSigHist;

  double          fXMax;
  double          fXMin;
  int             fNBins;

  double          fBgrLevel;

  TMu2eLimitsExample();
  ~TMu2eLimitsExample();

  int Init();

  int FindLimit(int Npe = 0);

  ClassDef(TMu2eLimitsExample,0) 

};

#endif
