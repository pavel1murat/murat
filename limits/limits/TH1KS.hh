#ifndef zzx_limits_TH1KS
#define zzx_limits_TH1KS

#include "TH1.h"

class TH1KS : public TH1F {
public:
 
  TH1KS();
  TH1KS(const char* Name, const char* Title, int NBins, float XLow, float XHigh);
  ~TH1KS();

  double KolmogorovTest_Local(const TH1* Hist, int NPExp=1000, Option_t* Opt="") const;

  ClassDef(TH1KS,0)
};


#endif
