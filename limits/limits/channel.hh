#ifndef __murat_limits_channel_h__
#define __murat_limits_channel_h__

#include "TNamed.h"
#include "TH1.h"
#include "TH2.h"

namespace murat {
class  channel : public TNamed {
public:
  
  float  kLumiSF_1B;	     // scale factor for 1-batch running mode 
  float  kLumiSF_2B;	     // scale factor for 2-batch running mode
  float  fCapture;           // muon capture fraction in Al
  float  fMusr;              // muon stopping rate, hide it here
  TH2F*  fTimeVsMom;	     // main histogram, sum of the two components (could be more, in principle)
  TH1F*  fHist;
  TString fHistDir;
// -----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  channel(const char* ChannelName, float SF_1B, float SF_2B, float ExtraSF = 1.);

  double GetIntegral(float PMin = -1, float PMax=1.e6, float TMin = -1., float TMax = 1.e6);

  virtual double FluctuateBackground(float PMin, float PMax, float TMin, float TMax);

  const char*  GetHistDir() { return fHistDir.Data(); }
  
  TH1D*  CreateMomHist (double TMin = -1, double TMax=1.e6);
  TH1D*  CreateTimeHist(double PMin = -1, double PMax=1.e6);
//-----------------------------------------------------------------------------
// create running integral (T_i, TMax) histogram with given timing cutoff for
// the time distribution within given (PMin,PMax) momentum band
//-----------------------------------------------------------------------------
  TH1D*  CreateTimeIntegralHist(double PMin = -1, double PMax=1.e6, double TMin = -1., double TMax = 1.e6);
//-----------------------------------------------------------------------------
// create running integral (T_i, TMax) histogram with given timing cutoff for
// the time distribution within given (PMin,PMax) momentum band
//-----------------------------------------------------------------------------
  TH1D* CreateMomIntegralHist(double PMin = -1, double PMax=1.e6, double TMin = -1., double TMax = 1.e6);
};
}
#endif
