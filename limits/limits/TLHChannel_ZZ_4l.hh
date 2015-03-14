#ifndef zzx_limits_TLHChannel_ZZ_4l
#define zzx_limits_TLHChannel_ZZ_4l

#include "TH1.h"
#include "TH2.h"

#include "TLHChannel.hh"


//-----------------------------------------------------------------------------
class TLHChannel_ZZ_4l: public TLHChannel {
public:

  double fData[10][2];			// ! M(ZZ), PT(ZZ)

  int    fHighMassOnly;                 // consider only high-mass events

public:

  TLHChannel_ZZ_4l(int HighMassOnly);

  ~TLHChannel_ZZ_4l();
					// overloads TLHChannel

  virtual int Init();
  virtual int PseudoExperiment(double& Likelihood);

  ClassDef(TLHChannel_ZZ_4l,0) 

};

#endif
