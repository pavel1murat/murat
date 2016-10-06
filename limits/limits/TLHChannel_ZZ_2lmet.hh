#ifndef zzx_limits_TLHChannel_ZZ_llvv
#define zzx_limits_TLHChannel_ZZ_llvv

#include "TLHChannel.hh"

//-----------------------------------------------------------------------------
class TLHChannel_ZZ_2lmet: public TLHChannel {

public:

  TLHChannel_ZZ_2lmet();
  ~TLHChannel_ZZ_2lmet();
					// overloads TLHChannel

  virtual int PseudoExperiment(double& Likelihood);

  ClassDef(TLHChannel_ZZ_2lmet,0) 
};

#endif
