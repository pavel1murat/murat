#ifndef zzx_limits_TLHChannel_ZZ_2l2j
#define zzx_limits_TLHChannel_ZZ_2l2j

#include "TLHChannel.hh"


//-----------------------------------------------------------------------------
class TLHChannel_ZZ_2l2j: public TLHChannel {

public:

  TLHChannel_ZZ_2l2j();
  ~TLHChannel_ZZ_2l2j();
					// overloads TLHChannel

  virtual int PseudoExperiment(double& Likelihood);


  ClassDef(TLHChannel_ZZ_2l2j,0) 

};

#endif
