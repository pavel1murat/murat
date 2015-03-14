#ifndef zzx_limits_zzx_limits
#define zzx_limits_zzx_limits

class TMu2eLimits;

int mu2e_limits(int    Mode,
		double CL          = 0.95,
		double XMin        = 0, 
		double XMax        = -1, 
		int    NExp        = 100, 
		int    PrintPxFlag = 0); 

extern TMu2eLimits* gLMu2e;

#endif
