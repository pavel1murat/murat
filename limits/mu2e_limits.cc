//-----------------------------------------------------------------------------
// there are 10 mass points
//-----------------------------------------------------------------------------
#include "limits/TMu2eLimits.hh"
#include "limits/mu2e_limits.hh"

TMu2eLimits* gLMu2e;

namespace {
};


//-----------------------------------------------------------------------------
// example: mu2e_limits(2,0.9, 181,325.,100)
// 
// CL : confidence level , i.e. 0.90 or 0.95
//-----------------------------------------------------------------------------
int mu2e_limits(int         Mode  ,
                double      CL    ,
		double      XMin  , 
		double      XMax  , 
		int         NExp  , 
		int         PrintPxFlag) {

  TMu2eLimits::channel_data_t ch_data[] = {
//------------------------------------------------------------------------------------------
// define channels, channel with the name=0 - the last one, use "trk_1" histograms
// 
//  ana init    bit      name           Module    HistSet HistBin  HistName  xmin xmax Rebin
//------------------------------------------------------------------------------------------
    // { 0,  0,   1<<0, "mu2e"      , "Mu2eLimits" ,  "trk" ,   1,   "p"   ,   101.,   106., 1},
    { 0,  0,      0,      0      , 0            ,  0     ,   0,        0,    0.,    -1.,  0}
  };

  if (XMin < XMax) {
    for (int i=0; ch_data[i].name != 0; i++) {
      ch_data[i].xmin = XMin;
      ch_data[i].xmax = XMax;
    }
  }

  gLMu2e = new TMu2eLimits(Mode,ch_data,"sig0");

  gLMu2e->BayesCL(CL,NExp,PrintPxFlag);

  return 0;
}


