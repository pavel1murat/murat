/////////////////////////////////////////////////////////////////////////////

#ifndef __murat_limits_mumem_hh__
#define __murat_limits_mumem_hh__

#include "Stntuple/scripts/stn_catalog.hh"

#include "Stntuple/scripts/plot_data.hh"
#include "Stntuple/scripts/plot_hist_1D.C"
#include "Stntuple/scripts/plot_hist_2D.C"
#include "Stntuple/alg/TFeldmanCousinsA.hh"

#include "murat/limits/analysis.hh"

namespace murat {
  
class mumem : public murat::analysis {
public:
  double  fSFSignal;            // try to rescale to account for tan(dip) cut
  double  fSFCosmics;		//
  double  fPMin;		// signal window for plotting histograms
  double  fPMax;
  double  fTMin;
  double  fTMax;
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
  mumem(int Mode = 1, const char* Name = "mu-e-", const char* Title = "mu- --> e- conversion"); 
//-----------------------------------------------------------------------------
// other functions
//-----------------------------------------------------------------------------
  int  init_channels(int Mode);
  void plot          (int Figure, int Plot = -1);
  int  scan_pmin     (double P0=103.85, double PMax=105.0, double TMin=700, double TMax = 1700);
  int  scan_tmin     (double P0=103.85, double PMax=105.0, double TMin=700, double TMax = 1700);
  int  scan_tmax     (double P0=103.85, double PMax=105.0, double TMin=700, double TMax = 1700);
//-----------------------------------------------------------------------------
// overloaded function of the base class
//-----------------------------------------------------------------------------
  virtual double FluctuateBackground();
};

}
#endif
