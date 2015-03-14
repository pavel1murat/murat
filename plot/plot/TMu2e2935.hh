///////////////////////////////////////////////////////////////////////////////
// plots for mu2e-2935
// -------------------------
// directory with the histogram files: by default: $WORK_DIR/results
// can be redefined with "zzx.HistDir" in .rootrc
// 
// figures:
// ---------
// Fig.   1:  H(300) --> ZZ mass fit
//
////////////////////////////////////////////////////////////////////////////////
#include <stdlib.h>
#include "TArrow.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TEnv.h"
#include "TH2.h"
#include "TInterpreter.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TPad.h"
#include "TEnv.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TPaveLabel.h"
#include "TSystem.h"
#include "TString.h"
#include "TText.h"
#include "Stntuple/val/val/stntuple_val_functions.hh"

#include "murat/plot/plot/TPlotNote.hh"

class TMu2e2935: public TPlotNote {
public:
//-----------------------------------------------------------------------------
// data members
// 1. files:
//-----------------------------------------------------------------------------
  TString   me000001_tcalm003;   // conversions only, Mau7b
  TString   me000201_tcalm002;   // fully mixed events, vane-based, Mau7b
  TString   me000301_tcalm003;   // conversions only, Mau8-uniform-in-DS-downstream
  TString   me000401_tcalm003;   // conversions only, Mau8-full field

  int       fIPlane1;
  double    fStep;
  int       fNPoints;
  double    fZMax;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
public:

  TMu2e2935 (int PlotMode=kNoteMode, int BlessingMode=0);
  ~TMu2e2935();

  void  remake_plots();
//-----------------------------------------------------------------------------
// overloaded methods of TPlotCdfNote
//-----------------------------------------------------------------------------
  virtual void plot        (int Figure, const char* CanvasName = 0);
  virtual void get_filename(int Figure, char* Filename);

  virtual const char*  GetFiguresDir();

  void calculate_eff_vs_z2(const char* Filename, 
			   int         FirstBin, int     NPoints, 
			   float*     X       , float* Ex  , 
			   float*     Eff     , float* err ,
			   float*     eff2    , float* Err2);

  ClassDef(TMu2e2935,0)
  
};

