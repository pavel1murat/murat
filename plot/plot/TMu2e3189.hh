///////////////////////////////////////////////////////////////////////////////
// plots for mu2e-3189
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

class TMu2e3189: public TPlotNote {
public:
//-----------------------------------------------------------------------------
// data members
// 1. files:
//-----------------------------------------------------------------------------
  TString dio_00_050_80c  ;
  TString dio_50_105_50b  ;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
public:

  TMu2e3189 (int PlotMode=kNoteMode, int BlessingMode=0);
  ~TMu2e3189();

  void  remake_plots();
//-----------------------------------------------------------------------------
// overloaded methods of TPlotCdfNote
//-----------------------------------------------------------------------------
  virtual void plot        (int Figure, const char* CanvasName = 0);
  virtual void get_filename(int Figure, char* Filename);

  virtual const char*  GetFiguresDir();

  ClassDef(TMu2e3189,0)
  
};

