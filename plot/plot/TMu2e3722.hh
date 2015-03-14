///////////////////////////////////////////////////////////////////////////////
// plots for mu2e-3722
// -------------------------
// directory with the histogram files: by default: $WORK_DIR/results
// can be redefined with "mu2e.HistDir" in .rootrc
// 
// figures: see description in TMu2e3722.cc
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
#include "TGraphAsymmErrors.h"
#include "Stntuple/val/val/stntuple_val_functions.hh"

#include "murat/plot/plot/TPlotNote.hh"

class TMu2e3722: public TPlotNote {

  enum {
    kDio   = 0,
    kEpr   = 1,
    kEnu   = 2,
    kEph   = 3
  };

public:
//-----------------------------------------------------------------------------
// data members
// 1. files:
//-----------------------------------------------------------------------------
  char f_conv01[200];
  char f_cosm01[200];
  char f_data01[200];
  char f_dio01 [200];
  char f_dio02 [200];

  char f_080   [200];
  char f_087   [200];
  char f_090   [200];
  char f_095   [200];
  char f_100   [200];
  char f_105   [200];

  double    fNConvEle; // expected number of CE
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
public:

  TMu2e3722 (int PlotMode=kNoteMode, int BlessingMode=0);
  ~TMu2e3722();

  void  remake_plots();

  void make_graphs(double* Data, TGraphAsymmErrors  **G1, TGraphAsymmErrors **G2);
  void eff_vs_mom (TGraphErrors** G);
//-----------------------------------------------------------------------------
// resolution fit
//-----------------------------------------------------------------------------
  static double f_crystal_ball(double* X, double* P);
  static double f_eff_mom     (double* X, double* P);
//-----------------------------------------------------------------------------
// overloaded methods of TPlotCdfNote
//-----------------------------------------------------------------------------
  virtual void plot        (int Figure, const char* CanvasName = 0);
  virtual void get_filename(int Figure, char* Filename);

  virtual const char*  GetFiguresDir();

  ClassDef(TMu2e3722,0)
  
};

