///////////////////////////////////////////////////////////////////////////////
// plots for mu2e-Physics
// -------------------------
// directory with the histogram files: by default: $WORK_DIR/results
// can be redefined with "mu2e.HistDir" in .rootrc
// 
// figures: see description in TMu2ePhysics.cc
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

class TMu2ePhysics: public TPlotNote {

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
  char f_conv[200];
  char f_cosm[200];
  char f_data[200];
  char f_dio [200];
  char f_rpc [200];

  double    fNConvEle; // expected number of CE
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
public:

  TMu2ePhysics (int PlotMode=kNoteMode, int BlessingMode=0);
  ~TMu2ePhysics();

  void  remake_plots();
//-----------------------------------------------------------------------------
// overloaded methods of TPlotCdfNote
//-----------------------------------------------------------------------------
  virtual void plot        (int Figure, const char* CanvasName = 0);
  virtual void get_filename(int Figure, char* Filename);

  virtual const char*  GetFiguresDir();

  ClassDef(TMu2ePhysics,0)
  
};

