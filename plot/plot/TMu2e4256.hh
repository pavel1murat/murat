///////////////////////////////////////////////////////////////////////////////
// plots for mu2e-4256
// -------------------------
// directory with the histogram files: by default: $WORK_DIR/results
// can be redefined with "mu2e.HistDir" in .rootrc
// 
// figures: see description in TMu2e4256.cc
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

class TMu2e4256: public TPlotNote {

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
  char f_e00s1212[200];
  char f_e00s1412[200];
  char f_e00s1512[200];
  char f_e00s1612[200];

  char f_m00s1212[200];
  char f_m00s1412[200];
  char f_m00s1512[200];
  char f_m00s1612[200];
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
public:

  TMu2e4256 (int PlotMode=kNoteMode, int BlessingMode=0);
  ~TMu2e4256();

  void  remake_plots();
//-----------------------------------------------------------------------------
// overloaded methods of TPlotCdfNote
//-----------------------------------------------------------------------------
  virtual void plot        (int Figure, const char* CanvasName = 0);
  virtual void get_filename(int Figure, char* Filename);

  virtual const char*  GetFiguresDir();
//-----------------------------------------------------------------------------
// additional functions
// MuonRecoEff - reconstruction efficiency at zero background occupancy
//-----------------------------------------------------------------------------
  void create_llhr_rejection_graph(const char* FnEle      , 
				   const char* FnMuo      , 
				   double      MuonRecoEff, 
				   TGraph*&    Graph      );

  ClassDef(TMu2e4256,0)
  
};

