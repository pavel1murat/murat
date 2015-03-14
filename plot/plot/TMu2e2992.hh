///////////////////////////////////////////////////////////////////////////////
// plots for mu2e-2992
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

class TMu2e2992: public TPlotNote {
public:
//-----------------------------------------------------------------------------
// data members
// 1. files:
//-----------------------------------------------------------------------------
  TString e0000001_tcalm002;
  TString cb000101_tcalm002;  // cosmics MC - all tracks |D0| < 20cm, |Z0| < 1m
  TString m0000001_tcalm002;
  TString me000001_tcalm002;
  TString cb000401_tcalm004;
  TString cb000301_cosmicAnalyzer;//resimulation of the dataset cb000201

  TString egun0001_tcalm002;
  TString egun0101_tcalm002;
  TString egun0201_tcalm002;
  TString egun0401_tcalm002;

  TString mgun0001_tcalm002;
  TString mgun0101_tcalm002;
  TString mgun0201_tcalm002;
  TString mgun0401_tcalm002;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
public:

  TMu2e2992 (int PlotMode=kNoteMode, int BlessingMode=0);
  ~TMu2e2992();

  void  remake_plots();
//-----------------------------------------------------------------------------
// overloaded methods of TPlotCdfNote
//-----------------------------------------------------------------------------
  virtual void plot        (int Figure, const char* CanvasName = 0);
  virtual void get_filename(int Figure, char* Filename);

  virtual const char*  GetFiguresDir();

  ClassDef(TMu2e2992,0)
  
};

