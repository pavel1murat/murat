///////////////////////////////////////////////////////////////////////////////
// plots for mu2e-3280
// -------------------------
// directory with the histogram files: by default: $WORK_DIR/results
// can be redefined with "zzx.HistDir" in .rootrc
// 
// figures: see description in TMu2e3280.cc
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

class TMu2e3280: public TPlotNote {

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
  TString bgr_tcalm002_301[4];
  TString bgr_tcalm002_dev[4];
  TString bgr_tcalm002_nab[4];

  TString bgr_tcalm005_301[4];
  TString bgr_tcalm005_dev[4];
  TString bgr_tcalm005_nab[4];

  TString fBgrName[4];

  double  fQExp   [4];			// expected N/microbunch

  double  fQGen301[4];			// N(generated), v3_0_1
  double  fQGenDev[4];			// N(generated), new PA
  double  fQGenNab[4];			// N(generated), new PA, no neutron absorber

  TH1F*   fHistRC[2];                   // radial distributions of the crystals
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
public:

  TMu2e3280 (int PlotMode=kNoteMode, int BlessingMode=0);
  ~TMu2e3280();

  void  remake_plots();

  double q_background(TString* HistSet, int IBgr);

  TH1F*  occupancy(TString HistSet[], const char* HistName, int IBgr, double* NGen, 
		   const char* HistRE, const char* HistR);
  
  void  plot_energy_vs_r(TPad* Pad, TString HistSet1[], TString HistSet2[], 
			 int IBgr, const char* HistRE, const char* HistR, double EMax = -1.);

  void  plot_nhits_vs_r (TPad* Pad, TString HistSet1[], TString HistSet2[], 
			 int IBgr, const char* HistRN, const char* HistR, double NMax = -1.);

//-----------------------------------------------------------------------------
// overloaded methods of TPlotCdfNote
//-----------------------------------------------------------------------------
  virtual void plot        (int Figure, const char* CanvasName = 0);
  virtual void get_filename(int Figure, char* Filename);

  virtual const char*  GetFiguresDir();

  ClassDef(TMu2e3280,0)
  
};

