///////////////////////////////////////////////////////////////////////////////
// plots for pet-0001
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

class TPet0001: public TPlotNote {
public:
//-----------------------------------------------------------------------------
// data members
// 1. files:
//-----------------------------------------------------------------------------
  TString   fPet002;         // vacuum phantom, point source
  TString   fPet003;         // vacuum phantom, distributed source
  TString   fPet004;         // water  phantom, point source
  TString   fPet005;         // water  phantom, distributed source
  TString   fPet006;         // 32x05  water phantom, single event, distributed source
  TString   fPet007;         // 32x05  water phantom, dose = 1 mCi, distributed source
  TString   fPet008;         // 32x05  water phantom, dose = 1 mCi, distributed source
  TString   fPet009;         // 32x05  water phantom, dose = 1 mCi, distributed source
  TString   fPet010;         // 32x05  water phantom, dose = 1 mCi, distributed source

//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
public:

  TPet0001 (int PlotMode=kNoteMode, int BlessingMode=0);
  ~TPet0001();

  void  remake_plots();
//-----------------------------------------------------------------------------
// overloaded methods of TPlotCdfNote
//-----------------------------------------------------------------------------
  virtual void plot        (int Figure, const char* CanvasName = 0);
  virtual void get_filename(int Figure, char* Filename);

  virtual const char*  GetFiguresDir();

  ClassDef(TPet0001,0)
  
};

