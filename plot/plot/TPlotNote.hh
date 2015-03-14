///////////////////////////////////////////////////////////////////////////////
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



class TPlotNote {
public:
//-----------------------------------------------------------------------------
// data members
//-----------------------------------------------------------------------------
  enum { 
    kNoteMode  = 1,
    kTalkMode  = 2,
    kPaperMode = 3
  };
//-----------------------------------------------------------------------------
// data members
//-----------------------------------------------------------------------------
  TString       fWorkDir;
  TString       fFiguresDir;
  int           fPlotMode;        // 1:note 2:talk
  int           fBlessingMode;    // 1:plot for blessing

  TCanvas       *fCanvas;               // !
  TPad          *fP1;                   // !
  int           fDebugBit[100];
//-----------------------------------------------------------------------------
// histograms
//-----------------------------------------------------------------------------
  TH1           *fH1[1000];		// temp hist
  TH2           *fH2[1000];		// temp hist
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
public:

  TPlotNote (int PlotMode=kNoteMode, int BlessingMode=0);
  virtual ~TPlotNote();

  int  GetDebugBit(int I) { return fDebugBit[I]; }

  void SetDebugBit(int I, int Val) { fDebugBit[I] = Val; }

  virtual void plot        (int Figure, const char* CanvasName=0);
  virtual void print       (int Figure, const char* Filename=0, const char* Dir=0);
  virtual void get_filename(int Figure, char* Filename);

  virtual const char* GetFiguresDir();

  int  DrawPaveLabelNDC(TPaveLabel*& Label, 
			const char*  Text , 
			double       XMin , 
			double       YMin , 
			double       XMax , 
			double       YMax ,
			int          Font = 52);

  ClassDef(TPlotNote,0)
  
};

