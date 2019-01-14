//-----------------------------------------------------------------------------
// a 'dataset' here is a histogram file plus data fields for plotting attributes
//-----------------------------------------------------------------------------
#ifndef __murat_scripts_dataset_hh__
#define __murat_scripts_dataset_hh__

struct dataset_t {

  TString  fName   ;			// dataset name
  TString  fJobName;			// job     name
  TString  fFn     ;			// full name of the histogram file
  TString  fLabel  ;			// label to appear on a plot
  TString  fHistName;                   // redefines the histogram name
  int      fLineColor;                  // these are utility fields to be used as needed
  int      fLineWidth;
  int      fFillColor;                  // these are utility fields to be used as needed
  int      fFillStyle;                  // these are utility fields to be used as needed
  int      fMarkerStyle;
  float    fMarkerSize;
  int      fMarkerColor;
  float    fXMin;
  float    fXMax;
  float    fYMin;
  float    fYMax;
  int      fRebin;
  int      fYLogScale;
  long int fNPOT;			// number of generated events (for MC dataset)
  long int fScale;
  TString  fPlotName;
  TString  fPlotLabel;
  TString  fXAxisTitle;
  float    fLegendXMin;
  float    fLegendYMin;
  float    fLegendXMax;
  float    fLegendYMax;
  TCanvas* fCanvas;
  TH1*     fHist;
  TString  fOutputFn;

  dataset_t() {
    fName        = "";
    fJobName     = "";
    fFn          = "";
    fLabel       = "";
    fHistName    = "";
    fLineColor   = -1;
    fLineWidth   =  1;
    fMarkerSize  = -1;
    fMarkerStyle = -1;
    fMarkerColor = -1;
    fFillColor   = -1;
    fFillStyle   = -1;
    fXMin        =  0;
    fXMax        = -1;
    fYMin        =  0;
    fYMax        = -1;
    fRebin       = -1;
    fYLogScale   =  0;
    fNPOT        = -1;
    fScale       = -1;
    fPlotName    = "";
    fPlotLabel   = "";
    fXAxisTitle  = "";
    fLegendXMin  = -1;
    fLegendYMin  = -1;
    fLegendXMax  = -1;
    fLegendYMax  = -1;
  }
  
};

#endif
