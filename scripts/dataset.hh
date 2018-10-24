//-----------------------------------------------------------------------------
// description of a histogram file
//-----------------------------------------------------------------------------
#ifndef __murat_scripts_dataset_hh__
#define __murat_scripts_dataset_hh__

struct dataset_t {

  TString  fName  ;			// dataset name
  TString  fFn    ;			// full name of the histogram file
  TString  fLabel ;			// label to appear on a plot
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

  dataset_t() {
    fName        = "";
    fFn          = "";
    fLabel       = "";
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
  }
  
};

#endif
