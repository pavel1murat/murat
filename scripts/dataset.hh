//-----------------------------------------------------------------------------
// description of a histogram file
//-----------------------------------------------------------------------------
#ifndef __murat_scripts_dataset_hh__
#define __murat_scripts_dataset_hh__

struct dataset_t {

  TString  fName  ;			// dataset name
  TString  fFn    ;			// full name of the histogram file
  TString  fLabel ;			// label to appear on a plot
  int      fLineColor;
  int      fMarkerStyle;
  float    fMarkerSize;
  int      fMarkerColor;
  float    fXMin;
  float    fXMax;
  float    fYMin;
  float    fYMax;
  long int fNPOT;			// number of generated events (for MC dataset)

  dataset_t() {
    fName        = "";
    fFn          = "";
    fLabel       = "";
    fLineColor   = -1;
    fMarkerSize  = -1;
    fMarkerStyle = -1;
    fMarkerColor = -1;
    fXMin        =  0;
    fXMax        = -1;
    fYMin        =  0;
    fYMax        = -1;
    fNPOT        = -1;
  }
  
};

#endif
