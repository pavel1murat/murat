///////////////////////////////////////////////////////////////////////////////
// base class for plot notes
////////////////////////////////////////////////////////////////////////////////
#include "plot/TPlotNote.hh"
#include "TPaveLabel.h"
#include "TArrow.h"
#include "TObjArray.h"
#include "TObjString.h"

ClassImp(TPlotNote)

//_____________________________________________________________________________
TPlotNote::TPlotNote(int PlotMode, int BlessingMode) {

  //  const char* hist_dir;
  
  fWorkDir      = gSystem->Getenv("WORK_DIR");
  fPlotMode     = PlotMode;
  fBlessingMode = BlessingMode;

  for (int i=0; i<100; i++) {
    fDebugBit[i] = 0;
  }
}


//_____________________________________________________________________________
TPlotNote::~TPlotNote() {
}


//-----------------------------------------------------------------------------
// determine name of the .eps and .gif files in the ./figures directory
// today's default set is 'frr_03', all the plots should be coming from it
//-----------------------------------------------------------------------------
void TPlotNote::get_filename(int Figure, char* Filename) {

  if      (Figure == -1) sprintf(Filename,"fig_%i.%s",Figure,"undefined");
  else {
//-----------------------------------------------------------------------------
// undefined, just figure
//-----------------------------------------------------------------------------
    sprintf(Filename,"fig_%i",Figure);
  }
}


//-----------------------------------------------------------------------------
// 'Font' - ROOT font number
//-----------------------------------------------------------------------------
int TPlotNote::DrawPaveLabelNDC(TPaveLabel*& Label, 
				   const char*  Text , 
				   double       XMin , 
				   double       YMin , 
				   double       XMax , 
				   double       YMax ,
				   int          Font ) {
  Label = new TPaveLabel();

  Label->SetLabel(Text);

  Label->SetTextFont(Font);
  Label->SetTextSize(0.8);
  Label->SetBorderSize(0);
  Label->SetFillStyle(0);

  Label->SetX1NDC(XMin); 
  Label->SetX2NDC(XMax); 
  Label->SetY1NDC(YMin);
  Label->SetY2NDC(YMax);

  Label->Draw();

  return 0;
}

//_____________________________________________________________________________
void TPlotNote::plot(Int_t Figure, const char* CanvasName) {
}


//_____________________________________________________________________________
const char* TPlotNote::GetFiguresDir() {
  return gSystem->Getenv("WORK_DIR");
}


//_____________________________________________________________________________
void TPlotNote::print(int         Figure,
			 const char* Filename, 
			 const char* Dir) {
  // print canvas in .gif and .eps

  char dir[200], filename[200], canvas_name[200];

  sprintf(canvas_name,"fig_%i",Figure);

  if (Dir == 0) {
    if (fPlotMode == kNoteMode ) {
      if (fBlessingMode == 0) sprintf(dir,"%s/note"       ,GetFiguresDir());
      else                    sprintf(dir,"%s/public_note",GetFiguresDir());

    }

    if (fPlotMode == kTalkMode ) sprintf(dir,"%s/talk" ,GetFiguresDir());
    if (fPlotMode == kPaperMode) sprintf(dir,"%s/paper",GetFiguresDir());
  }
  else { 
    strcpy(dir,Dir);
  }

  fCanvas = (TCanvas*) gROOT->GetListOfCanvases()->FindObject(canvas_name);

  if (fCanvas) {
    if (Filename != 0) {
      fCanvas->Print(Form("%s/%s.eps",dir,Filename));
      fCanvas->Print(Form("%s/%s.gif",dir,Filename));
    }
    else {
      this->get_filename(Figure,filename);
      fCanvas->Print(Form("%s/%s.eps",dir,filename));
      fCanvas->Print(Form("%s/%s.gif",dir,filename));
    }
  }
  else {
//-----------------------------------------------------------------------------
// an attempt to print a not-yet plotted figure
//-----------------------------------------------------------------------------
    printf("ERROR: figure %i has not yet been plotted\n",Figure);
  }
}

