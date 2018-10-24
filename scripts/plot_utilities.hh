//
#ifndef __murat_scripts_plot_utilities__
#define __murat_scripts_plot_utilities__

#include "TPaveStats.h"
#include "TLatex.h"

#include "murat/scripts/dataset.hh"
//-----------------------------------------------------------------------------
void set_draw_options(TH1* Hist, dataset_t* Ds) {
  if (Ds->fLineColor   > 0) Hist->SetLineColor(Ds->fLineColor);
  if (Ds->fMarkerSize  > 0) Hist->SetMarkerSize (Ds->fMarkerSize);
  if (Ds->fMarkerColor > 0) Hist->SetMarkerColor(Ds->fMarkerColor);
  if (Ds->fMarkerStyle > 0) Hist->SetMarkerStyle(Ds->fMarkerStyle);

  if (Ds->fXMin < Ds->fXMax) Hist->GetXaxis()->SetRangeUser(Ds->fXMin,Ds->fXMax);
  if (Ds->fYMin < Ds->fYMax) Hist->GetYaxis()->SetRangeUser(Ds->fYMin,Ds->fYMax);
}

//-----------------------------------------------------------------------------
// redraw the stat box, with the color of the histogram itself
//-----------------------------------------------------------------------------
void plot_stat_box(TH1* Hist, double X1, double Y1, double X2, double Y2) {
  TPaveStats* s = (TPaveStats*) Hist->GetListOfFunctions()->FindObject("stats");
  if (s != NULL) {
    s->SetLineColor(Hist->GetLineColor());
    s->SetTextColor(Hist->GetLineColor());
    s->SetX1NDC(X1); s->SetY1NDC(Y1); s->SetX2NDC(X2); s->SetY2NDC(Y2);
    s->Draw();
  }
  else {
    printf("ERROR: stat box is not defined for %s\n",Hist->GetName());
  }
}

//-----------------------------------------------------------------------------
// draw label in normalized coordinates (x and y range from 0 to 1
//-----------------------------------------------------------------------------
void draw_label_ndc(const char* Text, double X1, double Y1, double FontSize, int Font = 52) {

  TLatex* label = new TLatex(X1,Y1,Text);
  label->SetNDC(kTRUE);
  label->SetTextSize(FontSize);
  label->SetTextFont(Font);
  label->Draw();
}

//-----------------------------------------------------------------------------
// draw a label in absolute coordinates
//-----------------------------------------------------------------------------
void draw_label_abs(const char* Text, double X1, double Y1, double FontSize, int Font = 52) {

  TLatex* label = new TLatex(X1,Y1,Text);
  label->SetNDC(kFALSE);
  label->SetTextSize(FontSize);
  label->SetTextFont(Font);
  label->Draw();
}

//-----------------------------------------------------------------------------
void make_empty_plot(const char* figure_name, int Print = 0) {

  TCanvas* c = new TCanvas("c_empty","EmptyPlot",1500,800);

  draw_label_ndc("EMPTY PLOT :  - don't forget to make !",0.25,0.4,0.06,52);

  // TArrow* arr = new TArrow(0.9,40.e3,0.9,10.e3,0.015);

  // arr->SetLineWidth(2);
  // arr->Draw();

  if (Print == 1) c->Print(figure_name);
}


#endif
