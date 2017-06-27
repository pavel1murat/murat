//
#ifndef __murat_scripts_plot_utilities__
#define __murat_scripts_plot_utilities__

#include "TPaveStats.h"
#include "TLatex.h"

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

#endif
