//
#ifndef __murat_scripts_plot_utilities__
#define __murat_scripts_plot_utilities__

//-----------------------------------------------------------------------------
// redraw the stat box, witht eh color of the histogram itself
//-----------------------------------------------------------------------------
void plot_stat_box(TH1* Hist, double X1, double Y1, double X2, double Y2) {
  TPaveStats* s = (TPaveStats*) Hist->GetListOfFunctions()->FindObject("stats");
  s->SetLineColor(Hist->GetLineColor());
  s->SetTextColor(Hist->GetLineColor());
  s->SetX1NDC(X1); s->SetY1NDC(Y1); s->SetX2NDC(X2); s->SetY2NDC(Y2);
  s->Draw();
}

//-----------------------------------------------------------------------------
void draw_label_ndc(const char* Text, double X1, double Y1, double FontSize, double Font = 52) {

  TLatex* label = new TLatex(X1,Y1,Text);
  label->SetNDC();
  label->SetTextSize(FontSize);
  label->SetTextFont(Font);
  label->Draw();
}

#endif