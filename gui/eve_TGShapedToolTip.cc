///////////////////////////////////////////////////////////////////////////////
#include "TVirtualX.h"
#include "murat/gui/eve_TGShapedToolTip.hh"

ClassImp(TGShapedToolTip)

//______________________________________________________________________________
TGShapedToolTip::TGShapedToolTip(const char *pname, Int_t cx, Int_t cy, Int_t cw,
                             Int_t ch, Int_t tx, Int_t ty, Int_t th,
                             const char *col) :
   TGShapedFrame(pname, gClient->GetDefaultRoot(), 400, 300, kTempFrame |
                 kHorizontalFrame), fEc(0), fHist(0)
{
   // Shaped window constructor

   fTextX = tx; fTextY = ty; fTextH = th;
   if (col)
      fTextCol = col;
   else
      fTextCol = "0x000000";

   // create the embedded canvas
   if ((cx > 0) && (cy > 0) && (cw > 0) && (ch > 0)) {
      Int_t lhRight  = fWidth-cx-cw;
      Int_t lhBottom = fHeight-cy-ch;
      fEc = new TRootEmbeddedCanvas("ec", this, cw, ch, 0);
      AddFrame(fEc, new TGLayoutHints(kLHintsTop | kLHintsLeft, cx,
                                      lhRight, cy, lhBottom));
   }
   MapSubwindows();
   Resize();
   Resize(fBgnd->GetWidth(), fBgnd->GetHeight());
}

//______________________________________________________________________________
TGShapedToolTip::~TGShapedToolTip()
{
   // Destructor.

   if (fHist)
      delete fHist;
   if (fEc)
      delete fEc;
}

//______________________________________________________________________________
void TGShapedToolTip::CloseWindow()
{
   // Close shaped window.

   DeleteWindow();
}

//______________________________________________________________________________
void TGShapedToolTip::Refresh()
{
   // Redraw the window with current attributes.

   const char *str = fText.Data();
   char *string = strdup(str);
   Int_t nlines = 0, size = fTextH;
   TString fp = gEnv->GetValue("Root.TTFontPath", "");
   TString ar = fp + "/arial.ttf";
   char *s = strtok((char *)string, "\n");
   TImage *img = (TImage*)fImage->Clone("img");
   img->DrawText(fTextX, fTextY+(nlines*size), s, size, fTextCol, ar);
   while ((s = strtok(0, "\n"))) {
      nlines++;
      img->DrawText(fTextX, fTextY+(nlines*size), s, size, fTextCol, ar);
   }
   img->PaintImage(fId, 0, 0, 0, 0, 0, 0, "opaque");
   free(string);
   delete img;
   gVirtualX->Update();
}

//______________________________________________________________________________
void TGShapedToolTip::CreateCanvas(Int_t cx, Int_t cy, Int_t cw, Int_t ch)
{

   // create the embedded canvas
   Int_t lhRight  = fWidth-cx-cw;
   Int_t lhBottom = fHeight-cy-ch;
   fEc = new TRootEmbeddedCanvas("ec", this, cw, ch, 0);
   AddFrame(fEc, new TGLayoutHints(kLHintsTop | kLHintsLeft, cx,
                                   lhRight, cy, lhBottom));
   MapSubwindows();
   Resize();
   Resize(fBgnd->GetWidth(), fBgnd->GetHeight());
   if (IsMapped()) {
      Refresh();
   }
}

//______________________________________________________________________________
void TGShapedToolTip::CreateCanvas(Int_t cw, Int_t ch, TGLayoutHints *hints)
{
   // Create the embedded canvas.

   fEc = new TRootEmbeddedCanvas("ec", this, cw, ch, 0);
   AddFrame(fEc, hints);
   MapSubwindows();
   Resize();
   Resize(fBgnd->GetWidth(), fBgnd->GetHeight());
   if (IsMapped()) {
      Refresh();
   }
}

//______________________________________________________________________________
void TGShapedToolTip::SetHisto(TH1 *hist)
{
   // Set which histogram has to be displayed in the embedded canvas.

   if (hist) {
      if (fHist) {
         delete fHist;
         if (fEc)
            fEc->GetCanvas()->Clear();
      }
      fHist = (TH1 *)hist->Clone();
      if (fEc) {
         fEc->GetCanvas()->SetBorderMode(0);
         fEc->GetCanvas()->SetFillColor(10);
         fEc->GetCanvas()->cd();
         fHist->Draw();
         fEc->GetCanvas()->Update();
      }
   }
}

//______________________________________________________________________________
void TGShapedToolTip::SetText(const char *text)
{
   // Set which text has to be displayed.

   if (text) {
      fText = text;
   }
   if (IsMapped())
      Refresh();
}

//______________________________________________________________________________
void TGShapedToolTip::SetTextColor(const char *col)
{
   // Set text color.

   fTextCol = col;
   if (IsMapped())
      Refresh();
}

//______________________________________________________________________________
void TGShapedToolTip::SetTextAttributes(Int_t tx, Int_t ty, Int_t th,
                                        const char *col)
{
   // Set text attributes (position, size and color).

   fTextX = tx; fTextY = ty; fTextH = th;
   if (col)
      fTextCol = col;
   if (IsMapped())
      Refresh();
}

//______________________________________________________________________________
void TGShapedToolTip::Show(Int_t x, Int_t y, const char *text, TH1 *hist)
{
   // Show (popup) the shaped window at location x,y and possibly
   // set the text and histogram to be displayed.

   Move(x, y);
   MapWindow();

   if (text)
      SetText(text);
   if (hist)
      SetHisto(hist);
   // end of demo code -------------------------------------------
   if (fHist) {
      fEc->GetCanvas()->SetBorderMode(0);
      fEc->GetCanvas()->SetFillColor(10);
      fEc->GetCanvas()->cd();
      fHist->Draw();
      fEc->GetCanvas()->Update();
   }
   Refresh();
}

