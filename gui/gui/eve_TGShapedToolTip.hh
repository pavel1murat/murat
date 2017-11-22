/// \file
/// \ingroup tutorial_eve
/// Helper classes for the alice_esd_split.C demo.
///
/// \macro_code
///
/// \author Bertrand Bellenot
#ifndef __murat_gui_eve_TGShapedToolTip__
#define __murat_gui_eve_TGShapedToolTip__

#include "TApplication.h"
#include "TSystem.h"
#include "TGFrame.h"
#include "TGLayout.h"
#include "TGSplitter.h"
#include "TGLWidget.h"
#include "TEvePad.h"
#include "TGeoManager.h"
#include "TString.h"
#include "TGMenu.h"
#include "TGStatusBar.h"
#include "TGFileDialog.h"
#include "TGMsgBox.h"
#include "TGLPhysicalShape.h"
#include "TGLLogicalShape.h"
#include "HelpText.h"
#include "TClass.h"
#include "Riostream.h"
#include "TEnv.h"
#include "TGListTree.h"
#include "TOrdCollection.h"
#include "TArrayF.h"
#include "TGHtml.h"
#include "TPRegexp.h"

#include "TEveManager.h"
#include "TEveViewer.h"
#include "TEveBrowser.h"
#include "TEveProjectionManager.h"
#include "TEveProjectionAxes.h"
#include "TEveScene.h"
#include "TEveGeoNode.h"
#include "TEveEventManager.h"
#include "TEveTrack.h"
#include "TEveSelection.h"

#include "TRootEmbeddedCanvas.h"
#include "TGSplitFrame.h"
#include "TGLOverlayButton.h"
#include "TGLEmbeddedViewer.h"
#include "TGDockableFrame.h"
#include "TGShapedFrame.h"
#include "TGButton.h"
#include "TGTab.h"
#include "TEnv.h"

#include "TCanvas.h"
#include "TFormula.h"
#include "TF1.h"
#include "TH1F.h"

// const char *filetypes[] = {
//    "ROOT files",    "*.root",
//    "All files",     "*",
//    0,               0
// };

// const char *rcfiletypes[] = {
//    "All files",     "*",
//    0,               0
// };

////////////////////////////////////////////////////////////////////////////////
class TGShapedToolTip : public TGShapedFrame {

private:
   TGShapedToolTip(const TGShapedToolTip&); // Not implemented
   TGShapedToolTip& operator=(const TGShapedToolTip&); // Not implemented

protected:
   Int_t                 fTextX, fTextY, fTextH;
   TString               fTextCol;

   TRootEmbeddedCanvas  *fEc;       // embedded canvas for histogram
   TH1                  *fHist;     // user histogram
   TString               fText;     // info (as tool tip) text

   virtual void          DoRedraw() {}

public:
   TGShapedToolTip(const char *picname, Int_t cx=0, Int_t cy=0, Int_t cw=0,
                   Int_t ch=0, Int_t tx=0, Int_t ty=0, Int_t th=0,
                   const char *col="#ffffff");
   virtual ~TGShapedToolTip();

   virtual void   CloseWindow();
   void           CreateCanvas(Int_t cx, Int_t cy, Int_t cw, Int_t ch);
   void           CreateCanvas(Int_t cw, Int_t ch, TGLayoutHints *hints);
   TH1           *GetHisto() const { return fHist; }
   const char    *GetText() const { return fText.Data(); }
   void           Refresh();
   void           SetHisto(TH1 *hist);
   void           SetText(const char *text);
   void           SetTextColor(const char *col);
   void           SetTextAttributes(Int_t tx, Int_t ty, Int_t th, const char *col=0);
   void           Show(Int_t x, Int_t y, const char *text = 0, TH1 *hist = 0);

   ClassDef(TGShapedToolTip, 0) // Shaped composite frame
};

#endif
