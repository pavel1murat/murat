/// \file
/// \ingroup tutorial_eve
/// Helper classes for the alice_esd_split.C demo.
///
/// \macro_code
///
/// \author Bertrand Bellenot
#ifndef __murat_gui_eve_HtmlSummary_hh__
#define __murat_gui_eve_HtmlSummary_hh__

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

#include "murat/gui/eve_HtmlObjTable.hh"

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
class HtmlSummary {
public:                           // make them public for shorter code
   Int_t           fNTables;
   TOrdCollection *fObjTables;    // ->array of object tables
   TString         fHtml;         // output HTML string
   TString         fTitle;        // page title
   TString         fHeader;       // HTML header
   TString         fFooter;       // HTML footer

   void     MakeHeader();
   void     MakeFooter();

public:
   HtmlSummary(const char *title);
   virtual ~HtmlSummary();

   HtmlObjTable  *AddTable(const char *name, Int_t nfields, Int_t nvals,
                           Bool_t exp=kTRUE, Option_t *opt="");
   HtmlObjTable  *GetTable(Int_t at) const { return (HtmlObjTable *)fObjTables->At(at); }
   void           Build();
   void           Clear(Option_t *option="");
   void           Reset(Option_t *option="");
   TString        Html() const { return fHtml; }

   ClassDef(HtmlSummary, 0);
};

#endif
