/// \file
/// \ingroup tutorial_eve
/// Helper classes for the alice_esd_split.C demo.
///
/// \macro_code
///
/// \author Bertrand Bellenot

#ifndef __murat_gui_eve_SplitGlView_hh__
#define __murat_gui_eve_SplitGlView_hh__

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

#include "murat/gui/eve_HtmlObjTable.hh"
#include "murat/gui/eve_HtmlSummary.hh"
#include "murat/gui/eve_TGShapedToolTip.hh"

////////////////////////////////////////////////////////////////////////////////
class SplitGLView : public TGMainFrame {

public:
   enum EMyCommands {
      kFileOpen, kFileExit, kFileLoadConfig, kFileSaveConfig,
      kHelpAbout, kGLPerspYOZ, kGLPerspXOZ, kGLPerspXOY, kGLXOY,
      kGLXOZ, kGLZOY, kGLOrthoRotate, kGLOrthoDolly, kSceneUpdate,
      kSceneUpdateAll, kSummaryUpdate
   };

private:
   TEvePad               *fPad;           // pad used as geometry container
   TGSplitFrame          *fSplitFrame;    // main (first) split frame
   TGLEmbeddedViewer     *fViewer0;       // main GL viewer
   TGLEmbeddedViewer     *fViewer1;       // first GL viewer
   TGLEmbeddedViewer     *fViewer2;       // second GL viewer
   TGLEmbeddedViewer     *fActViewer;     // actual (active) GL viewer
   static HtmlSummary    *fgHtmlSummary;  // summary HTML table
   static TGHtml         *fgHtml;
   TGMenuBar             *fMenuBar;       // main menu bar
   TGPopupMenu           *fMenuFile;      // 'File' popup menu
   TGPopupMenu           *fMenuHelp;      // 'Help' popup menu
   TGPopupMenu           *fMenuCamera;    // 'Camera' popup menu
   TGPopupMenu           *fMenuScene;     // 'Scene' popup menu
   TGStatusBar           *fStatusBar;     // status bar
   TGShapedToolTip       *fShapedToolTip; // shaped tooltip
   Bool_t                 fIsEmbedded;

   TEveViewer            *fViewer[3];
   TEveProjectionManager *fRPhiMgr;
   TEveProjectionManager *fRhoZMgr;

public:
   SplitGLView(const TGWindow *p=0, UInt_t w=800, UInt_t h=600, Bool_t embed=kFALSE);
   virtual ~SplitGLView();

   void           ItemClicked(TGListTreeItem *item, Int_t btn, Int_t x, Int_t y);
   void           HandleMenu(Int_t id);
   void           OnClicked(TObject *obj);
   void           OnMouseIdle(TGLPhysicalShape *shape, UInt_t posx, UInt_t posy);
   void           OnMouseOver(TGLPhysicalShape *shape);
   void           OnViewerActivated();
   void           OpenFile(const char *fname);
   void           SwapToMainView(TGLViewerBase *viewer);
   void           ToggleOrthoRotate();
   void           ToggleOrthoDolly();
   void           UnDock(TGLViewerBase *viewer);
   void           LoadConfig(const char *fname);
   void           SaveConfig(const char *fname);
   static void    UpdateSummary();

   TEveProjectionManager *GetRPhiMgr() const { return fRPhiMgr; }
   TEveProjectionManager *GetRhoZMgr() const { return fRhoZMgr; }

   ClassDef(SplitGLView, 0)
};

#endif
