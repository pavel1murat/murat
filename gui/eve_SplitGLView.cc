/// \file
/// \ingroup tutorial_eve
/// Helper classes for the alice_esd_split.C demo.
///
/// \macro_code
///
/// \author Bertrand Bellenot
//-----------------------------------------------------------------------------
#include "TVirtualX.h"

#include "murat/gui/eve_SplitGLView.hh"
#include "murat/gui/eve_HtmlObjTable.hh"
#include "murat/gui/eve_HtmlSummary.hh"

const char *filetypes[] = {
   "ROOT files",    "*.root",
   "All files",     "*",
   0,               0
};

const char *rcfiletypes[] = {
   "All files",     "*",
   0,               0
};

TEveProjectionManager *gRPhiMgr = 0;
TEveProjectionManager *gRhoZMgr = 0;

ClassImp(SplitGLView)

HtmlSummary *SplitGLView::fgHtmlSummary = NULL;
TGHtml      *SplitGLView::fgHtml        = NULL;

//______________________________________________________________________________
SplitGLView::SplitGLView(const TGWindow *p, UInt_t w, UInt_t h, Bool_t embed) :
   TGMainFrame(p, w, h), fActViewer(0), fShapedToolTip(0), fIsEmbedded(embed)
{
   // Main frame constructor.

   TGSplitFrame *frm;
   TEveScene *s = 0;
   TGHorizontalFrame *hfrm;
   TGDockableFrame *dfrm;
   TGPictureButton *button;

   // create the "file" popup menu
   fMenuFile = new TGPopupMenu(gClient->GetRoot());
   fMenuFile->AddEntry("&Open...", kFileOpen);
   fMenuFile->AddSeparator();
   fMenuFile->AddEntry( "&Update Summary", kSummaryUpdate);
   fMenuFile->AddSeparator();
   fMenuFile->AddEntry("&Load Config...", kFileLoadConfig);
   fMenuFile->AddEntry("&Save Config...", kFileSaveConfig);
   fMenuFile->AddSeparator();
   fMenuFile->AddEntry("E&xit", kFileExit);

   // create the "camera" popup menu
   fMenuCamera = new TGPopupMenu(gClient->GetRoot());
   fMenuCamera->AddEntry("Perspective (Floor XOZ)", kGLPerspXOZ);
   fMenuCamera->AddEntry("Perspective (Floor YOZ)", kGLPerspYOZ);
   fMenuCamera->AddEntry("Perspective (Floor XOY)", kGLPerspXOY);
   fMenuCamera->AddEntry("Orthographic (XOY)", kGLXOY);
   fMenuCamera->AddEntry("Orthographic (XOZ)", kGLXOZ);
   fMenuCamera->AddEntry("Orthographic (ZOY)", kGLZOY);
   fMenuCamera->AddSeparator();
   fMenuCamera->AddEntry("Ortho allow rotate", kGLOrthoRotate);
   fMenuCamera->AddEntry("Ortho allow dolly",  kGLOrthoDolly);

   fMenuScene = new TGPopupMenu(gClient->GetRoot());
   fMenuScene->AddEntry("&Update Current", kSceneUpdate);
   fMenuScene->AddEntry("Update &All", kSceneUpdateAll);

   // create the "help" popup menu
   fMenuHelp = new TGPopupMenu(gClient->GetRoot());
   fMenuHelp->AddEntry("&About", kHelpAbout);

   // create the main menu bar
   fMenuBar = new TGMenuBar(this, 1, 1, kHorizontalFrame);
   fMenuBar->AddPopup("&File", fMenuFile, new TGLayoutHints(kLHintsTop |
                      kLHintsLeft, 0, 4, 0, 0));
   fMenuBar->AddPopup("&Camera", fMenuCamera, new TGLayoutHints(kLHintsTop |
                      kLHintsLeft, 0, 4, 0, 0));
   fMenuBar->AddPopup("&Scene", fMenuScene, new TGLayoutHints(kLHintsTop |
                      kLHintsLeft, 0, 4, 0, 0));
   fMenuBar->AddPopup("&Help", fMenuHelp, new TGLayoutHints(kLHintsTop |
                      kLHintsRight));

   AddFrame(fMenuBar, new TGLayoutHints(kLHintsTop | kLHintsExpandX));

   // connect menu signals to our menu handler slot
   fMenuFile->Connect("Activated(Int_t)", "SplitGLView", this,
                      "HandleMenu(Int_t)");
   fMenuCamera->Connect("Activated(Int_t)", "SplitGLView", this,
                        "HandleMenu(Int_t)");
   fMenuScene->Connect("Activated(Int_t)", "SplitGLView", this,
                       "HandleMenu(Int_t)");
   fMenuHelp->Connect("Activated(Int_t)", "SplitGLView", this,
                      "HandleMenu(Int_t)");

   if (fIsEmbedded && gEve) {
      // use status bar from the browser
      fStatusBar = gEve->GetBrowser()->GetStatusBar();
   }
   else {
      // create the status bar
      Int_t parts[] = {45, 15, 10, 30};
      fStatusBar = new TGStatusBar(this, 50, 10);
      fStatusBar->SetParts(parts, 4);
      AddFrame(fStatusBar, new TGLayoutHints(kLHintsBottom | kLHintsExpandX,
               0, 0, 10, 0));
   }

   // create eve pad (our geometry container)
   fPad = new TEvePad();
   fPad->SetFillColor(kBlack);

   // create the split frames
   fSplitFrame = new TGSplitFrame(this, 800, 600);
   AddFrame(fSplitFrame, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
   // split it once
   fSplitFrame->HSplit(434);
   // then split each part again (this will make four parts)
   fSplitFrame->GetSecond()->VSplit(266);
   fSplitFrame->GetSecond()->GetSecond()->VSplit(266);

   TGLOverlayButton *but1, *but2, *but3, *but4, *but5, *but6;
   // get top (main) split frame
   frm = fSplitFrame->GetFirst();
   frm->SetName("Main_View");

   // create (embed) a GL viewer inside
   fViewer0 = new TGLEmbeddedViewer(frm, fPad);
   but1 = new TGLOverlayButton(fViewer0, "Swap", 10.0, -10.0, 55.0, 16.0);
   but1->Connect("Clicked(TGLViewerBase*)", "SplitGLView", this, "SwapToMainView(TGLViewerBase*)");
   but2 = new TGLOverlayButton(fViewer0, "Undock", 70.0, -10.0, 55.0, 16.0);
   but2->Connect("Clicked(TGLViewerBase*)", "SplitGLView", this, "UnDock(TGLViewerBase*)");
   frm->AddFrame(fViewer0->GetFrame(), new TGLayoutHints(kLHintsExpandX |
                 kLHintsExpandY));
   // set the camera to perspective (XOZ) for this viewer
   fViewer0->SetCurrentCamera(TGLViewer::kCameraPerspXOZ);
   // connect signal we are interested to
   fViewer0->Connect("MouseOver(TGLPhysicalShape*)", "SplitGLView", this,
                      "OnMouseOver(TGLPhysicalShape*)");
   fViewer0->Connect("Activated()", "SplitGLView", this,
                      "OnViewerActivated()");
   fViewer0->Connect("MouseIdle(TGLPhysicalShape*,UInt_t,UInt_t)",
                      "SplitGLView", this,
                      "OnMouseIdle(TGLPhysicalShape*,UInt_t,UInt_t)");
   fViewer0->Connect("Clicked(TObject*)", "SplitGLView", this,
                      "OnClicked(TObject*)");
   fViewer[0] = new TEveViewer("SplitGLViewer[0]");
   fViewer[0]->SetGLViewer(fViewer0, fViewer0->GetFrame());
   fViewer[0]->IncDenyDestroy();
   if (fIsEmbedded && gEve) {
      fViewer[0]->AddScene(gEve->GetGlobalScene());
      fViewer[0]->AddScene(gEve->GetEventScene());
      gEve->GetViewers()->AddElement(fViewer[0]);
      s = gEve->SpawnNewScene("Rho-Z Projection");
      // projections
      fRhoZMgr = new TEveProjectionManager(TEveProjection::kPT_RhoZ);
      s->AddElement(fRhoZMgr);
      gEve->AddToListTree(fRhoZMgr, kTRUE);
      TEveProjectionAxes* a = new TEveProjectionAxes(fRhoZMgr);
      s->AddElement(a);
   }

   // get bottom left split frame
   frm = fSplitFrame->GetSecond()->GetFirst();
   frm->SetName("Bottom_Left");

   // create (embed) a GL viewer inside
   fViewer1 = new TGLEmbeddedViewer(frm, fPad);
   but3 = new TGLOverlayButton(fViewer1, "Swap", 10.0, -10.0, 55.0, 16.0);
   but3->Connect("Clicked(TGLViewerBase*)", "SplitGLView", this, "SwapToMainView(TGLViewerBase*)");
   but4 = new TGLOverlayButton(fViewer1, "Undock", 70.0, -10.0, 55.0, 16.0);
   but4->Connect("Clicked(TGLViewerBase*)", "SplitGLView", this, "UnDock(TGLViewerBase*)");
   frm->AddFrame(fViewer1->GetFrame(),new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

   // set the camera to orthographic (XOY) for this viewer
   fViewer1->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
   // connect signal we are interested to
   fViewer1->Connect("MouseOver(TGLPhysicalShape*)", "SplitGLView", this,"OnMouseOver(TGLPhysicalShape*)");
   fViewer1->Connect("Activated()"                 , "SplitGLView", this,"OnViewerActivated()");
   fViewer1->Connect("MouseIdle(TGLPhysicalShape*,UInt_t,UInt_t)","SplitGLView",this,"OnMouseIdle(TGLPhysicalShape*,UInt_t,UInt_t)");
   fViewer1->Connect("Clicked(TObject*)", "SplitGLView", this,"OnClicked(TObject*)");
   
   fViewer[1] = new TEveViewer("SplitGLViewer[1]");
   fViewer[1]->SetGLViewer(fViewer1, fViewer1->GetFrame());
   fViewer[1]->IncDenyDestroy();
   if (fIsEmbedded && gEve) {
      fRhoZMgr->ImportElements((TEveElement *)gEve->GetGlobalScene());
      fRhoZMgr->ImportElements((TEveElement *)gEve->GetEventScene());
      fViewer[1]->AddScene(s);
      gEve->GetViewers()->AddElement(fViewer[1]);
      gRhoZMgr = fRhoZMgr;

      s = gEve->SpawnNewScene("R-Phi Projection");
      // projections
      fRPhiMgr = new TEveProjectionManager(TEveProjection::kPT_RPhi);
      s->AddElement(fRPhiMgr);
      gEve->AddToListTree(fRPhiMgr, kTRUE);
      TEveProjectionAxes* a = new TEveProjectionAxes(fRPhiMgr);
      s->AddElement(a);
   }

   // get bottom center split frame
   frm = fSplitFrame->GetSecond()->GetSecond()->GetFirst();
   frm->SetName("Bottom_Center");

   // create (embed) a GL viewer inside
   fViewer2 = new TGLEmbeddedViewer(frm, fPad);
   but5 = new TGLOverlayButton(fViewer2, "Swap", 10.0, -10.0, 55.0, 16.0);
   but5->Connect("Clicked(TGLViewerBase*)", "SplitGLView", this, "SwapToMainView(TGLViewerBase*)");
   but6 = new TGLOverlayButton(fViewer2, "Undock", 70.0, -10.0, 55.0, 16.0);
   but6->Connect("Clicked(TGLViewerBase*)", "SplitGLView", this, "UnDock(TGLViewerBase*)");
   frm->AddFrame(fViewer2->GetFrame(), new TGLayoutHints(kLHintsExpandX |
                  kLHintsExpandY));

   // set the camera to orthographic (XOY) for this viewer
   fViewer2->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
   // connect signal we are interested to
   fViewer2->Connect("MouseOver(TGLPhysicalShape*)", "SplitGLView", this,
                      "OnMouseOver(TGLPhysicalShape*)");
   fViewer2->Connect("Activated()", "SplitGLView", this,
                      "OnViewerActivated()");
   fViewer2->Connect("MouseIdle(TGLPhysicalShape*,UInt_t,UInt_t)",
                      "SplitGLView", this,
                      "OnMouseIdle(TGLPhysicalShape*,UInt_t,UInt_t)");
   fViewer2->Connect("Clicked(TObject*)", "SplitGLView", this,
                      "OnClicked(TObject*)");
   fViewer[2] = new TEveViewer("SplitGLViewer[2]");
   fViewer[2]->SetGLViewer(fViewer2, fViewer2->GetFrame());
   fViewer[2]->IncDenyDestroy();
   if (fIsEmbedded && gEve) {
      fRPhiMgr->ImportElements((TEveElement *)gEve->GetGlobalScene());
      fRPhiMgr->ImportElements((TEveElement *)gEve->GetEventScene());
      fViewer[2]->AddScene(s);
      gEve->GetViewers()->AddElement(fViewer[2]);
      gRPhiMgr = fRPhiMgr;
   }

   // get bottom right split frame
   frm = fSplitFrame->GetSecond()->GetSecond()->GetSecond();
   frm->SetName("Bottom_Right");

   dfrm = new TGDockableFrame(frm);
   dfrm->SetFixedSize(kFALSE);
   dfrm->EnableHide(kFALSE);
   hfrm = new TGHorizontalFrame(dfrm);
   button= new TGPictureButton(hfrm, gClient->GetPicture("swap.png"));
   button->SetToolTipText("Swap to big view");
   hfrm->AddFrame(button);
   button->Connect("Clicked()","SplitGLView",this,"SwapToMainView(TGLViewerBase*=0)");
   fgHtmlSummary = new HtmlSummary("Alice Event Display Summary Table");
   fgHtml = new TGHtml(hfrm, 100, 100, -1);
   hfrm->AddFrame(fgHtml, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
   dfrm->AddFrame(hfrm, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
   frm->AddFrame(dfrm, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

   if (fIsEmbedded && gEve) {
      gEve->GetListTree()->Connect("Clicked(TGListTreeItem*, Int_t, Int_t, Int_t)",
         "SplitGLView", this, "ItemClicked(TGListTreeItem*, Int_t, Int_t, Int_t)");
   }

   fShapedToolTip = new TGShapedToolTip("Default.png", 120, 22, 160, 110,
                                        23, 115, 12, "#ffff80");
   Resize(GetDefaultSize());
   MapSubwindows();
   MapWindow();
   LoadConfig(".everc");
}

//______________________________________________________________________________
SplitGLView::~SplitGLView()
{
   // Clean up main frame...
   //Cleanup();

   fMenuFile->Disconnect  ("Activated(Int_t)", this, "HandleMenu(Int_t)");
   fMenuCamera->Disconnect("Activated(Int_t)", this, "HandleMenu(Int_t)");
   fMenuScene->Disconnect ("Activated(Int_t)", this, "HandleMenu(Int_t)");
   fMenuHelp->Disconnect  ("Activated(Int_t)", this, "HandleMenu(Int_t)");
   
   fViewer0->Disconnect   ("MouseOver(TGLPhysicalShape*)", this,"OnMouseOver(TGLPhysicalShape*)");
   fViewer0->Disconnect   ("Activated()"     , this, "OnViewerActivated()");
   fViewer0->Disconnect   ("MouseIdle(TGLPhysicalShape*,UInt_t,UInt_t)",this, "OnMouseIdle(TGLPhysicalShape*,UInt_t,UInt_t)");

   fViewer1->Disconnect("MouseOver(TGLPhysicalShape*)", this,"OnMouseOver(TGLPhysicalShape*)");
   fViewer1->Disconnect("Activated()", this, "OnViewerActivated()");
   fViewer1->Disconnect("MouseIdle(TGLPhysicalShape*,UInt_t,UInt_t)",this, "OnMouseIdle(TGLPhysicalShape*,UInt_t,UInt_t)");

   fViewer2->Disconnect("MouseOver(TGLPhysicalShape*)", this,"OnMouseOver(TGLPhysicalShape*)");
   fViewer2->Disconnect("Activated()", this, "OnViewerActivated()");
   fViewer2->Disconnect("MouseIdle(TGLPhysicalShape*,UInt_t,UInt_t)",this, "OnMouseIdle(TGLPhysicalShape*,UInt_t,UInt_t)");

   if (!fIsEmbedded) {
      delete fViewer[0];
      delete fViewer[1];
      delete fViewer[2];
   }
   delete fShapedToolTip;
   delete fMenuFile;
   delete fMenuScene;
   delete fMenuCamera;
   delete fMenuHelp;
   if (!fIsEmbedded)
      delete fMenuBar;
   delete fViewer0;
   delete fViewer1;
   delete fViewer2;
   delete fSplitFrame;
   delete fPad;
   if (!fIsEmbedded) {
      delete fStatusBar;
      gApplication->Terminate(0);
   }
}

//______________________________________________________________________________
void SplitGLView::HandleMenu(Int_t id)
{
   // Handle menu items.

   static TString rcdir(".");
   static TString rcfile(".everc");

   switch (id) {

      case kFileOpen:
         {
            static TString dir(".");
            TGFileInfo fi;
            fi.fFileTypes = filetypes;
            fi.fIniDir    = StrDup(dir);
            new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fi);
            if (fi.fFilename)
               OpenFile(fi.fFilename);
            dir = fi.fIniDir;
         }
         break;

      case kFileLoadConfig:
         {
            TGFileInfo fi;
            fi.fFileTypes = rcfiletypes;
            fi.fIniDir    = StrDup(rcdir);
            fi.fFilename  = StrDup(rcfile);
            new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fi);
            if (fi.fFilename) {
               rcfile = fi.fFilename;
               LoadConfig(fi.fFilename);
            }
            rcdir = fi.fIniDir;
         }
         break;

      case kFileSaveConfig:
         {
            TGFileInfo fi;
            fi.fFileTypes = rcfiletypes;
            fi.fIniDir    = StrDup(rcdir);
            fi.fFilename  = StrDup(rcfile);
            new TGFileDialog(gClient->GetRoot(), this, kFDSave, &fi);
            if (fi.fFilename) {
               rcfile = fi.fFilename;
               SaveConfig(fi.fFilename);
            }
            rcdir = fi.fIniDir;
         }
         break;

      case kFileExit:
         CloseWindow();
         break;

      case kGLPerspYOZ:
         if (fActViewer)
            fActViewer->SetCurrentCamera(TGLViewer::kCameraPerspYOZ);
         break;
      case kGLPerspXOZ:
         if (fActViewer)
            fActViewer->SetCurrentCamera(TGLViewer::kCameraPerspXOZ);
         break;
      case kGLPerspXOY:
         if (fActViewer)
            fActViewer->SetCurrentCamera(TGLViewer::kCameraPerspXOY);
         break;
      case kGLXOY:
         if (fActViewer)
            fActViewer->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
         break;
      case kGLXOZ:
         if (fActViewer)
            fActViewer->SetCurrentCamera(TGLViewer::kCameraOrthoXOZ);
         break;
      case kGLZOY:
         if (fActViewer)
            fActViewer->SetCurrentCamera(TGLViewer::kCameraOrthoZOY);
         break;
      case kGLOrthoRotate:
         ToggleOrthoRotate();
         break;
      case kGLOrthoDolly:
         ToggleOrthoDolly();
         break;

      case kSceneUpdate:
         if (fActViewer)
            fActViewer->UpdateScene();
         UpdateSummary();
         break;

      case kSceneUpdateAll:
         fViewer0->UpdateScene();
         fViewer1->UpdateScene();
         fViewer2->UpdateScene();
         UpdateSummary();
         break;

      case kSummaryUpdate:
         UpdateSummary();
         break;

      case kHelpAbout:
         {
#ifdef R__UNIX
            TString rootx;
# ifdef ROOTBINDIR
            rootx = ROOTBINDIR;
# else
            rootx = gSystem->Getenv("ROOTSYS");
            if (!rootx.IsNull()) rootx += "/bin";
# endif
            rootx += "/root -a &";
            gSystem->Exec(rootx);
#else
#ifdef WIN32
            new TWin32SplashThread(kTRUE);
#else
            char str[32];
            sprintf(str, "About ROOT %s...", gROOT->GetVersion());
            hd = new TRootHelpDialog(this, str, 600, 400);
            hd->SetText(gHelpAbout);
            hd->Popup();
#endif
#endif
         }
         break;

      default:
         break;
   }
}

//______________________________________________________________________________
void SplitGLView::OnClicked(TObject *obj)
{
   // Handle click events in GL viewer

   if (obj)
      fStatusBar->SetText(Form("User clicked on: \"%s\"", obj->GetName()), 1);
   else
      fStatusBar->SetText("", 1);
}

//______________________________________________________________________________
void SplitGLView::OnMouseIdle(TGLPhysicalShape *shape, UInt_t posx, UInt_t posy)
{
   // Slot used to handle "OnMouseIdle" signal coming from any GL viewer.
   // We receive a pointer on the physical shape in which the mouse cursor is
   // and the actual cursor position (x,y)

   Window_t wtarget;
   Int_t    x = 0, y = 0;

   static TH1F *h1f = 0;
   //   TFormula *form1 = new TFormula("form1","abs(sin(x)/x)");
   //   form1->Update(); // silent warning about unused variable...
   TF1 *sqroot = new TF1("sqroot","x*gaus(0) + [3]*form1",0,10);
   sqroot->SetParameters(10,4,1,20);
   if (h1f == 0)
      h1f = new TH1F("h1f","",50,0,10);
   h1f->Reset();
   h1f->SetFillColor(45);
   h1f->SetStats(0);
   h1f->FillRandom("sqroot",200);

   if (fShapedToolTip) {
      fShapedToolTip->UnmapWindow();
   }
   if (shape && shape->GetLogical() && shape->GetLogical()->GetExternal()) {
      // get the actual viewer who actually emitted the signal
      TGLEmbeddedViewer *actViewer = dynamic_cast<TGLEmbeddedViewer*>((TQObject*)gTQSender);
      // then translate coordinates from the root (screen) coordinates
      // to the actual frame (viewer) ones
      gVirtualX->TranslateCoordinates(actViewer->GetFrame()->GetId(),
               gClient->GetDefaultRoot()->GetId(), posx, posy, x, y,
               wtarget);
      // Then display our tooltip at this x,y location
      if (fShapedToolTip) {
         fShapedToolTip->Show(x+5, y+5, Form("%s\n     \n%s",
                              shape->GetLogical()->GetExternal()->IsA()->GetName(),
                              shape->GetLogical()->GetExternal()->GetName()), h1f);
      }
   }
}

//______________________________________________________________________________
void SplitGLView::OnMouseOver(TGLPhysicalShape *shape)
{
   // Slot used to handle "OnMouseOver" signal coming from any GL viewer.
   // We receive a pointer on the physical shape in which the mouse cursor is.

   // display information on the physical shape in the status bar
   if (shape && shape->GetLogical() && shape->GetLogical()->GetExternal())
      fStatusBar->SetText(Form("Mouse Over: \"%s\"",
         shape->GetLogical()->GetExternal()->GetName()), 0);
   else
      fStatusBar->SetText("", 0);
}

//______________________________________________________________________________
void SplitGLView::OnViewerActivated()
{
   // Slot used to handle "Activated" signal coming from any GL viewer.
   // Used to know which GL viewer is active.

   static Pixel_t green = 0;
   // set the actual GL viewer frame to default color
   if (fActViewer && fActViewer->GetFrame())
      fActViewer->GetFrame()->ChangeBackground(GetDefaultFrameBackground());

   // change the actual GL viewer to the one who emitted the signal
   // fActViewer = (TGLEmbeddedViewer *)gTQSender;
   fActViewer = dynamic_cast<TGLEmbeddedViewer*>((TQObject*)gTQSender);

   if (fActViewer == 0) {
      printf ("dyncast failed ...\n");
      return;
   }

   // get the highlight color (only once)
   if (green == 0) {
      gClient->GetColorByName("green", green);
   }
   // set the new actual GL viewer frame to highlight color
   if (fActViewer->GetFrame())
      fActViewer->GetFrame()->ChangeBackground(green);

   // update menu entries to match actual viewer's options
   if (fActViewer->GetOrthoXOYCamera()->GetDollyToZoom() &&
       fActViewer->GetOrthoXOZCamera()->GetDollyToZoom() &&
       fActViewer->GetOrthoZOYCamera()->GetDollyToZoom())
      fMenuCamera->UnCheckEntry(kGLOrthoDolly);
   else
      fMenuCamera->CheckEntry(kGLOrthoDolly);

   if (fActViewer->GetOrthoXOYCamera()->GetEnableRotate() &&
       fActViewer->GetOrthoXOYCamera()->GetEnableRotate() &&
       fActViewer->GetOrthoXOYCamera()->GetEnableRotate())
      fMenuCamera->CheckEntry(kGLOrthoRotate);
   else
      fMenuCamera->UnCheckEntry(kGLOrthoRotate);
}

//______________________________________________________________________________
void SplitGLView::OpenFile(const char *fname)
{
   // Open a Root file to display a geometry in the GL viewers.

   TString filename = fname;
   // check if the file type is correct
   if (!filename.EndsWith(".root")) {
      new TGMsgBox(gClient->GetRoot(), this, "OpenFile",
                   Form("The file \"%s\" is not a root file!", fname),
                   kMBIconExclamation, kMBOk);
      return;
   }
   // check if the root file contains a geometry
   if (TGeoManager::Import(fname) == 0) {
      new TGMsgBox(gClient->GetRoot(), this, "OpenFile",
                   Form("The file \"%s\" does't contain a geometry", fname),
                   kMBIconExclamation, kMBOk);
      return;
   }
   gGeoManager->DefaultColors();
   // delete previous primitives (if any)
   fPad->GetListOfPrimitives()->Delete();
   // and add the geometry to eve pad (container)
   fPad->GetListOfPrimitives()->Add(gGeoManager->GetTopVolume());
   // paint the geometry in each GL viewer
   fViewer0->PadPaint(fPad);
   fViewer1->PadPaint(fPad);
   fViewer2->PadPaint(fPad);
}

//______________________________________________________________________________
void SplitGLView::ToggleOrthoRotate()
{
   // Toggle state of the 'Ortho allow rotate' menu entry.

   if (fMenuCamera->IsEntryChecked(kGLOrthoRotate))
      fMenuCamera->UnCheckEntry(kGLOrthoRotate);
   else
      fMenuCamera->CheckEntry(kGLOrthoRotate);
   Bool_t state = fMenuCamera->IsEntryChecked(kGLOrthoRotate);
   if (fActViewer) {
      fActViewer->GetOrthoXOYCamera()->SetEnableRotate(state);
      fActViewer->GetOrthoXOYCamera()->SetEnableRotate(state);
      fActViewer->GetOrthoXOYCamera()->SetEnableRotate(state);
   }
}

//______________________________________________________________________________
void SplitGLView::ToggleOrthoDolly()
{
   // Toggle state of the 'Ortho allow dolly' menu entry.

   if (fMenuCamera->IsEntryChecked(kGLOrthoDolly))
      fMenuCamera->UnCheckEntry(kGLOrthoDolly);
   else
      fMenuCamera->CheckEntry(kGLOrthoDolly);
   Bool_t state = ! fMenuCamera->IsEntryChecked(kGLOrthoDolly);
   if (fActViewer) {
      fActViewer->GetOrthoXOYCamera()->SetDollyToZoom(state);
      fActViewer->GetOrthoXOZCamera()->SetDollyToZoom(state);
      fActViewer->GetOrthoZOYCamera()->SetDollyToZoom(state);
   }
}

//______________________________________________________________________________
void SplitGLView::ItemClicked(TGListTreeItem *item, Int_t, Int_t, Int_t)
{
   // Item has been clicked, based on mouse button do:

   static const TEveException eh("SplitGLView::ItemClicked ");
   TEveElement* re = (TEveElement*)item->GetUserData();
   if(re == 0) return;
   TObject* obj = re->GetObject(eh);
   if (obj->InheritsFrom("TEveViewer")) {
      TGLViewer *v = ((TEveViewer *)obj)->GetGLViewer();
      //v->Activated();
      if (v->InheritsFrom("TGLEmbeddedViewer")) {
         TGLEmbeddedViewer *ev = (TGLEmbeddedViewer *)v;
         gVirtualX->SetInputFocus(ev->GetGLWidget()->GetId());
      }
   }
}

//______________________________________________________________________________
void SplitGLView::LoadConfig(const char *fname)
{

   Int_t height, width;
   TEnv *env = new TEnv(fname);

   Int_t mainheight = env->GetValue("MainView.Height", 434);
   Int_t blwidth    = env->GetValue("Bottom.Left.Width", 266);
   Int_t bcwidth    = env->GetValue("Bottom.Center.Width", 266);
   Int_t brwidth    = env->GetValue("Bottom.Right.Width", 266);
   Int_t top_height = env->GetValue("Right.Tab.Height", 0);
   Int_t bottom_height = env->GetValue("Bottom.Tab.Height", 0);

   if (fIsEmbedded && gEve) {
      Int_t sel = env->GetValue("Eve.Selection", gEve->GetSelection()->GetPickToSelect());
      Int_t hi = env->GetValue("Eve.Highlight", gEve->GetHighlight()->GetPickToSelect());
      gEve->GetBrowser()->EveMenu(9+sel);
      gEve->GetBrowser()->EveMenu(13+hi);

      width  = env->GetValue("Eve.Width", (Int_t)gEve->GetBrowser()->GetWidth());
      height = env->GetValue("Eve.Height", (Int_t)gEve->GetBrowser()->GetHeight());
      gEve->GetBrowser()->Resize(width, height);
   }

   // top (main) split frame
   width = fSplitFrame->GetFirst()->GetWidth();
   fSplitFrame->GetFirst()->Resize(width, mainheight);
   // bottom left split frame
   height = fSplitFrame->GetSecond()->GetFirst()->GetHeight();
   fSplitFrame->GetSecond()->GetFirst()->Resize(blwidth, height);
   // bottom center split frame
   height = fSplitFrame->GetSecond()->GetSecond()->GetFirst()->GetHeight();
   fSplitFrame->GetSecond()->GetSecond()->GetFirst()->Resize(bcwidth, height);
   // bottom right split frame
   height = fSplitFrame->GetSecond()->GetSecond()->GetSecond()->GetHeight();
   fSplitFrame->GetSecond()->GetSecond()->GetSecond()->Resize(brwidth, height);

   fSplitFrame->Layout();

   if (fIsEmbedded && gEve) {
      width = ((TGCompositeFrame *)gEve->GetBrowser()->GetTabBottom()->GetParent())->GetWidth();
      ((TGCompositeFrame *)gEve->GetBrowser()->GetTabBottom()->GetParent())->Resize(width, bottom_height);
      width = ((TGCompositeFrame *)gEve->GetBrowser()->GetTabRight()->GetParent())->GetWidth();
      ((TGCompositeFrame *)gEve->GetBrowser()->GetTabRight()->GetParent())->Resize(width, top_height);
   }
}

//______________________________________________________________________________
void SplitGLView::SaveConfig(const char *fname)
{

   Int_t bottom_height = 0;
   Int_t top_height = 0;
   TGSplitFrame *frm;
   TEnv *env = new TEnv(fname);

   if (fIsEmbedded && gEve) {
      env->SetValue("Eve.Width", (Int_t)gEve->GetBrowser()->GetWidth());
      env->SetValue("Eve.Height", (Int_t)gEve->GetBrowser()->GetHeight());
   }
   // get top (main) split frame
   frm = fSplitFrame->GetFirst();
   env->SetValue("MainView.Height", (Int_t)frm->GetHeight());
   // get bottom left split frame
   frm = fSplitFrame->GetSecond()->GetFirst();
   env->SetValue("Bottom.Left.Width", (Int_t)frm->GetWidth());
   // get bottom center split frame
   frm = fSplitFrame->GetSecond()->GetSecond()->GetFirst();
   env->SetValue("Bottom.Center.Width", (Int_t)frm->GetWidth());
   // get bottom right split frame
   frm = fSplitFrame->GetSecond()->GetSecond()->GetSecond();
   env->SetValue("Bottom.Right.Width", (Int_t)frm->GetWidth());
   if (fIsEmbedded && gEve) {
      top_height = (Int_t)((TGCompositeFrame *)gEve->GetBrowser()->GetTabRight()->GetParent())->GetHeight();
      env->SetValue("Right.Tab.Height", top_height);
      bottom_height = (Int_t)((TGCompositeFrame *)gEve->GetBrowser()->GetTabBottom()->GetParent())->GetHeight();
      env->SetValue("Bottom.Tab.Height", bottom_height);

      env->SetValue("Eve.Selection", gEve->GetSelection()->GetPickToSelect());
      env->SetValue("Eve.Highlight", gEve->GetHighlight()->GetPickToSelect());
   }

   env->SaveLevel(kEnvLocal);
#ifdef R__WIN32
   if (!gSystem->AccessPathName(Form("%s.new", fname))) {
      gSystem->Exec(Form("del %s", fname));
      gSystem->Rename(Form("%s.new", fname), fname);
   }
#endif
}

//______________________________________________________________________________
void SplitGLView::SwapToMainView(TGLViewerBase *viewer)
{
   // Swap frame embedded in a splitframe to the main view (slot method).

   TGCompositeFrame *parent = 0;
   if (!fSplitFrame->GetFirst()->GetFrame())
      return;
   if (viewer == 0) {
      TGPictureButton *src = (TGPictureButton*)gTQSender;
      parent = (TGCompositeFrame *)src->GetParent();
      while (parent && !parent->InheritsFrom("TGSplitFrame")) {
         parent = (TGCompositeFrame *)parent->GetParent();
      }
   }
   else {
      TGCompositeFrame *src = ((TGLEmbeddedViewer *)viewer)->GetFrame();
      if (!src) return;
      TGLOverlayButton *but = (TGLOverlayButton *)((TQObject *)gTQSender);
      but->ResetState();
      parent = (TGCompositeFrame *)src->GetParent();
   }
   if (parent && parent->InheritsFrom("TGSplitFrame"))
      ((TGSplitFrame *)parent)->SwitchToMain();
}

//______________________________________________________________________________
void SplitGLView::UnDock(TGLViewerBase *viewer)
{
   // Undock frame embedded in a splitframe (slot method).

   TGCompositeFrame *src = ((TGLEmbeddedViewer *)viewer)->GetFrame();
   if (!src) return;
   TGLOverlayButton *but = (TGLOverlayButton *)((TQObject *)gTQSender);
   but->ResetState();
   TGCompositeFrame *parent = (TGCompositeFrame *)src->GetParent();
   if (parent && parent->InheritsFrom("TGSplitFrame"))
      ((TGSplitFrame *)parent)->ExtractFrame();
}

//______________________________________________________________________________
void SplitGLView::UpdateSummary()
{
   // Update summary of current event.

   TEveElement::List_i i;
   TEveElement::List_i j;
   Int_t k;
   TEveElement *el;
   HtmlObjTable *table;
   TEveEventManager *mgr = gEve ? gEve->GetCurrentEvent() : 0;
   if (mgr) {
      fgHtmlSummary->Clear("D");
      for (i=mgr->BeginChildren(); i!=mgr->EndChildren(); ++i) {
         el = ((TEveElement*)(*i));
         if (el->IsA() == TEvePointSet::Class()) {
            TEvePointSet *ps = (TEvePointSet *)el;
            TString ename  = ps->GetElementName();
            TString etitle = ps->GetElementTitle();
            if (ename.First('\'') != kNPOS)
               ename.Remove(ename.First('\''));
            etitle.Remove(0, 2);
            Int_t nel = atoi(etitle.Data());
            table = fgHtmlSummary->AddTable(ename, 0, nel);
         }
         else if (el->IsA() == TEveTrackList::Class()) {
            TEveTrackList *tracks = (TEveTrackList *)el;
            TString ename  = tracks->GetElementName();
            if (ename.First('\'') != kNPOS)
               ename.Remove(ename.First('\''));
            table = fgHtmlSummary->AddTable(ename.Data(), 5,
                     tracks->NumChildren(), kTRUE, "first");
            table->SetLabel(0, "Momentum");
            table->SetLabel(1, "P_t");
            table->SetLabel(2, "Phi");
            table->SetLabel(3, "Theta");
            table->SetLabel(4, "Eta");
            k=0;
            for (j=tracks->BeginChildren(); j!=tracks->EndChildren(); ++j) {
               Float_t p     = ((TEveTrack*)(*j))->GetMomentum().Mag();
               table->SetValue(0, k, p);
               Float_t pt    = ((TEveTrack*)(*j))->GetMomentum().Perp();
               table->SetValue(1, k, pt);
               Float_t phi   = ((TEveTrack*)(*j))->GetMomentum().Phi();
               table->SetValue(2, k, phi);
               Float_t theta = ((TEveTrack*)(*j))->GetMomentum().Theta();
               table->SetValue(3, k, theta);
               Float_t eta   = ((TEveTrack*)(*j))->GetMomentum().Eta();
               table->SetValue(4, k, eta);
               ++k;
            }
         }
      }
      fgHtmlSummary->Build();
      fgHtml->Clear();
      fgHtml->ParseText((char*)fgHtmlSummary->Html().Data());
      fgHtml->Layout();
   }
}
