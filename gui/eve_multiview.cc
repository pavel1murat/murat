///////////////////////////////////////////////////////////////////////////////
// cloned from the original code tutorials/eve/MultiView.C by Matevz Tadel
///////////////////////////////////////////////////////////////////////////////
#include "murat/gui/eve_multiview.hh"

// MultiView
//
// Structure encapsulating standard views: 3D, r-phi and rho-z.
// Includes scenes and projection managers.
//
// Should be used in compiled mode.

//---------------------------------------------------------------------------
// Constructor --- creates required scenes, projection managers
// and GL viewers.
//-----------------------------------------------------------------------------
Mu2eMultiView::Mu2eMultiView() {
  
  // Scenes
  //========
  
  fRPhiGeomScene  = gEve->SpawnNewScene("RPhi Geometry"  ,"Scene holding projected geometry   for the RPhi view.");
  fRhoZGeomScene  = gEve->SpawnNewScene("RhoZ Geometry"  ,"Scene holding projected geometry   for the RhoZ view.");
  fRPhiEventScene = gEve->SpawnNewScene("RPhi Event Data","Scene holding projected event-data for the RPhi view.");
  fRhoZEventScene = gEve->SpawnNewScene("RhoZ Event Data","Scene holding projected event-data for the RhoZ view.");
  

  // Projection managers
  //=====================

  fRPhiMgr = new TEveProjectionManager(TEveProjection::kPT_RPhi);
  gEve->AddToListTree(fRPhiMgr, kFALSE);
  {
    TEveProjectionAxes* a = new TEveProjectionAxes(fRPhiMgr);
    a->SetMainColor(kWhite);
    a->SetTitle("R-Phi");
    a->SetTitleSize(0.05);
    a->SetTitleFont(102);
    a->SetLabelSize(0.025);
    a->SetLabelFont(102);
    fRPhiGeomScene->AddElement(a);
  }

  fRhoZMgr = new TEveProjectionManager(TEveProjection::kPT_RhoZ);
  gEve->AddToListTree(fRhoZMgr, kFALSE);
  {
    TEveProjectionAxes* a = new TEveProjectionAxes(fRhoZMgr);
    a->SetMainColor(kWhite);
    a->SetTitle("Rho-Z");
    a->SetTitleSize(0.05);
    a->SetTitleFont(102);
    a->SetLabelSize(0.025);
    a->SetLabelFont(102);
    fRhoZGeomScene->AddElement(a);
  }

  // Viewers
  //=========

  TEveWindowSlot *slot = 0;
  TEveWindowPack *pack = 0;

  slot = TEveWindow::CreateWindowInTab(gEve->GetBrowser()->GetTabRight());
  pack = slot->MakePack();
  pack->SetElementName("Multi View");
  pack->SetHorizontal();
  pack->SetShowTitleBar(kFALSE);
  pack->NewSlot()->MakeCurrent();
  f3DView = gEve->SpawnNewViewer("3D View", "");
  f3DView->AddScene(gEve->GetGlobalScene());
  f3DView->AddScene(gEve->GetEventScene());
    
  pack = pack->NewSlot()->MakePack();
  pack->SetShowTitleBar(kFALSE);
  pack->NewSlot()->MakeCurrent();
  fRPhiView = gEve->SpawnNewViewer("RPhi View", "");
  fRPhiView->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
  fRPhiView->AddScene(fRPhiGeomScene);
  fRPhiView->AddScene(fRPhiEventScene);

  pack->NewSlot()->MakeCurrent();
  fRhoZView = gEve->SpawnNewViewer("RhoZ View", "");
  fRhoZView->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
  fRhoZView->AddScene(fRhoZGeomScene);
  fRhoZView->AddScene(fRhoZEventScene);
}

//-----------------------------------------------------------------------------
// empty so far, bu this will change
//-----------------------------------------------------------------------------
Mu2eMultiView::~Mu2eMultiView() {
}

//-----------------------------------------------------------------------------
void Mu2eMultiView::SetDepth(float Depth) {
    // Set current depth on all projection managers.
  fRPhiMgr->SetCurrentDepth(Depth);
  fRhoZMgr->SetCurrentDepth(Depth);
}

//---------------------------------------------------------------------------
  
void Mu2eMultiView::ImportGeomRPhi(TEveElement* el) {
  fRPhiMgr->ImportElements(el, fRPhiGeomScene);
}

void Mu2eMultiView::ImportGeomRhoZ(TEveElement* el) {
  fRhoZMgr->ImportElements(el, fRhoZGeomScene);
}

void Mu2eMultiView::ImportEventRPhi(TEveElement* el) {
  fRPhiMgr->ImportElements(el, fRPhiEventScene);
}

void Mu2eMultiView::ImportEventRhoZ(TEveElement* el) {
  fRhoZMgr->ImportElements(el, fRhoZEventScene);
}

//---------------------------------------------------------------------------

void Mu2eMultiView::DestroyEventRPhi() {
  fRPhiEventScene->DestroyElements();
}

void Mu2eMultiView::DestroyEventRhoZ() {
  fRhoZEventScene->DestroyElements();
}
