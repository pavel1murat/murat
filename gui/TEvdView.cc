#


#include "murat/gui/TEvdView.hh"
#include "TPad.h"
#include "TGLSAViewer.h"


// ClassImp(TEvdView)

namespace murat {
  
//-----------------------------------------------------------------------------
  TEvdView::TEvdView(const char* Name, int Type, int Index) : TNamed(Name,Name) {
    fListOfNodes = new TObjArray();
    fCanvas      = nullptr;
  }

//-----------------------------------------------------------------------------
  TEvdView::TEvdView(const char* Name, int Type, int Index, const char* Title) :
    TNamed(Name,Title)
  {
    fListOfNodes = new TObjArray();
    fCanvas      = nullptr;
  }

//-----------------------------------------------------------------------------
  TEvdView::~TEvdView() {
  }

//-----------------------------------------------------------------------------
  void TEvdView::Draw(Option_t* Opt) {

    int nnodes = GetNNodes();
    for (int i=0; i<nnodes; i++) {
      TEvdVisNode* node = GetNode(i);
      node->Draw(Opt);
    }
                                        // do it only once
    if (fCanvas == nullptr) {
      fCanvas = GetCurrentOpenGLCanvas();
      // the default name is always 'c1', to avoid auto-deletion
      fCanvas->SetName(Form("c_%s",GetName()));
    }

    auto viewer = (TGLSAViewer*) gPad->GetViewer3D();
    // viewer->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);  // XY projection
    // viewer->ResetCurrentCamera();
    viewer->RequestDraw();
  }

//-----------------------------------------------------------------------------
  TCanvas* TEvdView::GetCurrentOpenGLCanvas() {
    if (!gPad) {
      return nullptr;
    }

    return gPad->GetCanvas();
  }
}

