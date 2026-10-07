#

#include "TPad.h"
#include "TGLSAViewer.h"

#include "murat/gui/TEvdCaloView.hh"

// ClassImp(TEvdCaloView)

namespace murat {

//-----------------------------------------------------------------------------
  TEvdCaloView::TEvdCaloView(const char* Name, int Type, int Index):
    TEvdView(Name,Type,Index)
  {
  }

//-----------------------------------------------------------------------------  
  TEvdCaloView::TEvdCaloView(const char* Name, int Type, int Index, const char* Title):
    TEvdView(Name,Type,Index,Title)
  {
  }

//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
  TEvdCaloView::~TEvdCaloView() {
  }

//-----------------------------------------------------------------------------
// XY projection
//-----------------------------------------------------------------------------
  void TEvdCaloView::Draw(Option_t* Opt) {

    int nnodes = GetNNodes();
    for (int i=0; i<nnodes; i++) {
      TEvdVisNode* node = GetNode(i);
      node->Draw(Opt);
    }
                                        // do it only once
    if (fCanvas == nullptr) {
      fCanvas = GetCurrentOpenGLCanvas();
    }
                                        // set XY projection
    
    auto viewer = (TGLSAViewer*) gPad->GetViewer3D();
    viewer->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);  // XY projection
    viewer->ResetCurrentCamera();
    viewer->RequestDraw();
  }


}
