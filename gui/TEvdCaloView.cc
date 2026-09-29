#

#include "TPad.h"
#include "TGLSAViewer.h"

#include "murat/gui/TEvdCaloView.hh"

// ClassImp(TEvdCaloView)

namespace murat {

//-----------------------------------------------------------------------------
  TEvdCaloView::TEvdCaloView(int Type, int Index) {
  }

//-----------------------------------------------------------------------------  
  TEvdCaloView::TEvdCaloView(int Type, int Index, const char* Name, const char* Title) {
  }

//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
  TEvdCaloView::~TEvdCaloView() {
  }

//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
  void TEvdCaloView::Draw(Option_t* Opt) {

    int nnodes = GetNNodes();
    for (int i=0; i<nnodes; i++) {
      TEvdVisNode* node = GetNode(i);
      node->Draw(Opt);
    }

    auto viewer = (TGLSAViewer*) gPad->GetViewer3D();
    viewer->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);  // XY projection
    viewer->ResetCurrentCamera();
    viewer->RequestDraw();

  }


}
