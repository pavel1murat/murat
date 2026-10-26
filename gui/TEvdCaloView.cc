#

#include "TPad.h"
#include "TGLSAViewer.h"

#include "murat/gui/TEvdCaloView.hh"

// ClassImp(TEvdCaloView)

namespace murat {

//-----------------------------------------------------------------------------
  TEvdCaloView::TEvdCaloView(const char* Name, TGLViewer::ECameraType Type, int Index):
    TEvdView(Name,Type,Index)
  {
    fViewer->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);  // XY projection
  }

//-----------------------------------------------------------------------------  
  TEvdCaloView::TEvdCaloView(const char* Name, TGLViewer::ECameraType Type, int Index, const char* Title):
    TEvdView(Name,Type,Index,Title)
  {
    fViewer->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);  // XY projection
  }

//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
  TEvdCaloView::~TEvdCaloView() {
  }

//-----------------------------------------------------------------------------
// XY projection
//-----------------------------------------------------------------------------
  // void TEvdCaloView::Draw(Option_t* Opt) {

  //   int nnodes = GetNNodes();
  //   for (int i=0; i<nnodes; i++) {
  //     TEvdVisNode* node = GetNode(i);
  //     node->Draw(Opt);
  //   }
  //                                       // do it only once
  //   auto viewer = (TGLSAViewer*) gPad->GetViewer3D();
  //   if (fCanvas == nullptr) {
  //     fCanvas = GetCurrentOpenGLCanvas();
  //     viewer->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);  // XY projection
  //   //    viewer->ResetCurrentCamera();
  //   }
  //                                       // set XY projection
    
  //   viewer->RequestDraw();
  // }


}
