#

#include "TPad.h"
#include "TGLSAViewer.h"

#include "murat/gui/TEvdCrvView.hh"

// ClassImp(TEvdCrvView)

namespace murat {

//-----------------------------------------------------------------------------
  TEvdCrvView::TEvdCrvView(const char* Name, TGLViewer::ECameraType Type, int Index):
    TEvdView(Name,Type,Index)
  {
  }

//-----------------------------------------------------------------------------  
  TEvdCrvView::TEvdCrvView(const char* Name, TGLViewer::ECameraType Type, int Index, const char* Title):
    TEvdView(Name,Type,Index,Title)
  {
  }

//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
  TEvdCrvView::~TEvdCrvView() {
  }

// //-----------------------------------------------------------------------------
// // XY projection
// //-----------------------------------------------------------------------------
//   void TEvdCrvView::Draw(Option_t* Opt) {

//     int nnodes = GetNNodes();
//     for (int i=0; i<nnodes; i++) {
//       TEvdVisNode* node = GetNode(i);
//       node->Draw(Opt);
//     }
//                                         // do it only once

//     auto viewer = (TGLSAViewer*) gPad->GetViewer3D();
//     if (fCanvas == nullptr) {
//       fCanvas = GetCurrentOpenGLCanvas();
//       viewer->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);  // XY projection
//       viewer->ResetCurrentCamera();
//     }
//                                         // set XY projection
//     viewer->RequestDraw();
//   }


}
