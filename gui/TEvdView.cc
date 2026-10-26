///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <format>


#include "murat/gui/TEvdView.hh"
#include "TROOT.h"
// #include "TPad.h"
#include "TGLViewer.h"


// ClassImp(TEvdView)

namespace murat {
  
//-----------------------------------------------------------------------------
// different views should have different names, rely on that
//-----------------------------------------------------------------------------
  TEvdView::TEvdView(const char* Name, TGLViewer::ECameraType Type, int Index) :
    TNamed(Name,Name),
    fOpen(0),
    fCameraType(Type)
  {
    fListOfNodes = new TObjArray();
    gROOT->SetBatch(1);
    fCanvas = new TCanvas(Form("c_%s",Name));
    gROOT->SetBatch(0);
    fViewer = nullptr;

    // fCanvas->cd();
    // fCanvas->Update();

    // auto* v3 = gPad->GetViewer3D();
    // fViewer = dynamic_cast<TGLViewer*>(v3);

    // if (fViewer == nullptr) {
    //   std::cout << std::format("ERROR: fViewer==nullptr\n");
    //   // fViewer->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
    // }
  }

//-----------------------------------------------------------------------------
  TEvdView::TEvdView(const char* Name, TGLViewer::ECameraType Type, int Index, const char* Title) :
    TNamed(Name,Title),
    fOpen(0),
    fCameraType(Type)
  {
    fListOfNodes = new TObjArray();

    gROOT->SetBatch(1);
    fCanvas = new TCanvas(Form("c_%s",Name));
    gROOT->SetBatch(0);

    fViewer = nullptr;
    // if (fViewer == nullptr) {
    //   std::cout << std::format("ERROR: fViewer==nullptr\n");

    // //   fViewer->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
    // }
  }

//-----------------------------------------------------------------------------
  TEvdView::~TEvdView() {
  }

//-----------------------------------------------------------------------------
  void TEvdView::Draw(Option_t* Opt) {

    fCanvas->cd();

    int nnodes = GetNNodes();
    for (int i=0; i<nnodes; i++) {
      TEvdVisNode* node = GetNode(i);
      node->Draw(Opt);
    }

    gPad->Update();

    if (fViewer == nullptr) {
      TVirtualViewer3D* v3 = gPad->GetViewer3D();
      fViewer = dynamic_cast<TGLViewer*>(v3);
      fViewer->SetCurrentCamera(fCameraType);
    }

    fViewer->RequestDraw();
  }

//-----------------------------------------------------------------------------
// a view has to be opened explicitly, not open by default
//-----------------------------------------------------------------------------
  int TEvdView::Open() {

    // do it only once, there should ba a pad associated with the view
    if (fOpen == 1) return 0;
    fOpen = 1;

    fCanvas->cd();

    if (fViewer == nullptr) { 
      auto* v3 = gPad->GetViewer3D();
      fViewer = dynamic_cast<TGLViewer*>(v3);
      fViewer->SetCurrentCamera(fCameraType);
    }

    fViewer->RequestDraw();
    return 0;
  }

//-----------------------------------------------------------------------------
  void TEvdView::Update() {

    //     if (fOpen == 0) return;
    
    fCanvas->Modified();
    fCanvas->Update();
  }
}
