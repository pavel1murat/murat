#


#include "murat/gui/TEvdView.hh"


// ClassImp(TEvdView)

namespace murat {
  
  TEvdView::TEvdView(int Type, int Index) {
    fListOfNodes = new TObjArray();
  }

  TEvdView::TEvdView(int Type, int Index, const char* Name, const char* Title) {
    fListOfNodes = new TObjArray();
  }

  TEvdView::~TEvdView() {
  }

}
