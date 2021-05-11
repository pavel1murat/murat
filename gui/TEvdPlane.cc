///////////////////////////////////////////////////////////////////////////////

#include "murat/gui/TEvdPlane.hh"

ClassImp(murat::TEvdPlane)

namespace murat {
//-----------------------------------------------------------------------------
TEvdPlane::TEvdPlane(int I): TEveElementList(Form("plane_%02i",I),Form("plane_%02i",I)) {
  fNumber = I;
  for (int i=0; i<kNPanels; i++) {
    fPanel[i]         = new TEvdPanel(i);
    fPanel[i]->fPlane = this;
    AddElement(fPanel[i]);
  }
  SetRnrSelfChildren(false,true);
}  

}
