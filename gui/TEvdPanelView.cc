///////////////////////////////////////////////////////////////////////////////
// a view showing a single tracker panel, see murat/gui/gui/TEvdPanelView.hh
///////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <string>

#include "murat/gui/TEvdPanelView.hh"

namespace murat {

//-----------------------------------------------------------------------------
TEvdPanelView::TEvdPanelView(const char* Name, TGLViewer::ECameraType Type) :
  TEvdView(Name,Type),
  fPlane(-1),
  fPanel(-1)
{
  fPanelNode = new TPanelVisNode((std::string(Name)+"_node").c_str());
  AddNode(fPanelNode);
}

//-----------------------------------------------------------------------------
TEvdPanelView::~TEvdPanelView() {
}

//-----------------------------------------------------------------------------
int TEvdPanelView::SetPanel(int Plane, int Panel) {
  if ((Plane < 0) or (Plane >= kNPlanes)) return -1;
  if ((Panel < 0) or (Panel >= kNPanelsPerPlane)) return -1;

  fPlane = Plane;
  fPanel = Panel;
  fPanelNode->SetPanel(Plane,Panel);
  return 0;
}



}
