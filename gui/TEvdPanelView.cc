///////////////////////////////////////////////////////////////////////////////
// a view showing a single tracker panel, see murat/gui/gui/TEvdPanelView.hh
///////////////////////////////////////////////////////////////////////////////

#include "murat/gui/TEvdPanelView.hh"

namespace murat {

//-----------------------------------------------------------------------------
TEvdPanelView::TEvdPanelView(const char* Name, TGLViewer::ECameraType Type) :
  TEvdView(Name,Type),
  fPlane(-1),
  fPanel(-1)
{
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
  return 0;
}

//-----------------------------------------------------------------------------
// drawing of the panel is not implemented yet
//-----------------------------------------------------------------------------
void TEvdPanelView::Draw(Option_t* Opt) {
}

}
