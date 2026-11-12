///////////////////////////////////////////////////////////////////////////////
// a view showing a single tracker panel, selected by a combo hit (TComboHitEvd::ShowPanel)
// the drawing is not implemented yet, only the interface
///////////////////////////////////////////////////////////////////////////////
#ifndef __murat_gui_TEvdPanelView_hh__
#define __murat_gui_TEvdPanelView_hh__

#include "murat/gui/TEvdView.hh"
#include "murat/gui/TPanelVisNode.hh"

namespace murat {

class TEvdPanelView: public TEvdView {
public:
  enum {
    kNPlanes           = 36,
    kNPanelsPerPlane   = 6,
  };

protected:
  int                    fPlane;        // plane of the panel, -1 if no panel is selected
  int                    fPanel;        // panel within the plane (0-5), -1 if no panel is selected
  TPanelVisNode*         fPanelNode;    // ==NOT OWNED== the node drawing the panel and its hits

public:
  TEvdPanelView(const char* Name, TGLViewer::ECameraType Type);
  virtual ~TEvdPanelView();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  int           GetPlane () const { return fPlane; }
  int           GetPanel () const { return fPanel; }
  TPanelVisNode* GetPanelNode() { return fPanelNode; }
//-----------------------------------------------------------------------------
// select the panel to show, returns 0 if OK, -1 if the plane or the panel is out of range
//-----------------------------------------------------------------------------
  int           SetPanel (int Plane, int Panel);

//-----------------------------------------------------------------------------
// overloaded functions of TEvdView
//-----------------------------------------------------------------------------

};
}
#endif
