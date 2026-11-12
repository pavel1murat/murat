//-----------------------------------------------------------------------------
// vis node of the single-panel view: a panel and all combo hits of the panel
//-----------------------------------------------------------------------------
#ifndef __murat_gui_TPanelVisNode_hh__
#define __murat_gui_TPanelVisNode_hh__

#include "TObjArray.h"
#include "murat/gui/TEvdVisNode.hh"

class TComboHitBlock;
class TGeoVolume;
class TStrTrackBlock;

namespace murat {

class TPanelVisNode: public TEvdVisNode {
protected:
  TGeoVolume*     fAssembly{nullptr};      // ==NOT OWNED== the panel with its straws, in the panel frame
  int             fPlane{-1};
  int             fPanel{-1};
  int             fBuiltPlane{-2};         // the plane and panel fAssembly was built for
  int             fBuiltPanel{-2};
  TComboHitBlock* fChBlock{nullptr};       // ==NOT OWNED==
  TStrTrackBlock* fTrackBlock;
  TObjArray*      fListOfHits{nullptr};    // ==OWNED== TComboHitEvd objects of the panel

public:
  TPanelVisNode(const char* name = "PanelVisNode");
  virtual ~TPanelVisNode();

  void SetPanel(int Plane, int Panel) { fPlane = Plane; fPanel = Panel; }
  double Azimuth();                        // azimuth of the panel in the tracker, [rad]

  virtual int  InitEvent  () override;
  virtual bool Initialized() override;
  virtual void Draw(Option_t* Opt = "") override;
};
}
#endif
