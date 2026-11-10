//-----------------------------------------------------------------------------
// vis node of the single-panel view: a panel and all combo hits of the panel
//-----------------------------------------------------------------------------
#ifndef __murat_gui_TPanelVisNode_hh__
#define __murat_gui_TPanelVisNode_hh__

#include "TObjArray.h"
#include "murat/gui/TEvdVisNode.hh"

class TComboHitBlock;
class TStrTrackBlock;

namespace murat {

class TPanelVisNode: public TEvdVisNode {
protected:
  TGeoVolumeAssembly* fAssembly;
  int             fPlane{-1};
  int             fPanel{-1};
  TComboHitBlock* fChBlock{nullptr};       // ==NOT OWNED==
  TStrTrackBlock* fTrackBlock;
  TObjArray*      fListOfHits{nullptr};    // ==OWNED== TComboHitEvd objects of the panel

public:
  TPanelVisNode(const char* name = "PanelVisNode");
  virtual ~TPanelVisNode();

  void SetPanel(int Plane, int Panel) { fPlane = Plane; fPanel = Panel; }

  virtual int  InitEvent  () override;
  virtual bool Initialized() override;
  virtual void Draw(Option_t* Opt = "") override;
};
}
#endif
