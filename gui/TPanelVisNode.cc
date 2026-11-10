//-----------------------------------------------------------------------------
// vis node of the single-panel view
//-----------------------------------------------------------------------------
#include <format>
#include <string>

#include "TGeoMatrix.h"
#include "TGeoVolume.h"

#include "Stntuple/obj/TComboHit.hh"
#include "Stntuple/obj/TComboHitBlock.hh"

#include "murat/gui/TEvdManager.hh"
#include "murat/gui/TEvdTracker.hh"
#include "murat/gui/TComboHitEvd.hh"
#include "murat/gui/TPanelVisNode.hh"

namespace murat {

//-----------------------------------------------------------------------------
TPanelVisNode::TPanelVisNode(const char* Name) : TEvdVisNode(Name) {
  // build the geometry - that can be done in standalone mode
  fAssembly = new TGeoVolumeAssembly("panel");

  fListOfHits = new TObjArray();
  fListOfHits->SetOwner(kTRUE);
}

//-----------------------------------------------------------------------------
TPanelVisNode::~TPanelVisNode() {
  fListOfHits->Delete();
  delete fListOfHits;
}

//-----------------------------------------------------------------------------
// all combo hits of the panel: the same plane and panel as the selected hit
//-----------------------------------------------------------------------------
int TPanelVisNode::InitEvent() {
  auto vm = TEvdManager::Instance();
  fChBlock = (TComboHitBlock*) vm->FindDataBlock("ChBlock");
  if (fChBlock == nullptr) return -1;

  f_RunNumber    = fChBlock->RunNumber   ();
  f_EventNumber  = fChBlock->EventNumber ();
  f_SubrunNumber = fChBlock->SubrunNumber();

  fListOfHits->Delete();
  int nch = fChBlock->NHits();
  for (int i=0; i<nch; i++) {
    TComboHit* ch = fChBlock->Hit(i);
    if (ch->Sid().getPlane() != fPlane) continue;
    if (ch->Sid().getPanel() != fPanel) continue;
    fListOfHits->Add(new TComboHitEvd(ch));
  }
  return 0;
}

//-----------------------------------------------------------------------------
bool TPanelVisNode::Initialized() {
  auto vm = TEvdManager::Instance();
  auto cb = (TComboHitBlock*) vm->FindDataBlock("ChBlock");
  if (cb == nullptr) return false;

  return ((cb->EventNumber () == f_EventNumber ) and
          (cb->SubrunNumber() == f_SubrunNumber) and
          (cb->RunNumber   () == f_RunNumber   ));
}

//-----------------------------------------------------------------------------
// the panel volume, placed with the transform it has in the tracker, and the hits
//-----------------------------------------------------------------------------
void TPanelVisNode::Draw(Option_t* Opt) {
  auto vm = TEvdManager::Instance();
  TEvdTracker* tracker = dynamic_cast<TEvdTracker*>(vm->FindSubdetector("TRACKER"));

  if (tracker && (fPlane >= 0) && (fPanel >= 0)) {
    std::string panelName = std::format("TrackerPlane_{}_Panel_{}", fPlane, fPanel);

    TGeoNode* panelNode = nullptr;
    TObjArray* nodes = tracker->GetNodes();
    for (int i=0; i<nodes->GetEntriesFast(); i++) {
      TGeoNode* n = (TGeoNode*) nodes->UncheckedAt(i);
      if (panelName == n->GetVolume()->GetName()) { panelNode = n; break; }
    }

    if (panelNode) {
      // the view volume: the panel with the tracker transform, so that the hits line up
      static int nviews = 0;
      std::string viewName = std::format("panel_view_{}", nviews++);
      TGeoVolumeAssembly* view = new TGeoVolumeAssembly(viewName.c_str());
      view->AddNode(panelNode->GetVolume(), 1, new TGeoHMatrix(*panelNode->GetMatrix()));
      view->Draw(Opt);
    }
  }

  for (int i=0; i<fListOfHits->GetEntriesFast(); i++) {
    fListOfHits->UncheckedAt(i)->Draw(Opt);
  }
}

}
