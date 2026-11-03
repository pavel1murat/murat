#ifndef __murat_gui_TTrkVisNode_hh__
#define __murat_gui_TTrkVisNode_hh__

#include "TObject.h"
#include "TString.h"
#include "TObjArray.h"

#include "murat/gui/TEvdVisNode.hh"

class TStrawDigiBlock;
class TStrTrackBlock;
class TComboHitBlock;
class TStnTimeClusterBlock;

namespace murat {

class TTrkVisNode: public TEvdVisNode {

protected:
  TStrTrackBlock*       fTrackBlock     {nullptr}; // ==NOT OWNED==
  TStnTimeClusterBlock* fTcBlock        {nullptr}; // ==NOT OWNED==
  TComboHitBlock*       fChBlock        {nullptr}; // ==NOT OWNED==
  TStrawDigiBlock*      fSdBlock        {nullptr}; // ==NOT OWNED==
  TObjArray*            fListOfTracks   {nullptr}; // ==OWNED== TEvdTrajectory, one per track
  TObjArray*            fListOfTcHits   {nullptr}; // ==OWNED== TPolyLine3D, one per time cluster hit
  TObjArray*            fListOfTrackHits{nullptr}; // ==OWNED== TPolyLine3D, one per track hit

public:
					// ****** constructors and destructor
  TTrkVisNode(const char* name = "TrkVisNode");
  virtual ~TTrkVisNode();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
					// generic print callback. to be overloaded, if needed

  // virtual void        NodePrint(const void* Object, const char* ClassName) ;

  void SetTrackBlock(TStrTrackBlock*       Block) { fTrackBlock = Block; }
  void SetTcBlock   (TStnTimeClusterBlock* Block) { fTcBlock    = Block; }
  void SetChBlock   (TComboHitBlock*       Block) { fChBlock    = Block; }

				// called by TEvdManager::DisplayEvent. a must to overload

  virtual int         InitEvent  () override;
  virtual bool        Initialized() override;
//-----------------------------------------------------------------------------
// overloaded functions of TObject
//-----------------------------------------------------------------------------
  virtual void        Draw(Option_t* Opt = "") override;

  //  ClassDefOverride(TVisNode,0)
};
}
#endif
