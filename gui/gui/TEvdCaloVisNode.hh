#ifndef __murat_gui_TEvdCaloVisNode_hh__
#define __murat_gui_TEvdCaloVisNode_hh__

#include "TObject.h"
#include "TString.h"
#include "TGeoVolume.h"

#include "Stntuple/obj/TCaloHitBlock.hh"
#include "murat/gui/TEvdVisNode.hh"

namespace murat {
  
class TEvdCaloVisNode: public TEvdVisNode {
  
protected:
  TCaloHitBlock*      fCaloHitBlock;    // ==NOT OWNED==
  TGeoVolumeAssembly* fAssembly;
  
public:
					// ****** constructors and destructor
  TEvdCaloVisNode(const char* name = "TEvdCaloVisNode");
  virtual ~TEvdCaloVisNode();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
					// generic print callback. to be overloaded, if needed

  // virtual void        NodePrint(const void* Object, const char* ClassName) ;

  void SetCaloHitBlock(TCaloHitBlock* Block) { fCaloHitBlock = Block; }

				// called by TEvdManager::DisplayEvent. a must to overload

  virtual int         InitEvent();

  virtual void        Draw(Option_t* Opt = "");

  //  ClassDefOverride(TVisNode,0)
};
}
#endif
