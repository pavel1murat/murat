#ifndef __murat_gui_TEvdCaloNode_hh__
#define __murat_gui_TEvdCaloNode_hh__

#include "TObject.h"
#include "TString.h"
#include "TGeoVolume.h"

#include "Stntuple/obj/TCaloHitBlock.hh"
#include "murat/gui/TEvdVisNode.hh"

namespace murat {
  
class TEvdCaloNode: public TEvdVisNode {
  
protected:
  TCaloHitBlock*      fCaloHitBlock;    // ==NOT OWNED==
  TGeoVolumeAssembly* fAssembly;
  
public:
					// ****** constructors and destructor
  TEvdCaloNode(const char* name = "CaloNode");
  virtual ~TEvdCaloNode();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
					// generic print callback. to be overloaded, if needed

  // virtual void        NodePrint(const void* Object, const char* ClassName) ;

  void SetCaloHitBlock(TCaloHitBlock* Block) { fCaloHitBlock = Block; }

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
