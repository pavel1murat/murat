#ifndef __murat_gui_TEvdCrvNode_hh__
#define __murat_gui_TEvdCrvNode_hh__

#include "TObject.h"
#include "TString.h"
#include "TGeoVolume.h"

#include "Stntuple/obj/TCrvClusterBlock.hh"
#include "Stntuple/obj/TCrvDigiBlock.hh"
#include "Stntuple/obj/TCrvPulseBlock.hh"

#include "murat/gui/TEvdVisNode.hh"

namespace murat {
  
class TEvdCrvNode: public TEvdVisNode {
  
protected:
  TCrvDigiBlock*      fCrvdBlock;    // ==NOT OWNED==
  TCrvPulseBlock*     fCrvpBlock;    // ==NOT OWNED==
  TCrvClusterBlock*   fCrvcBlock;    // ==NOT OWNED==
  
  TGeoVolumeAssembly* fAssembly;
  
public:
					// ****** constructors and destructor
  TEvdCrvNode(const char* name = "CrvNode");
  virtual ~TEvdCrvNode();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
					// generic print callback. to be overloaded, if needed

  // virtual void        NodePrint(const void* Object, const char* ClassName) ;

  void SetCrvcBlock(TCrvClusterBlock* Block) { fCrvcBlock = Block; }
  void SetCrvdBlock(TCrvDigiBlock*    Block) { fCrvdBlock = Block; }
  void SetCrvpBlock(TCrvPulseBlock*   Block) { fCrvpBlock = Block; }

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
