#ifndef __murat_gui_TEvdDetNode_hh__
#define __murat_gui_TEvdDetNode_hh__

#include "TObject.h"
#include "TString.h"
#include "TGeoVolume.h"

#include "Stntuple/obj/TCaloDigiBlock.hh"
#include "Stntuple/obj/TCaloHitBlock.hh"
#include "Stntuple/obj/TStnClusterBlock.hh"

#include "Stntuple/obj/TCrvDigiBlock.hh"
#include "Stntuple/obj/TCrvPulseBlock.hh"
#include "Stntuple/obj/TCrvClusterBlock.hh"

#include "Stntuple/obj/TStrTrackBlock.hh"

#include "murat/gui/TEvdVisNode.hh"

namespace murat {
  
class TEvdDetNode: public TEvdVisNode {
  
protected:
  
  TGeoVolumeAssembly* fAssembly;

  TObjArray*          fListOfTracks;
  
public:
					// ****** constructors and destructor
  TEvdDetNode(const char* name = "DetNode");
  virtual ~TEvdDetNode();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
					// generic print callback. to be overloaded, if needed

                                        // virtual void        NodePrint(const void* Object, const char* ClassName) ;

                                        // called by TEvdManager::DisplayEvent. a must to overload
  int InitCalorimeter();
  int InitCrv();
  int InitTracker();
  int InitTracks ();
  
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
