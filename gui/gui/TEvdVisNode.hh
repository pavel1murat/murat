#ifndef __murat_gui_TEvdVisNode_hh__
#define __murat_gui_TEvdVisNode_hh__

#include "TObject.h"
#include "TString.h"

namespace murat {
  
class TEvdVisNode: public TObject {
protected:
  TString    fName;
  TObject*   fClosestObject;
  int        fDist;

  int        fDebugLevel;
public:
					// ****** constructors and destructor
  TEvdVisNode(const char* name = "TVisNode");
  virtual ~TEvdVisNode();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  TObject*            GetClosestObject() { return fClosestObject; }
  
  virtual const char* GetName() const override   { return fName.Data(); }

  int                 DebugLevel()       { return fDebugLevel; }

					// generic print callback. to be overloaded, if needed

  // virtual void        NodePrint(const void* Object, const char* ClassName) ;

				// called by TEvdManager::DisplayEvent. a must to overload

  virtual int         InitEvent() = 0;

  void                SetDebugLevel(int Level) { fDebugLevel = Level; }

  void                SetClosestObject(TObject* Obj, int Dist) {
    fClosestObject = Obj;
    fDist          = Dist;
  }

  //  ClassDefOverride(TVisNode,0)
};
}
#endif
