#ifndef __draw_calo_geometry__
#define __draw_calo_geometry__

#include <string>

#include "TObject.h"
#include "TObjString.h"
#include "TGeoVolume.h"
#include "TGeoManager.h"
#include "TMap.h"

#include "murat/gui/TEvdVisNode.hh"
#include "murat/gui/TEvdView.hh"
#include "murat/gui/TEvdSubdetector.hh"

class TStnDataBlock;

namespace murat {

//-----------------------------------------------------------------------------  
class TEvdManager : public TNamed {
public:
  TGeoVolume*  fTop;

  TGeoManager* fGeoManager;
  TObjArray*   fListOfViews;                // multiple views
  TObjArray*   fListOfNodes;                // multiple nodes
  TMap*        fListOfDataBlocks;           // data blocks are named
  
  TObjArray*   fListOfSubdetectors;         // each     view

  std::string  fGeometryFile;
  int          fDisplayCalorimeter;
  int          fDisplayCrv;
  int          fDisplayTracker;

private:
  TEvdManager();        // configuration FCL - in .rootrc (EvdManager.ConfigFcl)
  
public:
  static TEvdManager* Instance();

  TGeoManager* GetGeoManager() { return fGeoManager; }

  TEvdSubdetector* GetSubdetector(int I) { return (TEvdSubdetector*) fListOfSubdetectors->At(I); }

  int GetNNodes       () { return fListOfNodes->GetEntriesFast(); }
  int GetNViews       () { return fListOfViews->GetEntriesFast(); }
  int GetNSubdetectors() { return fListOfSubdetectors->GetEntriesFast(); }
  
  TEvdVisNode* GetNode(int I) { return (TEvdVisNode*) fListOfNodes->At(I); }
  TEvdView*    GetView(int I) { return (TEvdView*   ) fListOfViews->At(I); }
 
  void  AddSubdetector(TEvdSubdetector* Sd) ;

  TEvdSubdetector* FindSubdetector(const char* Name);

  int   AddView(TEvdView*    View);
  int   AddNode(TEvdVisNode* Node);

  void  AddDataBlock(const char* Name, TObject* Block) {
    fListOfDataBlocks->Add(new TObjString(Name),Block);
  }

  TObject* GetDataBlock(const char* Name) {
    TObjString key(Name);
    TObject* v = fListOfDataBlocks->GetValue(&key);
    return v;
  }

  int   DisplayEvent();
  
  int   InitGeometry();

  int   InitEvent();
  
  //  ClassDefOverride(murat::TEvdManager,0)
};
  
}
#endif

