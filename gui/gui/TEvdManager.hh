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
  ~TEvdManager();
  
public:
  static TEvdManager* Instance();

  TGeoManager* GetGeoManager() { return fGeoManager; }

  TEvdSubdetector* GetSubdetector(int I) { return (TEvdSubdetector*) fListOfSubdetectors->At(I); }

  int GetNNodes       () { return fListOfNodes->GetEntriesFast(); }
  int GetNViews       () { return fListOfViews->GetEntriesFast(); }
  int GetNSubdetectors() { return fListOfSubdetectors->GetEntriesFast(); }
  
  TEvdVisNode* GetNode(int I) { return (TEvdVisNode*) fListOfNodes->At(I); }
  TEvdView*    GetView(int I) { return (TEvdView*   ) fListOfViews->At(I); }
 
  TEvdSubdetector* FindSubdetector(const char* Name);
  TEvdView*        FindView       (const char* Name);

  int   AddNode       (TEvdVisNode*     Node);
  int   AddView       (TEvdView*        View);
  void  AddSubdetector(TEvdSubdetector* Sd) ;

  void  AddDataBlock(const char* Name, TObject* Block) {
    fListOfDataBlocks->Add(new TObjString(Name),Block);
  }

  TObject* GetDataBlock(const char* Name) {
    TObjString key(Name);
    TObject* v = fListOfDataBlocks->GetValue(&key);
    return v;
  }

  int   DisplayEvent();
  int   Update();                       // draw next event..
  
  int   InitGeometry();

  int   InitEvent();
  
  //  ClassDefOverride(murat::TEvdManager,0)
};
  
}
#endif

