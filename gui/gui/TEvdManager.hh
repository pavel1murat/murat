#ifndef __draw_calo_geometry__
#define __draw_calo_geometry__

#include <string>

#include "TObject.h"
#include "TGeoVolume.h"
#include "TGeoManager.h"

#include "murat/gui/TEvdVisNode.hh"
#include "murat/gui/TEvdView.hh"
#include "murat/gui/TEvdSubdetector.hh"

namespace murat {

//-----------------------------------------------------------------------------  
class TEvdManager : public TNamed {
public:
  TGeoVolume*  fTop;

  TGeoManager* fGeoManager;
  TObjArray*   fListOfViews;                // multiple views
  TObjArray*   fListOfNodes;                // multiple nodes
  
  TObjArray*   fListOfSubdetectors;         // each     view

  std::string  fGeometryFile;
  int          fDisplayCalorimeter;
  int          fDisplayCrv;
  int          fDisplayTracker;

private:
  TEvdManager(const char* Fcl);        // configuration file name
  
public:
  static TEvdManager* Instance(const char* Fcl = "");

  TGeoManager* GetGeoManager() { return fGeoManager; }

  TEvdSubdetector* GetSubdetector(int I) { return (TEvdSubdetector*) fListOfSubdetectors->At(I); }

  int GetNNodes() { return fListOfNodes->GetEntriesFast(); }
  int GetNViews() { return fListOfViews->GetEntriesFast(); }
  
  TEvdVisNode* GetNode(int I) { return (TEvdVisNode*) fListOfNodes->At(I); }
  TEvdView*    GetView(int I) { return (TEvdView*   ) fListOfViews->At(I); }
 
  void  AddSubdetector(TEvdSubdetector* Sd) ;

  TEvdSubdetector* FindSubdetector(const char* Name);

  int   AddView(TEvdView*    View);
  int   AddNode(TEvdVisNode* Node);

  int   InitGeometry();

  int   InitEvent();
  
  //  ClassDefOverride(murat::TEvdManager,0)
};
  
}
#endif

