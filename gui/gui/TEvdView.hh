///////////////////////////////////////////////////////////////////////////////
// a view 
///////////////////////////////////////////////////////////////////////////////
#ifndef __murat_gui_TEvdView_hh__
#define __murat_gui_TEvdView_hh__

#include "TMarker.h"
#include "TNamed.h"
#include "TObjArray.h"
#include "TGaxis.h"
#include "TVector3.h"
#include "TCanvas.h"

#include "TGLViewer.h"

#include "murat/gui/TEvdVisNode.hh"

namespace murat {
  
class TEvdView: public TNamed {
  //protected:
public:
  //  int                    fType;         // view type
  int                    fIndex;        // for calorimeter - 2 views, for example
  int                    fOpen;         // 1 if the window is open
  void*                  fMother;       // non-null. if in the local ref system of some object

  static int             fgDebugLevel;
   
  TMarker*               fCenter;

  TObjArray*             fListOfNodes;	// list of TEvdVisNode's
  TCanvas*               fCanvas;       // each view has a pad corresponging to it [NOT OWNED]

  TGLViewer::ECameraType fCameraType;
  TGLViewer*             fViewer;

protected:
  float     fTMin;
  float     fTMax;

public:
  TEvdView(const char* Name, TGLViewer::ECameraType Type, int Index = -1); 

  TEvdView(const char* Name, TGLViewer::ECameraType Type, int Index, const char* Title);

  virtual ~TEvdView();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  // int           Type () { return fType;  }
  int           Index() { return fIndex; }

  void*         GetMother()      { return fMother; }
  int           GetNNodes()      { return fListOfNodes->GetEntriesFast(); }
  TEvdVisNode*  GetNode  (int I) { return (TEvdVisNode*) fListOfNodes->UncheckedAt(I);   }
  TObjArray*    GetListOfNodes() { return fListOfNodes; }
  TGLViewer*    GetViewer()      { return fViewer; }
  
  void          AddNode(TEvdVisNode* Node) { fListOfNodes->Add(Node); }

  bool          IsOpen() { return (fCanvas != nullptr); }

  int           Open();

  void          Update();
//-----------------------------------------------------------------------------
// setters
//-----------------------------------------------------------------------------
  void          SetDebugLevel (int Level);               // *MENU*
  void          SetIndex      (int Index) { fIndex = Index; } 

  void          SetMother     (void* Mother) { fMother = Mother; };
  void          SetTimeWindow (float TMin, float TMax);  // *MENU* 
  //  void          SetType       (int Type ) { fType  = Type;  } 
//-----------------------------------------------------------------------------
// overloaded functions of TObject
//-----------------------------------------------------------------------------
  virtual void  Draw(Option_t* Opt = "") override;
  // virtual void  Paint               (Option_t* option = "") override;
  // virtual void  Print               (Option_t* option = "") const override;  // *MENU* 

  //  ClassDefOverride(TEvdView,0)
};
}
#endif
