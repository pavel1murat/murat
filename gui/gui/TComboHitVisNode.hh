///////////////////////////////////////////////////////////////////////////////
// TZ view prototyping
///////////////////////////////////////////////////////////////////////////////
#ifndef murat_gui_TComboHitVisNode_hh
#define murat_gui_TComboHitVisNode_hh

#include "Gtypes.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TPad.h"
#include "TArc.h"

#include "murat/gui/TComboHitData.hh"

#include "Stntuple/gui/TStnVisNode.hh"

// #ifndef __CINT__
// #include "RecoDataProducts/inc/StrawHitCollection.hh"
// #include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
// #include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
// #include "RecoDataProducts/inc/TimeCluster.hh"

// #else
// namespace mu2e {
//   class StrawHitCollection;
//   class StrawHitPositionCollection;
//   class StrawHitFlagCollection;
// };
//#endif

namespace mu2e {
  class Tracker;
}

namespace murat {
  
class TComboHitVisNode: public TStnVisNode {
public:
  enum {
    kPickHits     = 0,
    kPickTracks   = 1,
    kPickClusters = 2
  };
  
protected:
  //  TObjArray**    fListOfClusters;

  // const mu2e::StrawHitCollection**             fStrawHitColl;
  // const mu2e::StrawHitPositionCollection**     fStrawHitPosColl;  //
  // const mu2e::StrawHitFlagCollection**         fStrawHitFlagColl; //
  // const mu2e::TimeClusterCollection**          fTimeClusterColl;  //
 
  // TArc*         fArc;

  HitData_t*    fHitData; 		// his should be enough for initialization

  TObjArray*    fListOfHits;

  //  const mu2e::TimeCluster*  fTimePeak;

  Int_t         fDisplayBackgroundHits;
  Int_t         fTimeWindow;
  Int_t         fPickMode;
  double        fEventTime;

public:
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
  TComboHitVisNode();
  TComboHitVisNode(const char* Name, HitData_t* Data); 

  virtual ~TComboHitVisNode();

  int NHits() { return fListOfHits->GetEntriesFast(); }

  TComboHitData* Hit(int I) { return (TComboHitData*) fListOfHits->UncheckedAt(I); } 

  void  SetHitData(HitData_t* Data) { fHitData = Data; }
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  // const mu2e::StrawHitCollection* GetStrawHitColl() { 
  //   return *fStrawHitColl; 
  // }

  // const mu2e::StrawHitPositionCollection* GetStrawHitPosColl() { 
  //   return *fStrawHitPosColl;
  // }

  // const mu2e::StrawHitFlagCollection* GetStrawHitFlagColl() { 
  //   return *fStrawHitFlagColl;
  // }

  // const mu2e::PtrStepPointMCVectorCollection* GetMcPtrColl() { 
  //   return *fMcPtrColl;
  // }

  int DisplayBackgroundHits() { return fDisplayBackgroundHits; }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  // void SetStrawHitColl(const mu2e::StrawHitCollection** Coll) { 
  //   fStrawHitColl = Coll;
  // }

  // void SetStrawHitPosColl (const mu2e::StrawHitPositionCollection** Coll) { 
  //   fStrawHitPosColl = Coll;
  // }

  // void SetStrawHitFlagColl(const mu2e::StrawHitFlagCollection** Coll) { 
  //   fStrawHitFlagColl = Coll;
  // }

  // void SetMcPtrColl(const mu2e::PtrStepPointMCVectorCollection** Coll) { 
  //   fMcPtrColl = Coll;
  // }

  // void SetTimeClusterColl(const mu2e::TimeClusterCollection** Coll) { 
  //   fTimeClusterColl = Coll;
  // }

  void  SetPickMode   (Int_t Mode) { fPickMode    = Mode; }

  void  SetDisplayBackgroundHits(Int_t Mode) { fDisplayBackgroundHits = Mode; }

  //  virtual void  Draw    (Option_t* option = "");

  virtual int InitEvent();

  virtual void  PaintXY (Option_t* option = "") ;
  virtual void  PaintTZ (Option_t* option = "") ;

  //  virtual void  ExecuteEvent(Int_t event, Int_t px, Int_t py);

  // virtual Int_t DistancetoPrimitiveXY(Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveTZ(Int_t px, Int_t py) ;

  virtual void   Print(const char* Opt = "") const ; // **MENU**

  ClassDef(TComboHitVisNode,0)
};
}

#endif
