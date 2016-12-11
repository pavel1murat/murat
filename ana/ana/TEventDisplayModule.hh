//-----------------------------------------------------------------------------

#ifndef __murat_ana_TEventDisplayModule_hh__
#define __murat_ana_TEventDisplayModule_hh__

// ROOT includes
// ... libCore
#include <TApplication.h>
#include <TString.h>
#include <TSystem.h>
// ... libRIO
#include <TFile.h>
// ... libGui
#include <TGString.h>
#include <TGLabel.h>
#include <TGButton.h>
#include <TGButtonGroup.h>
#include <TGTextEntry.h>
#include <TGTextView.h>
#include <TGLayout.h>
#include <TGTab.h>
#include <TG3DLine.h>
// ... libGeom
#include <TGeoManager.h>
#include <TGeoTube.h>
#include <TGeoCompositeShape.h>
#include <TGeoBoolNode.h>
// ... libEG
#include <TParticle.h>
// ... libRGL
#include <TGLViewer.h>
// ... libEve
#include <TEveManager.h>
#include <TEveEventManager.h>
#include <TEveBrowser.h>
#include <TEveGeoNode.h>
#include <TEveViewer.h>
#include <TEveScene.h>
#include <TEveProjectionManager.h>
#include <TEveProjectionAxes.h>
#include <TEvePointSet.h>
#include <TEveTrack.h>
#include <TEveTrackPropagator.h>

#include "murat/ana/TEventDisplayUtils.hh"
#include "murat/ana/TMu2eEveMagField.hh"

#include "Stntuple/loop/TStnModule.hh"
#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TGenpBlock.hh"
#include "Stntuple/obj/TSimpBlock.hh"
#include "Stntuple/gui/TStnGeoManager.hh"

class TEventDisplayModule : public TStnModule {
protected:

  // Set by parameter set variables.
  
  int             fDisplayTracks;
  bool            drawHits_;
  double          hitMarkerSize_;
  double          fTrkMaxR;
  double          fTrkMaxZ;
  double          fTrkMaxStepSize;
  double          camRotateCenterH_;
  double          camRotateCenterV_;
  double          camDollyDelta_;
  
  TEventDisplayUtils*     fEvdUtils;
  TEveGeoShape*           fSimpleGeom;

  TStnGeoManager*         fGeoManager;
  TMu2eEveMagField*       fBField;
  TEveTrackPropagator*    fTrackPropagator;
  TEveTrackPropagator*    fBeamPropagator;

  TEveViewer              *fXYView;
  TEveViewer              *fRZView;
  TEveProjectionManager   *fXYMgr;
  TEveProjectionManager   *fRZMgr;
  TEveScene               *fDetXYScene;
  TEveScene               *fDetRZScene;
  TEveScene               *fEvtXYScene;
  TEveScene               *fEvtRZScene;
  
  TGTextEntry             *fTeRun,*fTeEvt;
  TGLabel                 *fTlRun,*fTlEvt;
  
  TEveTrackList           *fBeamList;
  TEveTrackList           *fTrackList;
  TEveElementList         *fHitsList;

  TStnTrackBlock*         fTrackBlock;
  TGenpBlock*             fGenpBlock;
  TSimpBlock*             fSimpBlock;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  public:
  TEventDisplayModule(const char* Name = "EventDisplay", const char* Title = "EventDisplay");
  ~TEventDisplayModule();

  virtual int  BeginJob()           ;
  virtual int  BeginRun()           ;
  virtual int  EndJob  ()           ;
  virtual int  Event   (int IEntry) ;
  
  void MakeNavPanel();
  
  //  void HideTop            (TGeoNode* node);
  // void HideNodesByName    (TGeoNode* node, const std::string& str, bool onOff);
  // void HideNodesByMaterial(TGeoNode* node, const std::string& mat, bool onOff);
  // void HideBuilding       (TGeoNode* node);
  
  //  void SetRecursiveColorTransp(TGeoVolume *Vol, Int_t Color, Int_t Transp);
  
  ClassDef(TEventDisplayModule,0)
};

#endif
