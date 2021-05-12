#ifndef TEvdManager_hh
#define TEvdManager_hh

#include "TObjArray.h"
#include "Stntuple/base/TVisManager.hh"

#include "TGDoubleSlider.h"
#include "TGButton.h"
#include "TGTextEntry.h"
#include "TGTextBuffer.h"

class TControlBar;
class TGMenuBar;
class TGPopupMenu;
class TGLayoutHints;
class TGMainFrame;

class TTrkTZView;
class TSubdetector;
class TExtrapolator;

class TEvdManager : public TVisManager {
public:
//-----------------------------------------------------------------------------
// command codes
//-----------------------------------------------------------------------------
  enum CommandIdentifiers {
    M_TRACKER_XY,
    M_TRACKER_RZ,
    M_TRACKER_ZT,
    M_CALORIMETER_XY,
    M_CRV_XY,
    M_EXIT,

    M_OPTION_EVENT_STATUS,

    M_HELP_CONTENTS,
    M_HELP_SEARCH,
    M_HELP_ABOUT
  };

  enum WidgetIdentities{
    TIMESLIDER_ID = 10,
    TIMELOW_DISP  = 11,
    TIMEHIGH_DISP = 12,
    UPDATER_BTN   = 13
  };

//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
protected:
  TGMainFrame*         fMain;
  TGMenuBar           *fMenuBar;	    // !
  TGPopupMenu         *fMenu;               // !
  TGPopupMenu         *fMenuHelp;	    // !

  TGLayoutHints       *fMenuBarLayout;	    // !
  TGLayoutHints       *fMenuBarItemLayout;  // !
  TGLayoutHints       *fMenuBarHelpLayout;  // !

  TGTextButton        *trkrBtnXY, *trkrBtnTZ;
  TGTextButton*        updaterBtn;
  TGDoubleHSlider     *timeWindowSlider;
  TGTextBuffer        *timeWindowLowBuff, *timeWindowHighBuff;
  TGTextEntry         *timeWindowLowDisp, *timeWindowHighDisp;
//-----------------------------------------------------------------------------
// vis. manager also holds a list of objects to be displayed.
// The list has to be the same for all the views
//-----------------------------------------------------------------------------
  TObjArray*          fListOfDetectors;
  TSubdetector*       fClosestSubdetector;

  TTrkTZView*         fTrkTZView;

  TExtrapolator*      fExtrapolator;

  //  const art::Event*   fEvent;

  int                 fMinStation;
  int                 fMaxStation;
  int                 fTimeCluster;
  int                 fDebugLevel;

  float               fTMin;
  float               fTMax;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:

  TEvdManager(const char* name = "EvdManager",	const char* title = "EvdManager");

  virtual ~TEvdManager();

  static  TEvdManager* Instance();

  //Interface Handlers

  virtual int InitGui(const char* Title);
  //  virtual int InitViews();

  void HandleButtons();
  void HandleSlider();
  void HandleText(); //char * text);
  
  TSubdetector*  GetClosestSubdetector() { return fClosestSubdetector; }
  TExtrapolator* GetExtrapolator() { return fExtrapolator; }
  
  TObjArray*     GetListOfDetectors() { return fListOfDetectors; }

  void          AddDetector(TObject* det) { fListOfDetectors->Add(det); }

  //  const art::Event* Event() { return fEvent; }
  
  int    MinStation() { return fMinStation; }
  int    MaxStation() { return fMaxStation; }
  int    TimeCluster() { return fTimeCluster; }

  double TMin() { return fTMin; }
  double TMax() { return fTMax; }
  
  void   GetTimeWindow(double& TMin, double& TMax) {
    TMin = fTMin;
    TMax = fTMax;
  }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  // void SetEvent(art::Event& Evt) { fEvent = &Evt; }

  void SetClosestSubdetector(TSubdetector* det) { fClosestSubdetector = det; }
  void SetExtrapolator(TExtrapolator*  x) { fExtrapolator = x; }

  void UpdateViews();

  virtual TCanvas*  NewCanvas(const char* Name,
			      const char* Title,
			      Int_t       SizeX,
			      Int_t       SizeY);
//-----------------------------------------------------------------------------
// different views
//-----------------------------------------------------------------------------
  virtual void  OpenView(TStnView* Mother, int Px1, int Py, int Px2, int Py2);

  // Int_t   OpenTrkXYView();
  Int_t   OpenTrkXYView(TStnView* Mother, Axis_t x1, Axis_t y1, Axis_t x2, Axis_t y2);
  
  Int_t   OpenTrkTZView();
  Int_t   OpenTrkTZView(TStnView* Mother, Axis_t x1, Axis_t y1, Axis_t x2, Axis_t y2);
  
  void    CloseWindow();

  ClassDef(TEvdManager, 0)
};
#endif
