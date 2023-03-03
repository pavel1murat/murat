///////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////
#include "TROOT.h"
#include "TMath.h"
#include "TApplication.h"
#include "TVirtualX.h"

#include "TGMenu.h"
#include "TGMsgBox.h"
#include "TGFrame.h"
#include "TGStatusBar.h"
#include "TGaxis.h"
#include "TText.h"
#include "TGDoubleSlider.h"
#include "TGButton.h"
#include "TGTextEntry.h"
#include "TGTextBuffer.h"
#include "TGLabel.h"

// #include "Stntuple/gui/TEvdMainFrame.hh"

#include "Stntuple/gui/TStnFrame.hh"

#include "murat/gui/TEvdManager.hh"

// ClassImp(TEvdManager)

//-----------------------------------------------------------------------------
TEvdManager::TEvdManager(const char* Name, const char* Title): TVisManager(Name, Title) {
  if (gROOT->IsBatch()) return;

  InitGui(Title);
//-----------------------------------------------------------------------------
// views
//-----------------------------------------------------------------------------
//  InitViews();
  fListOfDetectors = new TObjArray(10);
  fMinStation =  0;
  fMaxStation = 50;
					// by default, no timing constraints
  fTMin       = 0;
  fTMax       = 1.e5;

  //  fTitleNode
}

//_____________________________________________________________________________
TEvdManager::~TEvdManager() {

  if (!gROOT->IsBatch()) {

    // delete tracking views
    //    delete fTrkXYView;
    // delete fTrkTZView;

    delete fMenuBarHelpLayout;
    delete fMenuBarItemLayout;
    delete fMenu;

    delete fMenuBarLayout;
    delete fMenuBar;

    delete fMain;

    delete fListOfDetectors;
  }
}

//_____________________________________________________________________________
TEvdManager* TEvdManager::Instance() {
  if (TVisManager::fgInstance != NULL) {
    return (TEvdManager*) TVisManager::fgInstance;
  }
  else {
    return new TEvdManager("EvdManagerInstance","EvdManagerInstance");
  }
}

//_____________________________________________________________________________
void TEvdManager::HandleButtons() {
  // Handle different buttons.
  
  TGButton *btn = (TGButton *) gTQSender;
  int id = btn->WidgetId();
  
  switch (id) {
  // case kXYView:
  //   OpenTrkXYView();
  //   break;
  case TEvdManager::kTZ:
    OpenTrkTZView();
    break;
  // case TEvdManager::kCal:
  //   OpenCalView();
  //   break;
  // case TEvdManager::kCrv:
  //   OpenCrvView();
  //   break;
  case UPDATER_BTN:
    // for (unsigned int i = 0; i < 6; i++) {
    //   fCrvView[i]->SetTimeWindow(timeWindowSlider->GetMinPosition(), timeWindowSlider->GetMaxPosition());
    // }
    UpdateViews();
    break;
  default:
    printf("Unknown button clicked\n");
    break;
  }
}

//_____________________________________________________________________________
void TEvdManager::HandleSlider() {
  // Handle slider widget

  // Int_t id;
  // TGFrame *frm = (TGFrame *) gTQSender;
  // TGDoubleSlider *sd = (TGDoubleSlider *) frm;
  // id = sd->WidgetId();
  
  // switch (id) {
  //   //case TEvdManager::TIMESLIDER_ID:
  //   // Update text boxes with max and min values
    
  // case TIMESLIDER_ID:
  //   timeWindowLowDisp->SetText(boost::lexical_cast<std::string>((int) timeWindowSlider->GetMinPosition()).c_str());
  //   gClient->NeedRedraw(timeWindowLowDisp);
    
  //   timeWindowHighDisp->SetText(boost::lexical_cast<std::string>((int) timeWindowSlider->GetMaxPosition()).c_str());
  //   gClient->NeedRedraw(timeWindowHighDisp);
  //   break;
  // default:
  //   break;
  // }
}

//_____________________________________________________________________________
void TEvdManager::HandleText() {
  // Handle text entry widgets

  // TGTextEntry *te = (TGTextEntry *) gTQSender;
  // Int_t id = te->WidgetId();

  // float textBoxNum;

  // switch (id) {
  // case TIMELOW_DISP:
  //   try{
  //     textBoxNum = boost::lexical_cast<float>(timeWindowLowDisp->GetText());
  //     if (textBoxNum < 0 || textBoxNum > timeWindowSlider->GetMaxPosition())
  // 	timeWindowLowDisp->SetText(boost::lexical_cast<std::string>((int) timeWindowSlider->GetMinPosition()).c_str());
  //     else {
  // 	timeWindowSlider->SetPosition(textBoxNum, timeWindowSlider->GetMaxPosition());
  // 	timeWindowLowDisp->SetText(boost::lexical_cast<std::string>(textBoxNum).c_str());
  //     }
  //   }
  //   catch (boost::bad_lexical_cast &){
  //     timeWindowLowDisp->SetText(boost::lexical_cast<std::string>((int) timeWindowSlider->GetMinPosition()).c_str());
  //   }		
  //   break;
  // case TIMEHIGH_DISP:
  //   try {
  //     textBoxNum = boost::lexical_cast<float>(timeWindowHighDisp->GetText());
  //     if (textBoxNum > 1695 || textBoxNum < timeWindowSlider->GetMinPosition())
  // 	timeWindowHighDisp->SetText(boost::lexical_cast<std::string>((int) timeWindowSlider->GetMaxPosition()).c_str());
  //     else {
  // 	timeWindowSlider->SetPosition(timeWindowSlider->GetMinPosition(), textBoxNum);
  // 	timeWindowHighDisp->SetText(boost::lexical_cast<std::string>(textBoxNum).c_str());
  //     }
  //   }
  //   catch (boost::bad_lexical_cast &){
  //     timeWindowHighDisp->SetText(boost::lexical_cast<std::string>((int) timeWindowSlider->GetMaxPosition()).c_str());
  //   }
  //   break;
  // default:
  //   break;
  // }
}


//-----------------------------------------------------------------------------
int TEvdManager::InitGui(const char* Title) {
//   fMain = new  TEvdMainFrame(gClient->GetRoot(),200,100,kMainFrame | kVerticalFrame);
// //-----------------------------------------------------------------------------
// //  create menu bar
// //-----------------------------------------------------------------------------
//   fMenuBarLayout     = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX, 0, 0, 1, 1);
//   fMenuBarItemLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0);
//   fMenuBarHelpLayout = new TGLayoutHints(kLHintsTop | kLHintsRight);

//   fMenu = new TGPopupMenu(gClient->GetRoot());
//   fMenu->AddEntry("&Exit", M_EXIT);
//   fMenu->Associate(fMain);

//   fMenuHelp = new TGPopupMenu(gClient->GetRoot());
//   fMenuHelp->AddEntry("&Contents", M_HELP_CONTENTS);
//   fMenuHelp->AddEntry("&Search...", M_HELP_SEARCH);
//   fMenuHelp->AddSeparator();
//   fMenuHelp->AddEntry("&About", M_HELP_ABOUT);
//   fMenuHelp->Associate(fMain);

//   fMenuBar = new TGMenuBar(fMain, 1, 1, kHorizontalFrame);
//   fMenuBar->AddPopup("&Menu", fMenu, fMenuBarItemLayout);
//   fMenuBar->AddPopup("&Help", fMenuHelp, fMenuBarHelpLayout);

//   fMain->AddFrame(fMenuBar, fMenuBarLayout);

//   trkrBtnTZ = new TGTextButton(fMain, "Tracker TZ", TEvdManager::kTZ);
//   trkrBtnTZ->Connect("Clicked()", "TEvdManager", this, "HandleButtons()");
//   trkrBtnTZ->SetTextJustify(36);
//   trkrBtnTZ->SetMargins(0, 0, 0, 0);
//   trkrBtnTZ->SetWrapLength(-1);
//   trkrBtnTZ->MoveResize(16, 26, 98, 24);

  // trkrBtnRZ = new TGTextButton(fMain, "Tracker RZ", kRZView);
  // trkrBtnRZ->Connect("Clicked()", "TEvdManager", this, "HandleButtons()");
  // trkrBtnRZ->SetTextJustify(36);
  // trkrBtnRZ->SetMargins(0, 0, 0, 0);
  // trkrBtnRZ->SetWrapLength(-1);
  // trkrBtnRZ->MoveResize(16, 58, 98, 24);

  // calBtn = new TGTextButton(fMain, "Calorimeter", kCalView);
  // calBtn->Connect("Clicked()", "TEvdManager", this, "HandleButtons()");
  // calBtn->SetTextJustify(36);
  // calBtn->SetMargins(0, 0, 0, 0);
  // calBtn->SetWrapLength(-1);
  // calBtn->MoveResize(16, 90, 98, 24);

  // crvBtn = new TGTextButton(fMain, "CRV", kCrvView);
  // crvBtn->Connect("Clicked()", "TEvdManager", this, "HandleButtons()");
  // crvBtn->SetTextJustify(36);
  // crvBtn->SetMargins(0, 0, 0, 0);
  // crvBtn->SetWrapLength(-1);
  // crvBtn->MoveResize(16, 122, 98, 24);

//   timeWindowSlider = new TGDoubleHSlider(fMain, 100, kDoubleScaleBoth, TIMESLIDER_ID);
//   timeWindowSlider->SetRange(0, 1695);
//   timeWindowSlider->SetPosition(400, 1695);
//   timeWindowSlider->MoveResize(150, 45, 200, 20);
//   timeWindowSlider->Connect("PositionChanged()", "TEvdManager", this, "HandleSlider()");
	
//   timeWindowLowDisp = new TGTextEntry(fMain, timeWindowLowBuff = new TGTextBuffer(10), TIMELOW_DISP);
//   timeWindowLowBuff->AddText(0, "400");
//   timeWindowLowDisp->MoveResize(150, 70, 40, 20);
//   timeWindowLowDisp->Connect("ReturnPressed()", "TEvdManager", this, "HandleText()");

//   timeWindowHighDisp = new TGTextEntry(fMain, timeWindowHighBuff = new TGTextBuffer(10), TIMEHIGH_DISP);
//   timeWindowHighBuff->AddText(0, "1695");
//   timeWindowHighDisp->MoveResize(310, 70, 40, 20);
//   timeWindowHighDisp->Connect("ReturnPressed()", "TEvdManager", this, "HandleText()");

//   TGLabel *sliderLabelLow = new TGLabel(fMain, "0");
//   sliderLabelLow->SetTextJustify(36);
//   sliderLabelLow->SetMargins(0, 0, 0, 0);
//   sliderLabelLow->SetWrapLength(-1);
//   sliderLabelLow->MoveResize(140, 25, 30, 20);

//   TGLabel *sliderLabelHigh = new TGLabel(fMain, "1695");
//   sliderLabelHigh->SetTextJustify(36);
//   sliderLabelHigh->SetMargins(0, 0, 0, 0);
//   sliderLabelHigh->SetWrapLength(-1);
//   sliderLabelHigh->MoveResize(330, 25, 30, 20);

//   updaterBtn = new TGTextButton(fMain, "Update", UPDATER_BTN);
//   updaterBtn->Connect("Clicked()", "TEvdManager", this, "HandleButtons()");
//   updaterBtn->SetTextJustify(36);
//   updaterBtn->SetMargins(0, 0, 0, 0);
//   updaterBtn->SetWrapLength(-1);
//   updaterBtn->MoveResize(220, 120, 60, 20);
// //-----------------------------------------------------------------------------
// // final actions
// //-----------------------------------------------------------------------------
//   fMain->MapSubwindows();
//   fMain->Resize(fMain->GetDefaultSize());
//   fMain->Resize(400, 150);

//   fMain->SetWindowName(Title);
//   fMain->MapWindow();

  return 0;
}

// //-----------------------------------------------------------------------------
// int TEvdManager::InitViews() {
//   //  fTrkXYView       = new TTrkXYView();
//   // fTrkTZView       = new TTrkTZView();
//   return 0;
// }


//-----------------------------------------------------------------------------
int TEvdManager::GetViewID(const char* View) {

  TString view_id = View;
  view_id.ToLower();

  if      (view_id == "xy" ) return TEvdManager::kXY;
  else if (view_id == "rz" ) return TEvdManager::kRZ;
  else if (view_id == "tz" ) return TEvdManager::kTZ;
  else if (view_id == "cal") return TEvdManager::kCal;
  else if (view_id == "crv") return TEvdManager::kCrv;
  else if (view_id == "vst") return TEvdManager::kVST;
  else {
    printf("TEvdManager::%s: ERROR: unknown view type : %s\n",__func__,View);
    return -1;
  }
}

//-----------------------------------------------------------------------------
void TEvdManager::OpenView(const char* View) {

  TString view_id = View;
  view_id.ToLower();

  if      (view_id == "xy" ) OpenTrkXYView();
  else if (view_id == "tz" ) OpenTrkTZView();
  else {
    printf("TEvdManager::OpenView: ERROR: unknown view type : %s\n",View);
  }
}

//-----------------------------------------------------------------------------
TCanvas* TEvdManager::NewCanvas(const char* Name, const char* Title, int SizeX, int SizeY) {
  TStnFrame* win = new TStnFrame(Name, Title, this, 0, SizeX, SizeY);
  TCanvas*c = win->GetCanvas();
  DeclareCanvas(c);
  return c;
}


//_____________________________________________________________________________
Int_t TEvdManager::OpenTrkXYView() {
  // open new XY view of the detector with the default options

  int n = fListOfCanvases->GetSize();

  char name[100], title[100];

  sprintf(name, "xy_view_%i", n);
  sprintf(title, "XY view number %i", n);

  TStnFrame* win = new TStnFrame(name, title, this, TEvdManager::kXY, 740, 760);
  TCanvas* c = win->GetCanvas();
  fListOfCanvases->Add(c);

  TString name1(name);
  name1 += "_1";
  TPad* p1 = (TPad*) c->FindObject(name1);
  p1->Range(-1000., -1000., 1000., 1000.);
  p1->cd();

  TStnView* v = FindView(TEvdManager::kXY,-1);

  v->Draw();

  TString name_title(name);
  name1 += "_title";
  TPad* title_pad = (TPad*) c->FindObject(name_title);
  title_pad->cd();
  fTitleNode->Draw();

  c->Modified();
  c->Update();
  return 0;
}

//_____________________________________________________________________________
Int_t TEvdManager::OpenTrkXYView(TStnView* mother, Axis_t x1, Axis_t y1, Axis_t x2, Axis_t y2) {
	// open new XY view of the detector with the default options

  int n = fListOfCanvases->GetSize();

  char name[100], title[100];

  sprintf(name, "xy_view_%i", n);
  sprintf(title, "XY view number %i", n);

  // try to preserve the aspect ratio
  Int_t   xsize, ysize;

  xsize = 540;
  ysize = (Int_t) (xsize*TMath::Abs((y2 - y1) / (x2 - x1)) + 20);

  TStnFrame* win = new TStnFrame(name, title, this, TEvdManager::kXY, xsize, ysize);
  TCanvas* c = win->GetCanvas();
  fListOfCanvases->Add(c);

  TString name1(name);
  name1 += "_1";
  TPad* p1 = (TPad*) c->FindObject(name1);
  p1->Range(x1, y1, x2, y2);
  p1->cd();
  mother->Draw();

  TString name_title(name);
  name1 += "_title";
  TPad* title_pad = (TPad*) c->FindObject(name_title);
  title_pad->cd();
  fTitleNode->Draw();

  c->Modified();
  c->Update();
  return 0;
}

//-----------------------------------------------------------------------------
Int_t TEvdManager::OpenTrkTZView() {
  // open new TZ view of the detector with the default options

  int n = fListOfCanvases->GetSize();

  char name[100], title[100];

  sprintf(name,  "zt_view_%i", n);
  sprintf(title, "ZT view number %i", n);

  TStnFrame* win = new TStnFrame(name, title, this, TEvdManager::kXY, 1240, 760);
  TCanvas* c = win->GetCanvas();
  fListOfCanvases->Add(c);

  TString name1(name);
  name1 += "_1";
  TPad* p1 = (TPad*) c->FindObject(name1);
  p1->Range(-1600., 0., 1600., 1800.);
  p1->cd();

  TStnView* v = FindView(TEvdManager::kTZ,-1);

  v->Draw();

  TString name_title(name);
  name1 += "_title";
  TPad* title_pad = (TPad*) c->FindObject(name_title);
  title_pad->cd();
  if (fTitleNode) fTitleNode->Draw();

  c->Modified();
  c->Update();
  return 0;
}

//-----------------------------------------------------------------------------
Int_t TEvdManager::OpenTrkTZView(TStnView* Mother, Axis_t Z1, Axis_t T1, Axis_t Z2, Axis_t T2) {
	// open new XY view of the detector with the default options

  int n = fListOfCanvases->GetSize();

  char name[100], title[100];

  sprintf(name,  "zt_view_%i", n);
  sprintf(title, "ZT view number %i", n);

  // try to preserve the aspect ratio
  Int_t   xsize, ysize;

  xsize = 540;
  ysize = (Int_t) (xsize*TMath::Abs((T2 - T1) / (Z2 - Z1)) + 20);

  TStnFrame* win = new TStnFrame(name, title, this, TEvdManager::kTZ, xsize, ysize);
  TCanvas* c = win->GetCanvas();
  fListOfCanvases->Add(c);

  TString name1(name);
  name1 += "_1";
  TPad* p1 = (TPad*) c->FindObject(name1);
  p1->Range(Z1, T1, Z2, T2);
  p1->cd();
  Mother->Draw();

  TString name_title(name);
  name1 += "_title";
  TPad* title_pad = (TPad*) c->FindObject(name_title);
  title_pad->cd();
  if (fTitleNode) fTitleNode->Draw();

  c->Modified();
  c->Update();
  return 0;
}

//-----------------------------------------------------------------------------
void TEvdManager::OpenView(TStnView* Mother, int Px1, int Py1, int Px2, int Py2) {
  int vtype = Mother->Type();

  if      (vtype == TEvdManager::kXY  ) OpenTrkXYView(Mother,Px1,Py1,Px2,Py2);
  //  else if (vtype == TEvdManager::kRZ ) OpenTrkRZView(Mother,Px1,Py1,Px2,Py2);
  else if (vtype == TEvdManager::kTZ ) OpenTrkTZView(Mother,Px1,Py1,Px2,Py2);
  // else if (vtype == TEvdManager::kCal) OpenCalView  (Mother,Px1,Py1,Px2,Py2);
  // else if (vtype == TEvdManager::rv) OpenCrvView  (Mother,Px1,Py1,Px2,Py2);
  else {
    printf("TEvdManager::OpenView: ERROR: unknown view type : %i\n",vtype);
  }
}

//-----------------------------------------------------------------------------
void TEvdManager::UpdateViews() {
  TIter it(fListOfCanvases);
  while (TCanvas* c = (TCanvas*) it.Next()) {
    TIter it1(c->GetListOfPrimitives());
    while (TObject* o = it1.Next()) {
      if (o->InheritsFrom("TPad")) {
	TPad* pad = (TPad*) o;
	MarkModified(pad);
      }
    }
    c->Modified();
    c->Update();
  }
}


//_____________________________________________________________________________
void TEvdManager::CloseWindow() {
	// Called when window is closed via the window manager.

  delete this;
}

