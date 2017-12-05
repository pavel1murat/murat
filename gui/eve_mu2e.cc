///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
// #include "alice_esd_html_summary.C"

// Forward declarations.

class AliESDEvent;
class AliESDfriend;
class AliESDtrack;
class AliExternalTrackParam;

#include "TFile.h"
#include "TTree.h"
#include "TList.h"
#include "TEveTrans.h"
#include "TEveTrack.h"
#include "TEveGeoShape.h"
#include "TEveManager.h"
#include "TEveBrowser.h"
#include "TGFrame.h"
#include "TGButton.h"

#include "murat/gui/eve_mu2e.hh"
#include "murat/gui/eve_multiview.hh"
#include "murat/gui/eve_HtmlSummary.hh"

#include "murat/gui/TEvdTracker.hh"
#include "murat/gui/TEvdStrawHitHolder.hh"

void       make_gui();
void       load_event();
void       update_html_summary();

// void       alice_esd_read();

// TEveTrack* esd_make_track(TEveTrackPropagator* trkProp, Int_t index,
//                           AliESDtrack* at,
//                           AliExternalTrackParam* tp=0);

// Bool_t     trackIsOn(AliESDtrack* t, Int_t mask);
// void       trackGetPos(AliExternalTrackParam* tp, Double_t r[3]);
// void       trackGetMomentum(AliExternalTrackParam* tp, Double_t p[3]);
// Double_t   trackGetP(AliExternalTrackParam* tp);


// Configuration and global variables.

const char* esd_file_name = "http://root.cern.ch/files/alice_ESDs.root";
// Temporarily disable reading of ESD friend.
// There seems to be no way to get it working without AliRoot.
// const char* esd_friends_file_name =
//       "http://root.cern.ch/files/alice_ESDfriends.root";
const char* esd_friends_file_name = 0;

const char* esd_geom_file_name =
   "http://root.cern.ch/files/alice_ESDgeometry.root";

// For testing
// const char* esd_file_name         = "AliESDs.root";
// const char* esd_friends_file_name = "AliESDfriends.root";

namespace {
  // TFile *esd_file          = 0;
  // TFile *esd_friends_file  = 0;

  // TTree *esd_tree          = 0;

  // AliESDEvent  *esd        = 0;
  // TList        *esd_objs   = 0;
  // AliESDfriend *esd_friend = 0;
  
  Int_t esd_event_id       = 0; // Current event id.

  //  TEveTrackList *gTrackList = 0;
  
  TEvdTracker* _tracker     = 0;
  TEvdStrawHitHolder* _strawHitHolder(NULL);
  
  Mu2eMultiView* gMultiView = 0;

  extern HtmlSummary* fgHtmlSummary;
  //  TGHtml      *fgHtml        = 0;
}
/******************************************************************************/
// Initialization and steering functions
/******************************************************************************/

//______________________________________________________________________________
void run_eve_mu2e() {
   // Main function, initializes the application.
   //
   // 1. Load the auto-generated library holding ESD classes and
   //    ESD dictionaries.
   // 2. Open ESD data-files.
   // 3. Load cartoon geometry.
   // 4. Spawn simple GUI.
   // 5. Load first event.

  TEveManager::Create();

//-----------------------------------------------------------------------------
// load geometry and let gEve know about it
//-----------------------------------------------------------------------------
  {
  //   TFile* geom = TFile::Open(esd_geom_file_name, "CACHEREAD");
  //   if (!geom) return;
  //   TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) geom->Get("Gentle");
  //   gGeomGentle = TEveGeoShape::ImportShapeExtract(gse, 0);
  //   geom->Close();
  //   delete geom;

    _tracker = new TEvdTracker();
    _tracker->InitGeometry("trackerNumerology.txt");
    
    //    gEve->AddGlobalElement(_tracker);
    gEve->AddElement(_tracker);

    if (_strawHitHolder == NULL) {
      _strawHitHolder = new TEvdStrawHitHolder();
      _strawHitHolder->IncDenyDestroy();              // protect against destruction
      //-----------------------------------------------------------------------------
      // define transformations such that hits could be placed into the local reference
      // frame of the panel
      // tracker need to be constructed at this time
      //-----------------------------------------------------------------------------
      for (int is=0; is<20; is++) {
	for (int iplane=0; iplane<2; iplane++) {
	  for (int ip=0; ip<6; ip+=1) {
	    TEvdPanel* panel = _tracker->Panel(is,iplane,ip);
	    double phi = panel->fPhi;
	    
	    double zpanel = panel->Z();
	    
	    TEvdPanelStrawHitHolder* phh = _strawHitHolder->Panel(is,iplane,ip);
	    phh->RefMainTrans().SetPos(0,0,zpanel);
	    phh->RefMainTrans().RotatePF(1,2,phi);
	  }
	}
      }
    }
    gEve->AddElement(_strawHitHolder,NULL);
  }

//-----------------------------------------------------------------------------
// Standard multi-view 
//-----------------------------------------------------------------------------
  gMultiView = new Mu2eMultiView;

  // gMultiView->ImportGeomRPhi(gGeomGentle);
  // gMultiView->ImportGeomRhoZ(gGeomGentle);


   // // HTML summary view
   // //===================

   // fgHtmlSummary = new HtmlSummary("Alice Event Display Summary Table");
   // TEveWindowSlot* slot = TEveWindow::CreateWindowInTab(gEve->GetBrowser()->GetTabRight());
   // fgHtml = new TGHtml(0, 100, 100);
   // TEveWindowFrame *wf = slot->MakeFrame(fgHtml);
   // fgHtml->MapSubwindows();
   // wf->SetElementName("Summary");


   // Final stuff
   //=============

   gEve->GetBrowser()->GetTabRight()->SetTab(1);

   make_gui();

   load_event();

   gEve->Redraw3D(kTRUE); // Reset camera after the first event has been shown.
}

//-----------------------------------------------------------------------------
// this is the place where an event is getting read in
//-----------------------------------------------------------------------------
void load_event() {
   // Load event specified in global esd_event_id.
   // The contents of previous event are removed.

   printf("Loading event %d.\n", esd_event_id);

   _strawHitHolder->ReadHits("validation_640_0004_0205_hits.txt",_tracker);

   gEve->GetViewers()->DeleteAnnotations();

   // if (gTrackList) {
   //   gTrackList->DestroyElements();
   // }

   // esd_tree->GetEntry(esd_event_id);
   // esd_tree->Show();
   //   alice_esd_read();

   TEveElement* top = gEve->GetCurrentEvent();

   gMultiView->DestroyEventRPhi();
   gMultiView->ImportEventRPhi(top);

   gMultiView->DestroyEventRhoZ();
   gMultiView->ImportEventRhoZ(top);

   update_html_summary();

   gEve->Redraw3D(kFALSE, kTRUE);
}


/******************************************************************************/
// GUI
/******************************************************************************/

//______________________________________________________________________________
//
// EvNavHandler class is needed to connect GUI signals.

class EvNavHandler
{
public:
   void Fwd()
   {
      // if (esd_event_id < esd_tree->GetEntries() - 1) {
      //    ++esd_event_id;
      //    load_event();
      // } else {
      //    printf("Already at last event.\n");
      // }
   }
   void Bck()
   {
      // if (esd_event_id > 0) {
      //    --esd_event_id;
      //    load_event();
      // } else {
      //    printf("Already at first event.\n");
      // }
   }
};

//______________________________________________________________________________
void make_gui()
{
   // Create minimal GUI for event navigation.

   TEveBrowser* browser = gEve->GetBrowser();
   browser->StartEmbedding(TRootBrowser::kLeft);

   TGMainFrame* frmMain = new TGMainFrame(gClient->GetRoot(), 1000, 600);
   frmMain->SetWindowName("XX GUI");
   frmMain->SetCleanup(kDeepCleanup);

   TGHorizontalFrame* hf = new TGHorizontalFrame(frmMain);
   {

      TString icondir( Form("%s/icons/", gSystem->Getenv("ROOTSYS")) );
      TGPictureButton* b = 0;
      EvNavHandler    *fh = new EvNavHandler;

      b = new TGPictureButton(hf, gClient->GetPicture(icondir+"GoBack.gif"));
      hf->AddFrame(b);
      b->Connect("Clicked()", "EvNavHandler", fh, "Bck()");

      b = new TGPictureButton(hf, gClient->GetPicture(icondir+"GoForward.gif"));
      hf->AddFrame(b);
      b->Connect("Clicked()", "EvNavHandler", fh, "Fwd()");
   }
   frmMain->AddFrame(hf);

   frmMain->MapSubwindows();
   frmMain->Resize();
   frmMain->MapWindow();

   browser->StopEmbedding();
   browser->SetTabTitle("Event Control", 0);
}


/******************************************************************************/
// Code for reading AliESD and creating visualization objects
/******************************************************************************/

enum ESDTrackFlags {
   kITSin=0x0001,kITSout=0x0002,kITSrefit=0x0004,kITSpid=0x0008,
   kTPCin=0x0010,kTPCout=0x0020,kTPCrefit=0x0040,kTPCpid=0x0080,
   kTRDin=0x0100,kTRDout=0x0200,kTRDrefit=0x0400,kTRDpid=0x0800,
   kTOFin=0x1000,kTOFout=0x2000,kTOFrefit=0x4000,kTOFpid=0x8000,
   kHMPIDpid=0x20000,
   kEMCALmatch=0x40000,
   kTRDbackup=0x80000,
   kTRDStop=0x20000000,
   kESDpid=0x40000000,
   kTIME=0x80000000
};

// //______________________________________________________________________________
// void alice_esd_read()
// {
//    // Read tracks and associated clusters from current event.

//    AliESDRun    *esdrun = (AliESDRun*)    esd_objs->FindObject("AliESDRun");
//    TClonesArray *tracks = (TClonesArray*) esd_objs->FindObject("Tracks");

//    // This needs further investigation. Clusters not shown.
//    // esd_friend = (AliESDfriend*) esd_objs->FindObject("AliESDfriend");
//    // printf("Friend %p, n_tracks:%d\n",
//    //        esd_friend,
//    //        esd_friend->fTracks.GetEntries());

//    //-----------------------------------------------------------------------------
//    // add track list
//    //-----------------------------------------------------------------------------
//    if (gTrackList == 0) {
//       gTrackList = new TEveTrackList("ESD Tracks");
//       gTrackList->SetMainColor(6);
//       gTrackList->SetMarkerColor(kYellow);
//       gTrackList->SetMarkerStyle(4);
//       gTrackList->SetMarkerSize(0.5);

//       gEve->AddElement(gTrackList);
//    }

//    TEveTrackPropagator* trkProp = gTrackList->GetPropagator();
//    trkProp->SetMagField( 0.1 * esdrun->fMagneticField ); // kGaus to Tesla

//    for (Int_t n=0; n<tracks->GetEntriesFast(); ++n) {
//       AliESDtrack* at = (AliESDtrack*) tracks->At(n);

//       // If ITS refit failed, take track parameters at inner TPC radius.
//       AliExternalTrackParam* tp = at;
//       if (! trackIsOn(at, kITSrefit)) {
//          tp = at->fIp;
//       }

//       TEveTrack* track = esd_make_track(trkProp, n, at, tp);
//       track->SetAttLineAttMarker(gTrackList);
//       gTrackList->AddElement(track);

//       // This needs further investigation. Clusters not shown.
//       // if (frnd)
//       // {
//       //     AliESDfriendTrack* ft = (AliESDfriendTrack*) frnd->fTracks->At(n);
//       //     printf("%d friend = %p\n", ft);
//       // }
//    }

//    gTrackList->MakeTracks();
// }

//______________________________________________________________________________
// TEveTrack* esd_make_track(TEveTrackPropagator*   trkProp,
//                           Int_t                  index,
//                           AliESDtrack*           at,
//                           AliExternalTrackParam* tp)
// {
//    // Helper function creating TEveTrack from AliESDtrack.
//    //
//    // Optionally specific track-parameters (e.g. at TPC entry point)
//    // can be specified via the tp argument.

//    Double_t      pbuf[3], vbuf[3];
//    TEveRecTrack  rt;

//    if (tp == 0) tp = at;

//    rt.fLabel  = at->fLabel;
//    rt.fIndex  = index;
//    rt.fStatus = (Int_t) at->fFlags;
//    rt.fSign   = (tp->fP[4] > 0) ? 1 : -1;

//    trackGetPos(tp, vbuf);      rt.fV.Set(vbuf);
//    trackGetMomentum(tp, pbuf); rt.fP.Set(pbuf);

//    Double_t ep = trackGetP(at);
//    Double_t mc = 0.138; // at->GetMass(); - Complicated function, requiring PID.

//    rt.fBeta = ep/TMath::Sqrt(ep*ep + mc*mc);

//    TEveTrack* track = new TEveTrack(&rt, trkProp);
//    track->SetName(Form("TEveTrack %d", rt.fIndex));
//    track->SetStdTitle();

//    return track;
// }

//______________________________________________________________________________
// Bool_t trackIsOn(AliESDtrack* t, Int_t mask)
// {
//    // Check is track-flag specified by mask are set.

//    return (t->fFlags & mask) > 0;
// }

//______________________________________________________________________________
// void trackGetPos(AliExternalTrackParam* tp, Double_t r[3])
// {
//    // Get global position of starting point of tp.

//   r[0] = tp->fX; r[1] = tp->fP[0]; r[2] = tp->fP[1];

//   Double_t cs=TMath::Cos(tp->fAlpha), sn=TMath::Sin(tp->fAlpha), x=r[0];
//   r[0] = x*cs - r[1]*sn; r[1] = x*sn + r[1]*cs;
// }

//______________________________________________________________________________
// void trackGetMomentum(AliExternalTrackParam* tp, Double_t p[3])
// {
//    // Return global momentum vector of starting point of tp.

//    p[0] = tp->fP[4]; p[1] = tp->fP[2]; p[2] = tp->fP[3];

//    Double_t pt=1./TMath::Abs(p[0]);
//    Double_t cs=TMath::Cos(tp->fAlpha), sn=TMath::Sin(tp->fAlpha);
//    Double_t r=TMath::Sqrt(1 - p[1]*p[1]);
//    p[0]=pt*(r*cs - p[1]*sn); p[1]=pt*(p[1]*cs + r*sn); p[2]=pt*p[2];
// }

//______________________________________________________________________________
// Double_t trackGetP(AliExternalTrackParam* tp)
// {
//    // Return magnitude of momentum of tp.

//    return TMath::Sqrt(1.+ tp->fP[3]*tp->fP[3])/TMath::Abs(tp->fP[4]);
// }



//______________________________________________________________________________
void update_html_summary() {
   // Update summary of current event.

   // TEveElement::List_i i;
   // TEveElement::List_i j;
   // Int_t k;
   // TEveElement *el;
   // HtmlObjTable *table;
   // TEveEventManager *mgr = gEve ? gEve->GetCurrentEvent() : 0;
   // if (mgr) {
   //    fgHtmlSummary->Clear("D");
   //    for (i=mgr->BeginChildren(); i!=mgr->EndChildren(); ++i) {
   //       el = ((TEveElement*)(*i));
   //       if (el->IsA() == TEvePointSet::Class()) {
   //          TEvePointSet *ps = (TEvePointSet *)el;
   //          TString ename  = ps->GetElementName();
   //          TString etitle = ps->GetElementTitle();
   //          if (ename.First('\'') != kNPOS)
   //             ename.Remove(ename.First('\''));
   //          etitle.Remove(0, 2);
   //          Int_t nel = atoi(etitle.Data());
   //          table = fgHtmlSummary->AddTable(ename, 0, nel);
   //       }
   //       else if (el->IsA() == TEveTrackList::Class()) {
   //          TEveTrackList *tracks = (TEveTrackList *)el;
   //          TString ename  = tracks->GetElementName();
   //          if (ename.First('\'') != kNPOS)
   //             ename.Remove(ename.First('\''));
   //          table = fgHtmlSummary->AddTable(ename.Data(), 5,
   //                   tracks->NumChildren(), kTRUE, "first");
   //          table->SetLabel(0, "Momentum");
   //          table->SetLabel(1, "P_t");
   //          table->SetLabel(2, "Phi");
   //          table->SetLabel(3, "Theta");
   //          table->SetLabel(4, "Eta");
   //          k=0;
   //          for (j=tracks->BeginChildren(); j!=tracks->EndChildren(); ++j) {
   //             Float_t p     = ((TEveTrack*)(*j))->GetMomentum().Mag();
   //             table->SetValue(0, k, p);
   //             Float_t pt    = ((TEveTrack*)(*j))->GetMomentum().Perp();
   //             table->SetValue(1, k, pt);
   //             Float_t phi   = ((TEveTrack*)(*j))->GetMomentum().Phi();
   //             table->SetValue(2, k, phi);
   //             Float_t theta = ((TEveTrack*)(*j))->GetMomentum().Theta();
   //             table->SetValue(3, k, theta);
   //             Float_t eta   = ((TEveTrack*)(*j))->GetMomentum().Eta();
   //             table->SetValue(4, k, eta);
   //             ++k;
   //          }
   //       }
   //    }
   //    fgHtmlSummary->Build();
   //    fgHtml->Clear();
   //    fgHtml->ParseText((char*)fgHtmlSummary->Html().Data());
   //    fgHtml->Layout();
   // }
}
