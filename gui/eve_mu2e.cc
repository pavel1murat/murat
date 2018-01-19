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
#include "murat/gui/TEvdHelix.hh"

void       make_gui();
void       load_event(int IEvent);
void       update_html_summary();

// Configuration and global variables.

const char* esd_file_name = "http://root.cern.ch/files/alice_ESDs.root";
const char* esd_friends_file_name = 0;

const char* esd_geom_file_name =
   "http://root.cern.ch/files/alice_ESDgeometry.root";

namespace {
  
  Int_t esd_event_id       = 0; // Current event id.

  //  TEveTrackList *gTrackList = 0;
  
  TEvdTracker*        _tracker(NULL);
  TEvdStrawHitHolder* _strawHitHolder(NULL);
  TEveElementList*    _trackHolder(NULL);
 
  Mu2eMultiView* gMultiView = 0;

  extern HtmlSummary* fgHtmlSummary;
  //  TGHtml      *fgHtml        = 0;
}
/******************************************************************************/
// Initialization and steering functions
/******************************************************************************/

//______________________________________________________________________________
void run_eve_mu2e(int IEvent) {
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

   load_event(IEvent);

   gEve->Redraw3D(kTRUE); // Reset camera after the first event has been shown.
}


//-----------------------------------------------------------------------------
int eve_mu2e_read_tracks(const char* TracksFile) {
  
  if (_trackHolder == NULL) {
    _trackHolder = new TEveElementList("Tracks"); 
    
    _trackHolder->SetMainColor(6);
    //    _trackHolder->SetMarkerColor(kYellow);
    //    _trackHolder->SetMarkerStyle(4);
    //    _trackHolder->SetMarkerSize(0.5);

    gEve->AddElement(_trackHolder);
  }

  _trackHolder->DestroyElements();

  FILE* f = fopen(TracksFile,"r");

  if (f == NULL) {
    printf("ERROR: eve_mu2e_read_tracks can\'t open %s, BAIL OUT\n",TracksFile);
    return -1;
  }
  
  char   c[1000];
  float  z0, d0, omega, phi0, tandip, zmin(-1600.), zmax(1600.);

  while ((c[0]=getc(f)) != EOF) {
					// check if it is a comment line
    if (c[0] != '#') {
      ungetc(c[0],f);
      // read hit data 
      fscanf(f,"%f" ,&z0    );
      fscanf(f,"%f" ,&d0    );
      fscanf(f,"%f" ,&phi0  );
      fscanf(f,"%f" ,&omega );
      fscanf(f,"%f" ,&tandip);

      if (phi0 < 0) phi0 += 2*TMath::Pi();

      printf(" %8.3f %8.3f %8.3f %8.3f %8.5f\n",z0,d0,phi0,omega,tandip);
//-----------------------------------------------------------------------------
// Helix parameterization assumes that
// phi0 gives the direction of the particle at a point of closest approach
// Z0   is the z-coordinate a a point of closest approach
// omega is the signed 1/R, positive R corresponds to particle going counterclockwize
// tandip = r*Dphi/Dz
// d0 is sqrt(x0^2+y0^2)-R, y-coordinate of the point of closest approach in the
//    rotated coordinate system with Y-axis pointing from (0,0) to (x0,y0)
//-----------------------------------------------------------------------------
      TEvdHelix* helix = new TEvdHelix(z0,d0,phi0,omega,tandip,zmin,zmax);

      _trackHolder->AddElement(helix);

      int npt = kNStations*2*2*2; // nplanes*2*nlayers
      TEvePointSet* pset = new TEvePointSet(npt);

      TVector3  v;
      for (int ist=0; ist<kNStations; ist++) {
	for (int ipln=0; ipln<2; ipln++) {
	  for (int ip=0; ip<2; ip++) {
	    TEvdPanel* panel = _tracker->Panel(ist,ipln,ip);
	    for (int iw=0; iw<2; iw++) {
	      TEvdStraw* straw = panel->Straw(iw);
	      double z = straw->Z();
	      helix->GetPointAtZ(z,&v);
	      printf("helix point : ist: %2i %12.5f %12.5f %12.5f\n",ist, v.x(),v.y(),v.z());
	      pset->SetNextPoint(v.x(),v.y(),v.z());
	    }
	  }
	}
      }
      pset->SetMarkerColor(4);
      pset->SetMarkerSize(0.8);

      gEve->AddElement(pset);
    }
    fgets(c,1000,f);
  }

  fclose(f);
  
  
  return 0;
}
//-----------------------------------------------------------------------------
// this is the place where an event is getting read in
//-----------------------------------------------------------------------------
void load_event(int IEvent) {
   // Load event specified in global esd_event_id.
   // The contents of previous event are removed.

   printf("Loading event %d.\n", esd_event_id);

   //   _strawHitHolder->ReadHits("validation_640_0004_0205_hits.txt",_tracker);
   char HitsFile[200];
   sprintf(HitsFile,"validation_640_0004_%04i_hits.txt",IEvent);
   
   _strawHitHolder->ReadHits(HitsFile,_tracker);

   char TracksFile[200];
   sprintf(TracksFile,"validation_640_0004_%04i_tracks.txt",IEvent);

   eve_mu2e_read_tracks(TracksFile);
   
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
