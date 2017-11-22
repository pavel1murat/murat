//
//#include "TEve.h"
#include "TEveViewer.h"
#include "TEveManager.h"
#include "TEveGeoShape.h"
#include "TEveTrans.h"
#include "TEveScene.h"
#include "TRandom3.h"
#include "TGeoTube.h"
#include "TGeoMatrix.h"
#include "stdio.h"

#include "murat/gui/read_tracker_geometry.hh"
#include "murat/gui/TEvdPanel.hh"
#include "murat/gui/TEvdStraw.hh"
#include "murat/gui/TEvdStrawHit.hh"
#include "murat/gui/TEvdTracker.hh"
#include "murat/gui/TEvdStrawHitHolder.hh"

TEvdTracker*        _tracker;
TEvdStrawHitHolder* _strawHitHolder(NULL);

//-----------------------------------------------------------------------------
int read_hits(const char* HitsFile) {
  //  const char* fn = "validation_640_0003_0050_hits.txt";
  
  FILE* f = fopen(HitsFile,"r");

  if (f == NULL) {
    printf("ERROR: can\'t open %s, BAIL OUT\n",HitsFile);
    return -1;
  }

  _strawHitHolder->Clear();

  char   c[1000];
  int    index, straw_index, plane, panel, layer, straw_number, pdg_id, sim_id;
  int    delta_id, rad_ok, edep_ok;
  float  time, dt, edep, p, x, y, z, wdist, wres;

  while ((c[0]=getc(f)) != EOF) {
					// check if it is a comment line
    if (c[0] != '#') {
      ungetc(c[0],f);
      // read hit data 
      fscanf(f,"%i" ,&index );
      fscanf(f,"%i" ,&straw_index    );
      fscanf(f,"%i" ,&plane );
      fscanf(f,"%i" ,&panel );
      fscanf(f,"%i" ,&layer );
      fscanf(f,"%i" ,&straw_number );
      fscanf(f,"%f" ,&time );
      fscanf(f,"%f" ,&dt );
      fscanf(f,"%f" ,&edep );
      fscanf(f,"%f" ,&wdist );
      fscanf(f,"%f" ,&wres );
      fscanf(f,"%i" ,&pdg_id );
      fscanf(f,"%i" ,&sim_id );
      fscanf(f,"%f" ,&p );
      fscanf(f,"%f" ,&x );
      fscanf(f,"%f" ,&y );
      fscanf(f,"%f" ,&z );
      fscanf(f,"%i" ,&delta_id);
      fscanf(f,"%i" ,&rad_ok);
      fscanf(f,"%i" ,&edep_ok);

      printf(" %5i %5i %3i %3i %3i %3i %10i %10i",index,straw_index,plane, panel,layer,
	     straw_number, pdg_id, sim_id);
      printf(" %8.3f %8.3f %8.3f %9.3f %8.3f %8.3f",p,x,y,z, wdist, wres);
      printf(" %3i %3i %3i\n",delta_id, rad_ok, edep_ok);
      //-----------------------------------------------------------------------------
      // initialize the corresponding object
      // find the right panel
      //-----------------------------------------------------------------------------
      int station = plane / 2;
      int ipln    = plane % 2;
      TEvdPanelStrawHitHolder* phh = _strawHitHolder->Panel(station,ipln,panel);

      TEvdPanel* evd_panel = _tracker->Panel(station,ipln,panel);
      
      TEvdStrawHit* hit = new TEvdStrawHit();
      hit->Init(index,straw_index,plane,panel,layer,straw_number,delta_id,
		time,dt,edep,wdist,wres,x,y,z);
					// thiese are not always defined
      hit->SetPdgID(pdg_id);
      hit->SetSimID(pdg_id);

      if (delta_id >= 0) {
	hit->SetMainColor(4);
	hit->SetMainTransparency(0);
	hit->fErrorBars->SetMainColor(4);
					// leave misidentified CE hits bright
	if ((pdg_id == 11) && (p > 90.)) {
	  hit->SetMainColor(kGreen+3);
	  hit->SetMainTransparency(0);
	  hit->fErrorBars->SetMainColor(kGreen+3);
	}
      }

      // 'straw_number' - number within the panel (0:95)
      TEvdStraw* s  = evd_panel->Straw(straw_number);
      double zpanel = evd_panel->Z();
      double dz     = s->Z()-zpanel;
      
      hit->RefMainTrans().SetPos(s->Rho(),wdist,dz);
      hit->RefMainTrans().RotateLF(2,3,TMath::Pi()/2);
      hit->RefMainTrans().RotatePF(1,2,TMath::Pi()/2);

      hit->fErrorBars->RefMainTrans().SetPos(s->Rho(),wdist,dz);
      hit->fErrorBars->RefMainTrans().RotateLF(2,3,TMath::Pi()/2);
      hit->fErrorBars->RefMainTrans().RotatePF(1,2,TMath::Pi()/2);
     
      phh->AddElement(hit);
      phh->AddElement(hit->fErrorBars);
      //      break;   // *DEBUG*
    }
    fgets(c,1000,f);
  }

  fclose(f);
  
  return 0;
}

//-----------------------------------------------------------------------------
int read_tracker_geometry(const char* HitsFile) {

  const char* fn = "trackerNumerology.txt";
  // const char* fn = "a.txt";

  FILE* f = fopen(fn,"r");

  if (f == NULL) {
    printf("ERROR: can\'t open %s, BAIL OUT\n",fn);
    return -1;
  }

  TEveManager::Create();

  _tracker = new TEvdTracker();
  
  // read input file
  char   c[1000];

  int    station, plane, face, panel, layer, straw, straw_id;
  float  r, x, y, z, rho, half_length, phi, nx, ny;
  
  while ((c[0]=getc(f)) != EOF) {
					// check if it is a comment line
    if (c[0] != '#') {
      ungetc(c[0],f);
      // read channel number
      fscanf(f,"%i" ,&station   );
      fscanf(f,"%i" ,&plane    );
      fscanf(f,"%i" ,&face );
      fscanf(f,"%i" ,&panel );
      fscanf(f,"%i" ,&layer );
      fscanf(f,"%i" ,&straw );
      fscanf(f,"%i" ,&straw_id );
      fscanf(f,"%f" ,&r );
      fscanf(f,"%f" ,&x );
      fscanf(f,"%f" ,&y );
      fscanf(f,"%f" ,&z );
      fscanf(f,"%f" ,&rho );
      fscanf(f,"%f" ,&half_length );
      fscanf(f,"%f" ,&phi );
      fscanf(f,"%f" ,&nx );
      fscanf(f,"%f" ,&ny );

//      printf("  %3i %6i %5i %4i %5i %5i %10i %8.3f %10.3f %10.3f %10.3f %10.3f %10.3f %8.2f %8.4f %8.4f\n",
//	     station,plane,face, panel,
//	     layer,straw, straw_index,r,x,y,z,rho, half_length,phi,nx,ny);
      //-----------------------------------------------------------------------------
      // initialize the corresponding object
      //-----------------------------------------------------------------------------
      TEvdStation* s       = _tracker->fStation[station];
      int ip               = plane % 2;
      TEvdPlane* pln       = s->fPlane[ip];
      TEvdPanel* evd_panel = pln->fPanel[panel];

      if (straw == 0) {
	evd_panel->fNx = nx;
	evd_panel->fNy = ny;
	evd_panel->fPhi = atan2(ny,nx)-TMath::Pi();
	if (evd_panel->fPhi < -TMath::Pi()) evd_panel->fPhi += 2*TMath::Pi();
      }
      
      evd_panel->InitStraw(straw,straw_id,plane,panel,layer,rho,z,nx,ny,half_length);
    }
    fgets(c,1000,f);
  }
  fclose(f);
//-----------------------------------------------------------------------------
// create scene
//-----------------------------------------------------------------------------
  TColor::SetPalette(1, 0);
  gRandom = new TRandom3(0);

  TEveScene *scene = gEve->SpawnNewScene("MyScene", "Tracker");
  scene->SetHierarchical(kTRUE);

  //-----------------------------------------------------------------------------
  // position and rotate panels
  //-----------------------------------------------------------------------------
  for (int is=0; is<kNStations; is++) {
    TEvdStation* station = _tracker->fStation[is];
    for (int ipln=0; ipln<2; ipln++) {
      TEvdPlane* pln = station->Plane(ipln);
      for (int ip=0; ip<kNPanels; ip++) {
	TEvdPanel* p = pln->Panel(ip);
	p->InitGeometry();
      }
    }
  }
  //-----------------------------------------------------------------------------
  // at this point the geometry is initialized, need to create a view
  // add stations to the scene - this is to be done just once
  //-----------------------------------------------------------------------------
  // for (int is=0; is<20; is++) {
  //   for (int iplane=0; iplane<2; iplane++) {
  //     for (int ip=0; ip<6; ip+=1) {
  // 	TEvdPanel* panel = _tracker->fStation[is]->fPlane[iplane]->fPanel[ip];
  // 	scene->AddElement(panel);
  //     }
  //   }
  // }

  scene->AddElement(_tracker);
//-----------------------------------------------------------------------------
// read hits
//-----------------------------------------------------------------------------
  if (_strawHitHolder == NULL) {
    _strawHitHolder = new TEvdStrawHitHolder();
    _strawHitHolder->IncDenyDestroy();              // protect against destruction
    scene->AddElement(_strawHitHolder);
  //-----------------------------------------------------------------------------
  // define transformations such that hits could be placed into the local reference
  // frame of the panel
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

  read_hits(HitsFile);
  //-----------------------------------------------------------------------------
  // so far, hits are not displayed
  gEve->GetDefaultViewer()->AddScene(scene);
  gEve->Redraw3D(kTRUE);
  
  return 0;
}
