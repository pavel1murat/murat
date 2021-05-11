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
#include "murat/gui/TEvdHelix.hh"

namespace {
  TEveScene*          _scene;
  murat::TEvdTracker*        _tracker;
  murat::TEvdStrawHitHolder* _strawHitHolder(NULL);
  TEveElementList*    _trackHolder(NULL);
}


//-----------------------------------------------------------------------------
int read_tracks(const char* TracksFile) {
  //  const char* fn = "validation_640_0004_0205_tracks.txt";
  
  FILE* f = fopen(TracksFile,"r");

  if (f == NULL) {
    printf("ERROR: read_tracks can\'t open %s, BAIL OUT\n",TracksFile);
    return -1;
  }
  
  if (_trackHolder == NULL) {
    _trackHolder = new TEveElementList("Tracks"); 
    //    _trackHolder->SetLineWidth(2);
    _scene->AddElement(_trackHolder);
  }

  _trackHolder->DestroyElements();

  char   c[1000];
  float  x0, y0, z0, r, phi0, dphidz, zmin(-1600.), zmax(1600.);

  while ((c[0]=getc(f)) != EOF) {
					// check if it is a comment line
    if (c[0] != '#') {
      ungetc(c[0],f);
      // read hit data 
      fscanf(f,"%f" ,&x0 );
      fscanf(f,"%f" ,&y0 );
      fscanf(f,"%f" ,&z0 );
      fscanf(f,"%f" ,&r );
      fscanf(f,"%f" ,&phi0 );
      fscanf(f,"%f" ,&dphidz );

      if (phi0 < 0) phi0 += 2*TMath::Pi();

      printf(" %8.3f %8.3f %8.3f %8.3f %8.5f %8.5f\n",x0,y0,z0,r,phi0,dphidz);
      //-----------------------------------------------------------------------------
      // initialize the corresponding object
      // find the right panel
      //-----------------------------------------------------------------------------
      double omega  = 1./r;
      double d0     = sqrt(x0*x0+y0*y0)-r;
      double tandip = r*dphidz;

      TEvdHelix* helix = new TEvdHelix(z0,d0,phi0+TMath::Pi()/2,omega,tandip,zmin,zmax);

      _trackHolder->AddElement(helix);

      int npt = kNStations*2*2*2; // nplanes*2*nlayers
      TEvePointSet* pset = new TEvePointSet(npt);

      TVector3  v;
      for (int ist=0; ist<kNStations; ist++) {
	for (int ipln=0; ipln<2; ipln++) {
	  for (int ip=0; ip<2; ip++) {
	    murat::TEvdPanel* panel = _tracker->Panel(ist,ipln,ip);
	    for (int iw=0; iw<2; iw++) {
	      murat::TEvdStraw* straw = panel->Straw(iw);
	      double z = straw->Z();
	      helix->GetPointAtZ(z,&v);
	      printf("helix point : ist: %2i %12.5f %12.5f %12.5f\n",ist, v.x(),v.y(),v.z());
	      pset->SetNextPoint(v.x(),v.y(),v.z());
	    }
	  }
	}
      }
      pset->SetMarkerColor(4);
      pset->SetMarkerSize(1.0);

      _scene->AddElement(pset);
    }
    fgets(c,1000,f);
  }

  fclose(f);
  
  
  return 0;
}

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
      murat::TEvdPanelStrawHitHolder* phh = _strawHitHolder->Panel(station,ipln,panel);

      murat::TEvdPanel* evd_panel = _tracker->Panel(station,ipln,panel);
      
      murat::TEvdStrawHit* hit = new murat::TEvdStrawHit();
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
      murat::TEvdStraw* s  = evd_panel->Straw(straw_number);
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
int read_tracker_geometry(const char* HitsFile, const char* TracksFile) {

  const char* fn = "trackerNumerology.txt";
  // const char* fn = "a.txt";

  FILE* f = fopen(fn,"r");

  if (f == NULL) {
    printf("ERROR: can\'t open %s, BAIL OUT\n",fn);
    return -1;
  }

  TEveManager::Create();

  _tracker = new murat::TEvdTracker();
  _tracker->InitGeometry(fn);
//-----------------------------------------------------------------------------
// create scene
//-----------------------------------------------------------------------------
  TColor::SetPalette(1, 0);
  gRandom = new TRandom3(0);

  _scene = gEve->SpawnNewScene("MyScene", "Tracker");
  _scene->SetHierarchical(kTRUE);

  _scene->AddElement(_tracker);
//-----------------------------------------------------------------------------
// read hits
//-----------------------------------------------------------------------------
  if (_strawHitHolder == NULL) {
    _strawHitHolder = new murat::TEvdStrawHitHolder();
    _strawHitHolder->IncDenyDestroy();              // protect against destruction
    //-----------------------------------------------------------------------------
    // define transformations such that hits could be placed into the local reference
    // frame of the panel
    //-----------------------------------------------------------------------------
    for (int is=0; is<20; is++) {
      for (int iplane=0; iplane<2; iplane++) {
  	for (int ip=0; ip<6; ip+=1) {
	  murat::TEvdPanel* panel = _tracker->Panel(is,iplane,ip);
  	  double phi = panel->fPhi;

  	  double zpanel = panel->Z();

	  murat::TEvdPanelStrawHitHolder* phh = _strawHitHolder->Panel(is,iplane,ip);
  	  phh->RefMainTrans().SetPos(0,0,zpanel);
  	  phh->RefMainTrans().RotatePF(1,2,phi);
  	}
      }
    }
  }

  _scene->AddElement(_strawHitHolder);
 
  //  read_hits(HitsFile);

  _strawHitHolder->ReadHits(HitsFile,_tracker);
  
  read_tracks(TracksFile);
  //-----------------------------------------------------------------------------
  // so far, hits are not displayed
  gEve->GetDefaultViewer()->AddScene(_scene);
  gEve->Redraw3D(kTRUE);
  
  return 0;
}
