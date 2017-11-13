//
//#include "TEve.h"
#include "TEveViewer.h"
#include "TEveManager.h"
#include "TEveGeoShape.h"
#include "TEveTrans.h"
#include "TEveScene.h"
#include "TRandom3.h"
#include "TGeoTube.h"
#include "stdio.h"

#include "murat/gui/read_tracker_geometry.hh"
#include "murat/gui/TEvdPanel.hh"
#include "murat/gui/TEvdStraw.hh"
//-----------------------------------------------------------------------------
class TEvdPlane {
public:
  int         fNumber;
  TEvdPanel*  fPanel[6];

  TEvdPlane(int I = -1) {
    fNumber = I;
    for (int i=0; i<kNPanels; i++) fPanel[i] = new TEvdPanel(i);
  }
  
  TEvdPanel*  Panel(int I) { return fPanel[I]; }
};

//-----------------------------------------------------------------------------
class TEvdStation {
public:
  int         fNumber;
  TEvdPlane*  fPlane[2];

  TEvdStation(int I = -1) {
    fNumber = I;
    for (int i=0; i<2; i++) fPlane[i] = new TEvdPlane(i);
  }

  TEvdPlane*  Plane(int I) { return fPlane[I]; }
};

//-----------------------------------------------------------------------------  
class TEvdTracker {
public:
  TEvdStation*    fStation[kNStations];

  TEvdTracker() {
    for (int i=0; i<kNStations; i++) fStation[i] = new TEvdStation(i);
  }
};


TEvdTracker* _tracker;

//-----------------------------------------------------------------------------
int read_hits() {
  return 0;
}

//-----------------------------------------------------------------------------
int read_tracker_geometry() {

  const char* fn = "trackerNumerology.txt";
  // const char* fn = "a.txt";

    FILE* f = fopen(fn,"r");

  if (f == NULL) {
    printf("ERROR: can\'t open %s, BAIL OUT\n",fn);
  }

  TEveManager::Create();

  _tracker = new TEvdTracker();
  
  // read input file
  char   c[1000];

  int    station, plane, face, panel, layer, straw, straw_index;
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
      fscanf(f,"%i" ,&straw_index );
      fscanf(f,"%f" ,&r );
      fscanf(f,"%f" ,&x );
      fscanf(f,"%f" ,&y );
      fscanf(f,"%f" ,&z );
      fscanf(f,"%f" ,&rho );
      fscanf(f,"%f" ,&half_length );
      fscanf(f,"%f" ,&phi );
      fscanf(f,"%f" ,&nx );
      fscanf(f,"%f" ,&ny );

      printf("  %3i %6i %5i %4i %5i %5i %10i %8.3f %10.3f %10.3f %10.3f %10.3f %10.3f %8.2f %8.4f %8.4f\n",
	     station,plane,face, panel,
	     layer,straw, straw_index,r,x,y,z,rho, half_length,phi,nx,ny);
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
      }
      
      evd_panel->InitStraw(straw,straw_index,rho,z,nx,ny,half_length);
    }
    fgets(c,1000,f);
  }
  fclose(f);
//-----------------------------------------------------------------------------
// read hits
//-----------------------------------------------------------------------------
  read_hits();
//-----------------------------------------------------------------------------
// create scene
//-----------------------------------------------------------------------------
  TColor::SetPalette(1, 0);
  gRandom = new TRandom3(0);

  TEveScene *s = gEve->SpawnNewScene("MyScene", "Tracker");
  s->SetHierarchical(kTRUE);

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
  //-----------------------------------------------------------------------------
  
  //-----------------------------------------------------------------------------
  // draw one station - insufficient memory to display all
  //-----------------------------------------------------------------------------
  
  for (int is=0; is<20; is++) {

    for (int iplane=0; iplane<2; iplane++) {
      for (int ip=0; ip<6; ip+=1) {
	TEvdPanel* panel = _tracker->fStation[is]->fPlane[iplane]->fPanel[ip];
	s->AddElement(panel);
      }
    }
  }

  gEve->GetDefaultViewer()->AddScene(s);
  gEve->Redraw3D(kTRUE);
  
  return 0;
}
