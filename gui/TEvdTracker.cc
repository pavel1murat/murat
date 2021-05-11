///////////////////////////////////////////////////////////////////////////////

#include "TGeoTube.h"
#include "murat/gui/TEvdTracker.hh"

ClassImp(murat::TEvdTracker)

namespace murat {
//-----------------------------------------------------------------------------
TEvdTracker::TEvdTracker(): TEveElementList ("tracker") {
  for (int i=0; i<kNStations; i++) {
    fStation[i] = new TEvdStation(i);
    AddElement(fStation[i]);
  }
  
  // TGeoTube* tube = new TGeoTube(300,700,1550);
  // SetShape(tube);
  
  SetRnrSelfChildren(false,true);
}


//-----------------------------------------------------------------------------
// sample file: "trackerNumerology.txt" on murat06
// need TEveManager initialized before this call
//-----------------------------------------------------------------------------
int TEvdTracker::InitGeometry(const char* FileName) {
  FILE* f = fopen(FileName,"r");
  
  if (f == NULL) {
    printf("ERROR: TTracker::InitGeometry can\'t open %s, BAIL OUT\n",FileName);
    return -1;
  }

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
      TEvdStation* s       = fStation[station];
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

  // position and rotate panels
  //-----------------------------------------------------------------------------
  for (int is=0; is<kNStations; is++) {
    TEvdStation* station = fStation[is];
    for (int ipln=0; ipln<2; ipln++) {
      TEvdPlane* pln = station->Plane(ipln);
      for (int ip=0; ip<kNPanels; ip++) {
	TEvdPanel* p = pln->Panel(ip);
	p->InitGeometry();
      }
    }
  }
  return 0;
}
}
