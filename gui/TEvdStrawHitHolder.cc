//
#include "TEveTrans.h"

#include "murat/gui/TEvdStrawHitHolder.hh"
#include "murat/gui/TEvdStrawHit.hh"
#include "murat/gui/TEvdTracker.hh"

ClassImp(murat::TEvdPanelStrawHitHolder)
ClassImp(murat::TEvdPlaneStrawHitHolder)
ClassImp(murat::TEvdStationStrawHitHolder)
ClassImp(murat::TEvdStrawHitHolder)

namespace murat {
//-----------------------------------------------------------------------------
TEvdPanelStrawHitHolder::TEvdPanelStrawHitHolder(int Number) {
  fNumber = Number;
  SetElementNameTitle(Form("Panel_%02i",Number),Form("Pane; %02i Hits",Number));

  SetRnrSelfChildren(false,true);
}

//-----------------------------------------------------------------------------
TEvdPanelStrawHitHolder::~TEvdPanelStrawHitHolder() {
  DestroyElements();
}


//-----------------------------------------------------------------------------
void TEvdPanelStrawHitHolder::Clear(Option_t* Opt) {
  DestroyElements();
}

//-----------------------------------------------------------------------------
TEvdPlaneStrawHitHolder::TEvdPlaneStrawHitHolder(int Number) {

  fNumber = Number;
  SetElementNameTitle(Form("Plane_%02i",Number),Form("Plane %02i Hits",Number));

  for (int i=0; i<kNPanels; i++) {
    fPanelHitHolder[i] = new TEvdPanelStrawHitHolder();
    AddElement(fPanelHitHolder[i]);
  }

  SetRnrSelfChildren(false,true);
}

//-----------------------------------------------------------------------------
TEvdPlaneStrawHitHolder::~TEvdPlaneStrawHitHolder() {

  for (int i=0; i<kNPanels; i++) {
    fPanelHitHolder[i]->DestroyElements();
    delete fPanelHitHolder[i];
  }
}


//-----------------------------------------------------------------------------
void TEvdPlaneStrawHitHolder::Clear(Option_t* Opt) {

  for (int i=0; i<kNPanels; i++) {
    fPanelHitHolder[i]->Clear();
  }
}

//-----------------------------------------------------------------------------
TEvdStationStrawHitHolder::TEvdStationStrawHitHolder(int I) {

  fNumber = I;
  SetElementNameTitle(Form("StationHits_%02i",I),Form("Station Hits %02i",I));

  for (int i=0; i<2; i++) {
    fPlaneHitHolder[i] = new TEvdPlaneStrawHitHolder(i);
    AddElement(fPlaneHitHolder[i]);
  }

  SetRnrSelfChildren(false,true);
}

//-----------------------------------------------------------------------------
TEvdStationStrawHitHolder::~TEvdStationStrawHitHolder() {

  for (int i=0; i<2; i++) {
    fPlaneHitHolder[i]->DestroyElements();
    delete fPlaneHitHolder[i];
  }
}


//-----------------------------------------------------------------------------
void TEvdStationStrawHitHolder::Clear(Option_t* Opt) {

  for (int i=0; i<2; i++) {
    fPlaneHitHolder[i]->Clear();
  }
}

//-----------------------------------------------------------------------------
TEvdStrawHitHolder::TEvdStrawHitHolder() {

  for (int i=0; i<kNStations; i++) {
    fStationHitHolder[i] = new TEvdStationStrawHitHolder(i);
    AddElement(fStationHitHolder[i]);
  }
  
  SetElementNameTitle("StrawHits","Straw Hits");
  SetRnrSelfChildren(false,true);
}

//-----------------------------------------------------------------------------
TEvdStrawHitHolder::~TEvdStrawHitHolder() {

  for (int i=0; i<kNStations; i++) {
    fStationHitHolder[i]->DestroyElements();
    delete fStationHitHolder[i];
  }
}

//-----------------------------------------------------------------------------
void TEvdStrawHitHolder::Clear(Option_t* Opt) {

  for (int i=0; i<kNStations; i++) {
    fStationHitHolder[i]->Clear();
  }
}

//-----------------------------------------------------------------------------
int TEvdStrawHitHolder::ReadHits(const char* Filename, TEvdTracker* Tracker) {
  //  const char* fn = "validation_640_0003_0050_hits.txt";
  
  FILE* f = fopen(Filename,"r");

  if (f == NULL) {
    printf("ERROR: can\'t open %s, BAIL OUT\n",Filename);
    return -1;
  }

  Clear();

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
      TEvdPanelStrawHitHolder* phh = this->Panel(station,ipln,panel);

      murat::TEvdPanel* evd_panel = Tracker->Panel(station,ipln,panel);
      
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
}
