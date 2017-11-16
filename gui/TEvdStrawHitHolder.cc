//

#include "murat/gui/TEvdStrawHitHolder.hh"

ClassImp(TEvdPanelStrawHitHolder)
ClassImp(TEvdPlaneStrawHitHolder)
ClassImp(TEvdStationStrawHitHolder)
ClassImp(TEvdStrawHitHolder)

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
