///////////////////////////////////////////////////////////////////////////////
// May 04 2013 P.Murat
// 
// in 'XY' mode draw calorimeter clusters as circles with different colors 
// in 'Cal' mode draw every detail...
///////////////////////////////////////////////////////////////////////////////
#include "TVirtualX.h"
#include "TPad.h"
#include "TStyle.h"
#include "TVector3.h"
#include "TLine.h"
#include "TArc.h"
#include "TArrow.h"
#include "TBox.h"
// #include "art/Framework/Principal/Event.h"
// #include "art/Framework/Principal/Handle.h"

// #include "GeometryService/inc/GeometryService.hh"
// #include "GeometryService/inc/GeomHandle.hh"
// #include "ConditionsService/inc/ConditionsHandle.hh"
// #include "TrackerConditions/inc/StrawResponse.hh"

// #include "Stntuple/gui/TEvdTrack.hh"
// #include "Stntuple/gui/TTrkVisNode.hh"
// #include "Stntuple/gui/TEvdStraw.hh"
// #include "Stntuple/gui/TEvdStrawHit.hh"
// #include "Stntuple/gui/TEvdTrkStrawHit.hh"
// #include "Stntuple/gui/TEvdStation.hh"
// #include "Stntuple/gui/TEvdPanel.hh"
// #include "Stntuple/gui/TEvdPlane.hh"
// #include "Stntuple/gui/TEvdStrawTracker.hh"

#include "murat/gui/TEvdManager.hh"

// #include "RecoDataProducts/inc/StrawHitCollection.hh"

// #include "DataProducts/inc/StrawId.hh"
// #include "DataProducts/inc/XYZVec.hh"

#include "murat/gui/TComboHitVisNode.hh"

ClassImp(TComboHitVisNode)

//_____________________________________________________________________________
TComboHitVisNode::TComboHitVisNode() : TStnVisNode("unnamed_TComboHitVisNode") {
}

//_____________________________________________________________________________
TComboHitVisNode::TComboHitVisNode(const char* Name, HitData_t* Data): TStnVisNode(Name) {
  // fTracker    = new stntuple::TEvdStrawTracker(Tracker);
  // fTrackBlock = TrackBlock;

  fHitData    = Data;
  fEventTime  = 0;
  fTimeWindow = 1.e6;

  fListOfHits  = new TObjArray();
  //  fTimeCluster = NULL;
}

//_____________________________________________________________________________
TComboHitVisNode::~TComboHitVisNode() {
  // delete fArc;
  
  fListOfHits->Delete();
  delete fListOfHits;
}


//-----------------------------------------------------------------------------
int TComboHitVisNode::InitEvent() {

  //  const char* oname = "TComboHitVisNode::InitEvent";

  // mu2e::GeomHandle<mu2e::Tracker> ttHandle;
  // const mu2e::Tracker* tracker = ttHandle.get();

  // Tracker calibration object.
  //  mu2e::ConditionsHandle<mu2e::StrawResponse> srep = mu2e::ConditionsHandle<mu2e::StrawResponse>("ignored");

  // const mu2e::ComboHit              *hit;
  // const mu2e::StrawDigiMC           *hit_digi_mc;

  // stntuple::TEvdStrawHit                      *evd_straw_hit; 
  // const CLHEP::Hep3Vector           /**mid,*/ *w; 
  // const mu2e::Straw                 *straw; 

  // int                               n_straw_hits, display_hit, color, nl, ns; // , ipeak, ihit;
  // bool                              isFromConversion, intime;
  // double                            sigw, /*vnorm, v,*/ sigr; 
  // CLHEP::Hep3Vector                 vx0, vx1, vx2;
//-----------------------------------------------------------------------------
// first, clear the cached hit information from the previous event
//-----------------------------------------------------------------------------
  // stntuple::TEvdStation*   station;
  // stntuple::TEvdPlane*     plane;
  // stntuple::TEvdPanel*     panel;
  //  TEvdFace*      face;

  fListOfHits->Clear();

  for (int is=0; is<18; is++) {
    int n = fHitData->fStation[is].GetNHits();
    for (int i=0; i<n; i++) {
      TComboHitData* hit = fHitData->fStation[is].GetHit(i);
      fListOfHits->Add(hit);
    }
  }

  return 0;
}

//-----------------------------------------------------------------------------
// draw reconstructed tracks and straw hits
//-----------------------------------------------------------------------------
void TComboHitVisNode::PaintXY(Option_t* Option) {
}

//-----------------------------------------------------------------------------
// draw reconstructed tracks and straw hits
//-----------------------------------------------------------------------------
void TComboHitVisNode::PaintTZ(Option_t* Option) {

  int nhits = fListOfHits->GetEntries();
  for (int i=0; i<nhits; i++) {
    TComboHitData* hit = (TComboHitData*) fListOfHits->At(i);

    // float time  = hit->Time();

    //    if ((time >= tmin) && (time <= tmax)) {
    hit->Paint(Option);
    // }
  }

  gPad->Modified();
}

//_____________________________________________________________________________
Int_t TComboHitVisNode::DistancetoPrimitiveTZ(Int_t px, Int_t py) {

  static TVector3 global;
  global.SetXYZ(gPad->AbsPixeltoX(px),gPad->AbsPixeltoY(py),0);

  TObject* closest(nullptr);

  int  x1, y1, dx1, dy1, min_dist(9999), dist;

  int nhits = fListOfHits->GetEntries();
  for (int i=0; i<nhits; i++) {
    TComboHitData* hit = (TComboHitData*) fListOfHits->At(i);
    x1  = gPad->XtoAbsPixel(hit->Z());
    y1  = gPad->YtoAbsPixel(hit->T());
    dx1 = px-x1;
    dy1 = py-y1;

    dist  = (int) sqrt(dx1*dx1+dy1*dy1);
    if (dist < min_dist) {
      min_dist = dist;
      closest  = hit;
    }
  }

  SetClosestObject(closest,min_dist);

  return min_dist;
}

//-----------------------------------------------------------------------------
void TComboHitVisNode::Print(Option_t* Option) const {
  //TComboHitVisNode* node = (TComboHitVisNode*) this;

  printf("TComboHitVisNode name: %s\n",GetName());
  printf("nhits: %5i\n",fListOfHits->GetEntriesFast());
  for (int i=0; i<18; i++) {
    printf(" %3i",fHitData->fStation[i].fHits.GetEntries());
  }
  printf("\n");
}
