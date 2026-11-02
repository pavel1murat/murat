///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////

#include <format>
#include <unordered_map>
#include <vector>

#include "TVector3.h"
#include "TPolyLine3D.h"
#include "TPolyMarker3D.h"

#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/KalSegment.hh"

#include "Stntuple/obj/TComboHit.hh"
#include "Stntuple/obj/TComboHitBlock.hh"
#include "Stntuple/obj/TStnTimeClusterBlock.hh"
#include "Stntuple/obj/TStrTrackBlock.hh"

#include "murat/gui/TEvdManager.hh"
#include "murat/gui/TEvdTrajectory.hh"
#include "murat/gui/TTrkVisNode.hh"
#include "murat/gui/TEvdComboHit.hh"


// ClassImp(TTrkVisNode)

namespace {
//-----------------------------------------------------------------------------
// display window in Y, mm
//-----------------------------------------------------------------------------
  const double kYMin   = -1000.;
  const double kYMax   =  5000.;
//-----------------------------------------------------------------------------
// samples per KalSegment
//-----------------------------------------------------------------------------
  const int    kNSteps =  200;
//-----------------------------------------------------------------------------
// marker style of the track hits: 20 is a filled circle
//-----------------------------------------------------------------------------
  const int    kTrackHitMarker = 20;

//-----------------------------------------------------------------------------
// positions along the fitted track, in time order. this is the only place
// where the track type (currently a Kalman fit with kinematic line segments) matters
//-----------------------------------------------------------------------------
  int trackPoints(const TStrTrack* Track, std::vector<float>& X, std::vector<float>& Y, std::vector<float>& Z) {
    float const c_mm_per_ns(2.99792e2);

    float tracker_pos[3] = {-3904.,    0.0, 22935.0 };

    // most cosmic particles have ny < 0

    float t1     = Track->T0()+(kYMax-Track->Y0())/Track->Ny()/c_mm_per_ns;
    float t2     = Track->T0()+(kYMin-Track->Y0())/Track->Ny()/c_mm_per_ns;
    float t_step = (t2-t1)/kNSteps;

    float pos[3];
    int   npoints{0};
                                        // just in case
    X.clear();
    Y.clear();
    Z.clear();

    for (int i=0; i<kNSteps+1; i++) {
      float t  = t1+i*t_step;
      // pos is returned on the local tracker reference frame
      Track->GetPosAtT(t,pos);
      if ((pos[2] < -3000.) or (pos[2] > 3000.)) continue;
      X.emplace_back(pos[0]+tracker_pos[0]);
      Y.emplace_back(pos[1]+tracker_pos[1]);
      Z.emplace_back(pos[2]+tracker_pos[2]);
      npoints++;
    }

    return npoints;
  }

//-----------------------------------------------------------------------------
// a combo hit is drawn as a straight segment of length 2*Ures along the wire
// direction (Ux,Uy,0), centered on the hit position
//-----------------------------------------------------------------------------
//   TPolyLine3D* comboHitSegment(const TComboHit* Hit) {
//     float tracker_pos[3] = {-3904.,    0.0, 22935.0 };

//     double dx = Hit->Ures()*Hit->Ux();
//     double dy = Hit->Ures()*Hit->Uy();

//     TPolyLine3D* seg = new TPolyLine3D(2);
//     seg->SetPoint(0, Hit->X()+tracker_pos[0]-dx, Hit->Y()+tracker_pos[1]-dy, Hit->Z()+tracker_pos[2]);
//     seg->SetPoint(1, Hit->X()+tracker_pos[0]+dx, Hit->Y()+tracker_pos[1]+dy, Hit->Z()+tracker_pos[2]);
//     seg->SetLineWidth(2.);
//     return seg;
//   }
}

namespace murat {

//-----------------------------------------------------------------------------
  TTrkVisNode::TTrkVisNode(const char* Name) :TEvdVisNode(Name) {

    // auto vm = TEvdManager::Instance();  // has to be initialized at this point

    fListOfTracks      = new TObjArray();
    fListOfClusterHits = new TObjArray();
    fListOfTrackHits   = new TObjArray();
  }

//-----------------------------------------------------------------------------
  TTrkVisNode::~TTrkVisNode() {
    fListOfTracks->Delete();
    delete fListOfTracks;

    fListOfClusterHits->Delete();
    delete fListOfClusterHits;

    fListOfTrackHits->Delete();
    delete fListOfTrackHits;
  }

//-----------------------------------------------------------------------------
  int TTrkVisNode::InitEvent() {
    int rc(0);
    //    float tracker_pos[3] = {-3904.,    0.0, 22935.0 };

    auto vm = TEvdManager::Instance();  // has to be initialized at this point

    // keep the data block names fixed
    fChBlock    = (TComboHitBlock*)       vm->FindDataBlock("ChBlock");
    fTcBlock    = (TStnTimeClusterBlock*) vm->FindDataBlock("TcBlock");
    fTrackBlock = (TStrTrackBlock*)       vm->FindDataBlock("TrackBlock");
    if ((fChBlock == nullptr) or (fTcBlock == nullptr) or (fTrackBlock == nullptr)) {
      return -1;
    }

    f_RunNumber    = fTrackBlock->RunNumber   ();
    f_EventNumber  = fTrackBlock->EventNumber ();
    f_SubrunNumber = fTrackBlock->SubrunNumber();

    fListOfTracks->Delete();
    fListOfClusterHits->Delete();
    fListOfTrackHits->Delete();

    // combo hits are located using their original index stored as the unique ID
    std::unordered_map<uint32_t,TComboHit*> hit_map;
    int nch = fChBlock->NHits();
    for (int i=0; i<nch; i++) {
      TComboHit* hit = fChBlock->Hit(i);
      hit_map[hit->GetUniqueID()] = hit;
    }

//-----------------------------------------------------------------------------
// time cluster TEvdComboHits
//-----------------------------------------------------------------------------
    TStnLinkBlock* links = fTcBlock->ListOfChLinks();
    int ntc = fTcBlock->NTimeClusters();
    for (int itc=0; itc<ntc; itc++) {
      int nlinks = links->NLinks(itc);
      for (int j=0; j<nlinks; j++) {
        auto it = hit_map.find(links->Index(itc,j));
        if (it == hit_map.end()) continue;
        // create evd representation of the combo hit
        TEvdComboHit* evd_ch = new TEvdComboHit(it->second);
        fListOfClusterHits->Add(evd_ch);
      }
    }

//-----------------------------------------------------------------------------
// tracks and track hits
//-----------------------------------------------------------------------------
    int ntrk = fTrackBlock->NTracks();
    for (int itrk=0; itrk<ntrk; itrk++) {
      TStrTrack* trk = fTrackBlock->Track(itrk);

      std::vector<float> x,y,z;

      int npoints = trackPoints(trk, x,y,z);

      std::string name = std::format("track_{}", itrk);
      TEvdTrajectory* traj = new TEvdTrajectory(name.c_str(), trk, npoints,x.data(),y.data(),z.data());
      fListOfTracks->Add(traj);

      // track hits: a dot at the position of each combo hit of the track
      std::vector<TVector3> hit_points;
      for (const auto& seed : *trk->fHitSeeds) {
        auto it = hit_map.find(seed.index());
        if (it == hit_map.end()) continue;
        TComboHit* ch = it->second;
        TEvdComboHit* evd_ch = new TEvdComboHit(ch,1);
        fListOfTrackHits->Add(evd_ch);
      }
    }

    return rc;
  }

//-----------------------------------------------------------------------------
  bool TTrkVisNode::Initialized() {
    bool initialized(false);

    auto vm = TEvdManager::Instance();
    auto tb = (TStrTrackBlock*) vm->FindDataBlock("TrackBlock");

    if (tb == nullptr) {
      return initialized;
    }

    if ((tb->EventNumber () == f_EventNumber) and
        (tb->SubrunNumber() == f_SubrunNumber) and
        (tb->RunNumber   () == f_RunNumber   )     ) {
      initialized = true;
    }
    return initialized;
  }

//-----------------------------------------------------------------------------
  void TTrkVisNode::Draw(Option_t* Opt) {
    TObjArray* lists[] = { fListOfClusterHits, fListOfTracks, fListOfTrackHits };

    for (TObjArray* list : lists) {
      int n = list->GetEntriesFast();
      for (int i=0; i<n; i++) {
        list->At(i)->Draw(Opt);
      }
    }
  }
}
