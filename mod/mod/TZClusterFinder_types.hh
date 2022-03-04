#ifndef __CalPatRec_TZClusterFinder_types_hh__
#define __CalPatRec_TZClusterFinder_types_hh__

namespace art {
  class Event;
};

namespace fhicl {
  class ParameterSet;
};

#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/RecoDataProducts/inc/StereoHit.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/TrackerGeom/inc/Straw.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"

namespace mu2e {
  class Panel;
  class SimParticle;
//-----------------------------------------------------------------------------
// delta-electron seed: structure within the station
// doesn't own anything, no need to delete any pinters
//-----------------------------------------------------------------------------
  namespace TZClusterFinderTypes {
    enum {
      kNStations      = StrawId::_nplanes/2,   // number of tracking stations
      kNFaces         = StrawId::_nfaces*2 ,   // N(faces) per station (4)
      kNPanelsPerFace = StrawId::_npanels/2    // = 3
    };
//-----------------------------------------------------------------------------
// intersection of the two hit wires
//-----------------------------------------------------------------------------
    struct Intersection_t {
      double     x;			// x-coordinate of the intersection point
      double     y;			// y-coordinate of the intersection point
      double     t1;			// distanCe from the center of the 1st wire
      double     t2;			// distance from the center of the 2nd wire
      double     wd1;                   // distance btw the 1st hit and the intersection
      double     wd2;		        // distance btw the 2nd hit and the intersection
    };

    struct TcCandidate_t;
    
    struct HitData_t {
      const ComboHit*         fHit;
      float                   fSigW;     // cached resolution along the wire

      HitData_t(const ComboHit* Hit, /*const StrawHitPosition* Pos, const Straw* aStraw,*/ float SigW) {
	fHit         = Hit; 
	fSigW        = SigW; 
      }
    };
//-----------------------------------------------------------------------------
// MC diagnostics structure
//-----------------------------------------------------------------------------
    struct McPart_t {
      const SimParticle*            fSim;        // type undefined here
      std::vector<const HitData_t*> fListOfHits;
      int                           fFirstStation;
      int                           fLastStation;
      int                           fID;	           // SimParticle::id().asInt()
      int                           fPdgID;
      int                           fNHitsCE;
      int                           fNHitsDelta;	   // number of hits associated with reconstructed delta electrons
      float                         fTime;
      float                         fHitTMin;          // min and max times of the straw hits
      float                         fHitTMax;
      float                         fStartMom;

      McPart_t(const SimParticle* Sim = NULL) { 
	fSim          = Sim; 
	fID           = -1;
	fPdgID        = 0;
	fNHitsDelta   = 0;
	fNHitsCE      = 0;
	fFirstStation = 999;
	fLastStation  = -1;
	fStartMom     = -1;
	fTime         = 1.e6; 		// at initialization, make it absurd
	fHitTMin      = 1.e6;
	fHitTMax      = -1.e6;
      }

      ~McPart_t() { fListOfHits.clear(); }

      int NHits() { return fListOfHits.size(); }

      float Momentum () { return fStartMom; }
      float Time     () { return fTime; }
      float HitDt    () { return fHitTMax-fHitTMin; }
    }; 
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
    struct PanelZ_t {
      int                              fNHits  ;    // total number of ComboHits
      std::vector<HitData_t>           fHitData;
      const Panel*                     fPanel;      // backward pointer to the tracker panel
      double                           wx;          // direction cosines of the wires, assumed to be all the same
      double                           wy;
      double                           phi;         // phi angle of the wire
      double                           z;           // 
    }; 
//-----------------------------------------------------------------------------
// cliuster of hits within the station
//-----------------------------------------------------------------------------
    class StationTcc_t {
    public:
      int                            fNumber;		// number within the station
      int                            fStation;		// station
      int                            fType;             // defines indices of the two faces used for pre-seed seach
					                // 0: used in recovery
      //      int                            fFaceProcessed[kNFaces];
      PanelZ_t*                      panelz        [kNFaces];
      std::vector<const HitData_t*>  hitlist       [kNFaces];

      StationTcc_t() {
	fNumber           = -1;
	fStation          = -1;
	fType             =  0;

	for (int face=0; face<kNFaces; face++) {
	  panelz        [face] = NULL;
	}
      }
//------------------------------------------------------------------------------
// dont need a copy constructor
//------------------------------------------------------------------------------
      ~StationTcc_t() {}
    };
//-----------------------------------------------------------------------------
// future time cluster
//-----------------------------------------------------------------------------
    struct TCCandidate_t {
    public:
      int                   fNumber;
      int                   fFirstStation;
      int                   fLastStation;
      int                   fNComboHits;
      float                 fTime;

      TCCandidate_t() {
	fNumber = -1;
      }

      int   NComboHits() { return fNComboHits; }
      float Time      () { return fTime;       }

      int   StrawHitIndex(int I) { printf("TCCandidate_t::StrawHitIndex : define me! \n"); return -1; } 
    };
//-----------------------------------------------------------------------------
// data structure passed to the histogramming routine
//-----------------------------------------------------------------------------
    struct Data_t {
      const art::Event*             event;
      const Tracker*                tracker;
      std::vector<TCCandidate_t>    tccHolder;                   // list of time cluster candidates
      PanelZ_t                      oTracker   [kNStations][kNFaces][kNPanelsPerFace];
      int                           stationUsed[kNStations];
      const ComboHitCollection*     chcol;
      int                           debugLevel;	     // printout level
    };

//-----------------------------------------------------------------------------
// finally, utility functions
//-----------------------------------------------------------------------------
    int findIntersection(const HitData_t* Hd1, const HitData_t* Hd2, Intersection_t* Result);
  }
}
#endif
