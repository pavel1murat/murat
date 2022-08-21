///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef murat_ana_TCalAnaModule_hh
#define murat_ana_TCalAnaModule_hh

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

#include "Stntuple/loop/TStnModule.hh"

#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TStnClusterBlock.hh"
#include "Stntuple/obj/TCalDataBlock.hh"
#include "Stntuple/obj/TGenpBlock.hh"
#include "Stntuple/obj/TSimpBlock.hh"

#include "Stntuple/base/TStnArrayI.hh"

#include "Stntuple/geom/TDiskCalorimeter.hh"

#include "Stntuple/alg/TStnTrackID.hh"
#include "Stntuple/alg/TEmuLogLH.hh"
#include "Stntuple/geom/TStnCrystal.hh"

namespace murat {
class TCalAnaModule: public TStnModule {
public:
  enum { kNDisks = 2 } ;
//-----------------------------------------------------------------------------
//  histograms
//-----------------------------------------------------------------------------
  struct CrystalHist_t {
    TH1F*    fDiskID;		       // per crystal hit
    TH1F*    fEnergy  [kNDisks];
    TH1F*    fTime    [kNDisks];
    TH1F*    fNHits   [kNDisks];
    TH1F*    fRadius  [kNDisks];
    TH1F*    fRadiusWE[kNDisks];
    TH1F*    fE700    [kNDisks];
    TH1F*    fT700    [kNDisks];
    TH1F*    fN700    [kNDisks];
    TH1F*    fR700    [kNDisks];
    TH1F*    fRWE700  [kNDisks];
  };

  struct CalHitHist_t {
    TH1F*    fDiskID;		       // per crystal hit
    TH1F*    fEnergy;
    TH1F*    fTime;
  };

  struct ClusterHist_t {
    TH1F*    fDiskID;
    TH1F*    fEnergy;
    TH1F*    fT0;
    TH1F*    fRow;
    TH1F*    fCol;
    TH1F*    fX;
    TH1F*    fY;
    TH1F*    fZ;
    TH1F*    fR;
    TH1F*    fNCr0;			// all clustered
    TH1F*    fNCr1;			// above 1MeV
    TH1F*    fYMean;
    TH1F*    fZMean;
    TH1F*    fSigY;
    TH1F*    fSigZ;
    TH1F*    fSigR;
    TH1F*    fFrE1;
    TH1F*    fFrE2;
    TH1F*    fSigE1;
    TH1F*    fSigE2;
  };

  struct EventHist_t {
    TH1F*    fRv;			// MC truth information
    TH1F*    fZv;
    TH1F*    fEleCosTh;
    TH1F*    fNClusters;
    TH1F*    fNCl20;
    TH1F*    fNCl50;
    TH1F*    fNCl70;
    TH1F*    fEMax;			// energy of the first reco cluster
    TH1F*    fNGenp;                    // N(particles in GENP block)
    TH1F*    fNHits       [kNDisks];     // total number of hits per disk 
    TH1F*    fETot        [kNDisks];     // total energy deposited in the disk 
    TH1F*    fECaloTot;                   // total energy
 
					// calorimeter hit histograms

    TH2F*    fEHitVsR     [kNDisks];     // hit energy vs radius
    TH1D*    fECrVsR      [kNDisks];     // total energy_per_crystal/event vs radius

    TH1D*    fNHitsVsR    [kNDisks];            //
    TH1D*    fNCrystalsVsR[kNDisks];            //

    TH1F*    fNHitsTot;
    TH1F*    fNHitCrystalsTot;
  };

  struct GenpHist_t {
    TH1F*    fPdgCode[2];		// same distribution in different scale
    TH1F*    fGenID;			// 
    TH1F*    fZ0;			// 
    TH1F*    fT0;			// 
    TH1F*    fR0;			// 
    TH1F*    fP;			// 
    TH1F*    fCosTh;			// 
  };
// 					// histograms for the simulated CE
//   struct SimpHist_t {
//     TH1F*    fPdgCode;
//     TH1F*    fMomTargetEnd;
//     TH1F*    fMomTrackerFront;
//     TH1F*    fNStrawHits;
//   };
//-----------------------------------------------------------------------------
  enum { kNEventHistSets   = 100 };
  enum { kNClusterHistSets = 100 };
  enum { kNCalHitHistSets  =  10 };
  enum { kNCrystalHistSets = 100 };
  enum { kNGenpHistSets    = 100 };
  enum { kNSimpHistSets    = 100 };

  struct Hist_t {
    TH1F*          fCrystalR[kNDisks];	          // crystal radius
    EventHist_t*   fEvent   [kNEventHistSets];
    ClusterHist_t* fCluster [kNClusterHistSets];
    CalHitHist_t*  fCalHit  [kNCalHitHistSets];
    CrystalHist_t* fCrystal [kNCrystalHistSets];
    GenpHist_t*    fGenp    [kNGenpHistSets];
//     SimpHist_t*    fSimp    [kNSimpHistSets];
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used
  TStnTrackBlock*   fTrackBlock;
  TStnClusterBlock* fClusterBlock;
  TCalDataBlock*    fCalDataBlock;
  TGenpBlock*       fGenpBlock;
//   TSimpBlock*       fSimpBlock;
					// histograms filled
  Hist_t            fHist;

  TGenParticle*     fElectron;
  TSimParticle*     fSimp;
  double            fEleE;		// electron energy

  int               fNClusters;
  int               fNCl20;
  int               fNCl50;
  int               fNCl70;
  int               fNCalHits;
  int               fNGenp;		// N(generated particles)

  TStnTrack*        fTrack;
  TStnCluster*      fCluster;

  TDiskCalorimeter* fDiskCalorimeter;

  double            fMinT0;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TCalAnaModule(const char* name="CalAna", const char* title="CalAna");
  ~TCalAnaModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Hist_t*            GetHist        () { return &fHist;        }
  TStnTrackBlock*    GetTrackBlock  () { return fTrackBlock;   }
  TStnClusterBlock*  GetClusterBlock() { return fClusterBlock; }
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// overloaded methods of TStnModule
//-----------------------------------------------------------------------------
  int     BeginJob();
  int     BeginRun();
  int     Event   (int ientry);
  int     EndJob  ();
//-----------------------------------------------------------------------------
// other methods
//-----------------------------------------------------------------------------
  void    BookCrystalHistograms (CrystalHist_t* Hist, const char* Folder);
  void    BookClusterHistograms (ClusterHist_t* Hist, const char* Folder);
  void    BookEventHistograms   (EventHist_t*   Hist, const char* Folder);
  void    BookCalHitHistograms  (CalHitHist_t*  Hist, const char* Folder);
  void    BookGenpHistograms    (GenpHist_t*    Hist, const char* Folder);
  //  void    BookSimpHistograms    (SimpHist_t*    Hist, const char* Folder);

  void    FillCrystalHistograms  (CrystalHist_t* Hist, TStnCrystal*  Crystal);
  void    FillCalHitHistograms   (CalHitHist_t*  Hist, TCalHitData*  Hit    );
  void    FillClusterHistograms  (ClusterHist_t* Hist, TStnCluster*  Cluster);
  void    FillEventHistograms    (EventHist_t*   Hist, double        EMin , double TMin);
  void    FillGenpHistograms     (GenpHist_t*    Hist, TGenParticle* Genp   );
//   void    FillSimpHistograms     (SimpHist_t*    Hist, TSimParticle* Simp   );

  void    BookHistograms();
  void    FillHistograms();


  void    Debug();
//-----------------------------------------------------------------------------
// test
//-----------------------------------------------------------------------------
  void    Test001();

  ClassDef(murat::TCalAnaModule,0)
};
}
#endif
