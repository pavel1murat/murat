///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef murat_ana_TCosmicsAnaModule_hh
#define murat_ana_TCosmicsAnaModule_hh

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TStnClusterBlock.hh"
#include "Stntuple/obj/TCalDataBlock.hh"
#include "Stntuple/obj/TStrawHitBlock.hh"
#include "Stntuple/obj/TStepPointMCBlock.hh"
#include "Stntuple/obj/TGenpBlock.hh"
#include "Stntuple/obj/TSimpBlock.hh"

#include "Stntuple/obj/TCrvClusterBlock.hh"
#include "Stntuple/obj/TCrvPulseBlock.hh"

#include "Stntuple/base/TStnArrayI.hh"

#include "Stntuple/geom/TDiskCalorimeter.hh"
#include "Stntuple/geom/TStnCrystal.hh"

#include "Stntuple/alg/TStnTrackID.hh"
#include "Stntuple/alg/TEmuLogLH.hh"

#include "murat/ana/AnaDefs.hh"
#include "murat/ana/HistBase_t.h"
#include "murat/ana/TrackPar_t.hh"
#include "murat/ana/SimPar_t.hh"

#include "murat/ana/TAnaModule.hh"

namespace murat {
class TCosmicsAnaModule: public TAnaModule {
public:
//-----------------------------------------------------------------------------
// track and sim particle additional parameters
//-----------------------------------------------------------------------------
  enum { kNDisks        =   2 } ;
  enum { kNTrackBlocks  =   8 } ; // dem:dmm:dep:dmp:uem:umm:uep:ump
  enum { kMaxTrackID    =  10 } ;
  enum { kMaxNErrors    = 100 } ;
//-----------------------------------------------------------------------------
//  histograms
//-----------------------------------------------------------------------------
  struct CaloHist_t : public HistBase_t {
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

  struct ClusterHist_t : public HistBase_t {
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

  struct TrackEffHist_t : public HistBase_t {
    TH1F*    fPtMc;			// denominator
    TH1F*    fPtReco;			// numerator
  };

  struct Error_t {
    int fNReports;
    int fMaxNReports;
  };

//-----------------------------------------------------------------------------
  enum { kNEventHistSets   =  100 };
  enum { kNTrackHistSets   = 1000 };
  enum { kNClusterHistSets =  100 };
  enum { kNCaloHistSets    =  100 };
  enum { kNGenpHistSets    =  100 };
  enum { kNSimpHistSets    =  100 };
  enum { kNTrackIDHistSets =   10 };

  struct Hist_t {
    TH1F*                fCrystalR[2];	          // crystal radius
    EventHist_t*         fEvent  [kNEventHistSets];
    TrackHist_t*         fTrack  [kNTrackHistSets];
    ClusterHist_t*       fCluster[kNClusterHistSets];
    CaloHist_t*          fCalo   [kNCaloHistSets];
    GenpHist_t*          fGenp   [kNGenpHistSets];
    SimpHist_t*          fSimp   [kNSimpHistSets];
    TStnTrackID::Hist_t* fTrackID[kNTrackIDHistSets];
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used
  TStnTrackBlock*    fTrackBlockDem;
  TStnTrackBlock*    fTrackBlockDmm;
  TStnTrackBlock*    fTrackBlockDep;
  TStnTrackBlock*    fTrackBlockDmp;
  TStnTrackBlock*    fTrackBlockUem;
  TStnTrackBlock*    fTrackBlockUmm;
  TStnTrackBlock*    fTrackBlockUep;
  TStnTrackBlock*    fTrackBlockUmp;

  TStnTrackBlock*    fTrackBlock[kNTrackBlocks]; // same thing, cached pointers ...

  TStnClusterBlock*  fClusterBlock;
  TCalDataBlock*     fCalDataBlock;
  TStrawHitBlock*    fStrawHitBlock;
  TGenpBlock*        fGenpBlock;
  TSimpBlock*        fSimpBlock;
  TStepPointMCBlock* fVDetBlock;

  TCrvClusterBlock*  fCrvClusterBlock;
  TCrvPulseBlock*    fCrvPulseBlock;
					// additional track parameters (assume ntracks < 20)

  TrackPar_t        fTrackPar[kNTrackBlocks][20];
  SimPar_t          fSimPar;
					// histograms filled
  Hist_t            fHist;
					// cut values
  double            fPtMin;

  TGenParticle*     fParticle;		// electron or muon

  //   double            fEleE;		// electron energy

  int               fCalorimeterType;

  int               fNClusters;
  int               fNTracks       [kNTrackBlocks];	// for a given track block
  int               fNGoodTracks   [kNTrackBlocks];
  int               fNMatchedTracks[kNTrackBlocks];

  int               fNTrkUNeg;
  int               fNTrkDNeg;
  int               fNTrkDPos;
  int               fNTrkUPos;
  int               fNTrkUpstream;
  int               fNStrawHits;
  int               fNCalHits;
  int               fNGenp;		// N(generated particles)

  int               fNHyp;
  int               fBestHyp[10];
  int               fFillDioHist;
					// fTrackNumber[i]: track number, 
					// corresponding to OBSP particle #i
					// or -1
  TStnArrayI        fTrackNumber;

  TStnTrack*        fTrack;
  TStnCluster*      fCluster;

  TDiskCalorimeter* fDiskCalorimeter;

  int               fBestID;

  Error_t           fError[kMaxNErrors];

  TString           fTrackBlockName;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TCosmicsAnaModule(const char* name="CosmicsAna", const char* title="CosmicsAna");
  ~TCosmicsAnaModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Hist_t*            GetHist        () { return &fHist;         }
  TStnTrackBlock*    GetTrackBlock  () { return fTrackBlockDem; }
  TStnClusterBlock*  GetClusterBlock() { return fClusterBlock;  }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void               SetTrackBlockName(const char* Name) { fTrackBlockName = Name; }
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
  void    BookCaloHistograms    (HistBase_t*  Hist, const char* Folder);
  void    BookClusterHistograms (HistBase_t*  Hist, const char* Folder);

  //  void    BookTrackIDHistograms (TStnTrackID::Hist_t* Hist, const char* Folder);

  void    FillCaloHistograms     (HistBase_t*  Hist, TStnCrystal*  Crystal);
  void    FillClusterHistograms  (HistBase_t*  Hist, TStnCluster*  Cluster);

  void    BookHistograms();
  void    FillHistograms();

  void    Debug();
//-----------------------------------------------------------------------------
// test
//-----------------------------------------------------------------------------
  void    Test001();

  ClassDef(TCosmicsAnaModule,0)
};
}

#endif
