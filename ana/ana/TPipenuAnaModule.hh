///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef murat_ana_TPipenuAnaModule_hh
#define murat_ana_TPipenuAnaModule_hh

#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TStnTrackSeedBlock.hh"
#include "Stntuple/obj/TStnHelixBlock.hh"
#include "Stntuple/obj/TTrackStrawHitBlock.hh"
#include "Stntuple/obj/TStnClusterBlock.hh"
#include "Stntuple/obj/TCalDataBlock.hh"
#include "Stntuple/obj/TStrawHitBlock.hh"
#include "Stntuple/obj/TStepPointMCBlock.hh"
#include "Stntuple/obj/TGenpBlock.hh"
#include "Stntuple/obj/TSimpBlock.hh"

#include "Stntuple/base/TStnArrayI.hh"

#include "Stntuple/geom/TDiskCalorimeter.hh"
#include "Stntuple/geom/TStnCrystal.hh"

//-----------------------------------------------------------------------------
// among other things, includes definitions of various histogram structures - 
// TrackHist_t etc
//-----------------------------------------------------------------------------
#include "murat/ana/TAnaModule.hh"

namespace murat {
class TPipenuAnaModule: public murat::TAnaModule {
public:
  enum { kNDisks = 2 };
//-----------------------------------------------------------------------------
//  histograms
//-----------------------------------------------------------------------------
  struct CaloHist_t {
    TH1F*    fDiskID;                  // per crystal hit
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
    TH1F*    fNCr0;                     // all clustered
    TH1F*    fNCr1;                     // above 1MeV
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


  struct TrackEffHist_t {
    TH1F*    fPtMc;                     // denominator
    TH1F*    fPtReco;                   // numerator
  };
//-----------------------------------------------------------------------------
  enum { kNEventHistSets         = 100 };
  enum { kNTrackHistSets         = 400 };
  enum { kNHelixHistSets         = 100 };
  enum { kNTrackStrawHitHistSets =  10 };
  enum { kNClusterHistSets       = 100 };
  enum { kNCaloHistSets          = 100 };
  enum { kNGenpHistSets          = 100 };
  enum { kNSimpHistSets          = 100 };

  struct Hist_t {
    TH1F*                  fCrystalR[2];                  // crystal radius
    CaloHist_t*            fCalo   [kNCaloHistSets   ];
    ClusterHist_t*         fCluster[kNClusterHistSets];
    EventHist_t*           fEvent  [kNEventHistSets  ];
    GenpHist_t*            fGenp   [kNGenpHistSets   ];
    HelixHist_t*           fHelix  [kNHelixHistSets  ];
    SimpHist_t*            fSimp   [kNSimpHistSets   ];
    TrackHist_t*           fTrack  [kNTrackHistSets  ];
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
  TString              fTrackBlockName;      // 
  TString              fTrackStrawHitBlockName; // 
                                             // pointers to the data blocks used
  TStnTrackBlock*      fTrackBlock;
  TStnTrackSeedBlock*  fTrackSeedBlock;
  TStnHelixBlock*      fHelixBlock;
  TTrackStrawHitBlock* fTrackStrawHitBlock;
  TStnClusterBlock*    fClusterBlock;
  TCalDataBlock*       fCalDataBlock;
  TStrawHitBlock*      fStrawHitBlock;
  TStepPointMCBlock*   fVDetBlock;
  TGenpBlock*          fGenpBlock;
  TSimpBlock*          fSimpBlock;
                                        // additional parameters (assume nhelices/ntracks < 20)
  HelixPar_t           fHelixPar[20];
  TrackPar_t           fTrackPar[20];

  SimPar_t             fSimPar;
                                        // histograms filled
  Hist_t               fHist;

  //  TGenParticle*        fParticle;       // electron or muon
  int                  fDirection;         // 1:downstream, -1:upstream  [direction of the particle]

  //  TSimParticle*        fSimp;
  double               fEleE;           // electron energy

  int                  fCalorimeterType;

  int                  fNClusters;
  int                  fNCalPatRec;
  int                  fNMatchedTracks;
  int                  fNStrawHits;
  int                  fNCalHits;
  int                  fNGenp;          // N(generated particles)

  int                  fNHyp;
  int                  fBestHyp[10];
                                        // fTrackNumber[i]: track number, 
                                        // corresponding to OBSP particle #i
                                        // or -1
  TStnArrayI           fTrackNumber;

  int                  fBestID;

  TStnTrack*           fTrack;
  TStnCluster*         fCluster;

  TDiskCalorimeter*    fDiskCalorimeter;

  double               fMinETrig;       // minimal energy of the cluster to trigger on
                                        // Tcm - track-cluster matching
  double               fMinDtTcm;
  double               fMaxDtTcm;

  double               fLumWt;
  int                  fApplyCorrections;  // 0: do not apply momentum ant DT corrections
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TPipenuAnaModule(const char* name="PipenuAna", const char* title="PipenuAna");
  ~TPipenuAnaModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  // Hist_t*            GetHist        () { return &fHist;        }
  // TStnTrackBlock*    GetTrackBlock  () { return fTrackBlock;   }
  // TStnClusterBlock*  GetClusterBlock() { return fClusterBlock; }

  // TStnTrackID*       GetTrackID(int I) { return fTrackID[I];   }
  // TEmuLogLH*         GetLogLH       () { return fLogLH;        }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void               SetDirection    (int Dir  ) { fDirection     = Dir  ; }

  void               SetTrackBlockName         (const char* Name) { fTrackBlockName         = Name; }
  void               SetTrackStrawHitBlockName (const char* Name) { fTrackStrawHitBlockName = Name; }

  void               SetApplyCorrections(int YesNo) { fApplyCorrections = YesNo; }
//-----------------------------------------------------------------------------
// overloaded methods of TStnModule
//-----------------------------------------------------------------------------
  virtual int     BeginJob();
  // virtual int     BeginRun();   // use TAnaModule::BeginRun()
  virtual int     Event   (int ientry);
  virtual int     EndJob  ();
//-----------------------------------------------------------------------------
// other methods
//-----------------------------------------------------------------------------
  void    BookCaloHistograms      (CaloHist_t*    Hist, const char* Folder);
  void    BookClusterHistograms   (ClusterHist_t* Hist, const char* Folder);
  void    FillCaloHistograms      (CaloHist_t*    Hist, TStnCrystal*  Crystal);
  void    FillClusterHistograms   (ClusterHist_t* Hist, TStnCluster*  Cluster);

  void    FillEfficiencyHistograms(TStnTrackBlock* TrackBlock, TStnTrackID* TrackID, int HistSet);

  void    BookHistograms();
  void    FillHistograms();

  void    Debug();
//-----------------------------------------------------------------------------
// test
//-----------------------------------------------------------------------------
  void    Test001();

  ClassDef(murat::TPipenuAnaModule,0)
};
}
#endif
