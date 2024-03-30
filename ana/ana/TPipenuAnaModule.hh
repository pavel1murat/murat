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

  struct TimeClusterHist_t {
    TH1F*    fNsh;
    TH1F*    fNch;
    TH1F*    fT0;
    TH1F*    fT0Err;
  };

  struct TrackEffHist_t {
    TH1F*    fPtMc;                     // denominator
    TH1F*    fPtReco;                   // numerator
  };

  struct PipenuHist_t {
    TH2F*    fPVsT0;                    //
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
  enum { kNTimeClusterHistSets   =  10 };
  enum { kNPipenuHistSets        = 400 };

  struct Hist_t {
    TH1F*                  fCrystalR   [2];                  // crystal radius
    CaloHist_t*            fCalo       [kNCaloHistSets   ];
    ClusterHist_t*         fCluster    [kNClusterHistSets];
    EventHist_t*           fEvent      [kNEventHistSets  ];
    GenpHist_t*            fGenp       [kNGenpHistSets   ];
    HelixHist_t*           fHelix      [kNHelixHistSets  ];
    SimpHist_t*            fSimp       [kNSimpHistSets   ];
    TimeClusterHist_t*     fTimeCluster[kNTimeClusterHistSets];
    TrackHist_t*           fTrack      [kNTrackHistSets  ];
    PipenuHist_t*          fPipenu     [kNPipenuHistSets ];
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
  TString               fTrackBlockName;      // 
  TString               fTrackStrawHitBlockName; // 
                                             // pointers to the data blocks used
  TStnTrackBlock*       fTrackBlock;
  TStnTrackSeedBlock*   fTrackSeedBlock;
  TStnHelixBlock*       fHelixBlock;
  TStnTimeClusterBlock* fTimeClusterBlock;
  // TTrackStrawHitBlock* fTrackStrawHitBlock;
  TStnClusterBlock*     fClusterBlock;
  TCalDataBlock*        fCalDataBlock;
  // TStrawHitBlock*      fStrawHitBlock;
  TStepPointMCBlock*    fVDetBlock;
  TGenpBlock*           fGenpBlock;
  TSimpBlock*           fSimpBlock;
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

  int                  fEventPassedSelections; // 1: event passed analysis selections; 0: event didn't pass
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TPipenuAnaModule(const char* name="PipenuAna", const char* title="PipenuAna");
  ~TPipenuAnaModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
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
  void    BookCaloHistograms       (CaloHist_t*        Hist, const char* Folder);
  void    BookTimeClusterHistograms(TimeClusterHist_t* Hist, const char* Folder);
  void    BookPipenuHistograms     (PipenuHist_t*      Hist, const char* Folder);

  void    FillCaloHistograms       (CaloHist_t*        Hist, TStnCrystal*     Crystal);
  void    FillTimeClusterHistograms(TimeClusterHist_t* Hist, TStnTimeCluster* Tc     , double Weight = 1.);

  void    FillPipenuHistograms    (PipenuHist_t*  Hist, 
                                   TStnTrack*     Trk, 
				   TrackPar_t*    Tp, 
				   SimPar_t*      SimPar,
				   double         Weight = 1.);

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
