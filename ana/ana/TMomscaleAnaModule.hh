///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef murat_ana_TMomscaleAnaModule_hh
#define murat_ana_TMomscaleAnaModule_hh

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

#include "Stntuple/alg/TStnTrackID.hh"

#include "Stntuple/base/TStnArrayI.hh"

#include "Stntuple/geom/TDiskCalorimeter.hh"
#include "Stntuple/geom/TStnCrystal.hh"

//-----------------------------------------------------------------------------
// among other things, includes definitions of various histogram structures - 
// TrackHist_t etc
//-----------------------------------------------------------------------------
#include "murat/ana/TAnaModule.hh"

namespace murat {
class TMomscaleAnaModule: public murat::TAnaModule {
public:
  enum { kNDisks = 2 };
//-----------------------------------------------------------------------------
//  histograms
//-----------------------------------------------------------------------------
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

  struct MomscaleHist_t {
    TH2F*    fEpVsP;                   //
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
  enum { kNMomscaleHistSets        = 400 };
  enum { kNTrackIDHistSets       =  10 };

  struct Hist_t {
    ClusterHist_t*         fCluster    [kNClusterHistSets];
    EventHist_t*           fEvent      [kNEventHistSets  ];
    GenpHist_t*            fGenp       [kNGenpHistSets   ];
    HelixHist_t*           fHelix      [kNHelixHistSets  ];
    SimpHist_t*            fSimp       [kNSimpHistSets   ];
    TimeClusterHist_t*     fTimeCluster[kNTimeClusterHistSets];
    TrackHist_t*           fTrack      [kNTrackHistSets  ];
    MomscaleHist_t*          fMomscale     [kNMomscaleHistSets ];
    TStnTrackID::Hist_t*   fTrackID    [kNTrackIDHistSets];
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
  TStnClusterBlock*     fClusterBlock;
  TCalDataBlock*        fCalDataBlock;
  TStepPointMCBlock*    fVDetBlock;
  TGenpBlock*           fGenpBlock;
  TSimpBlock*           fSimpBlock;
                                        // additional parameters (assume nhelices/ntracks < 20)
  HelixPar_t           fHelixPar[20];
  TrackPar_t           fTrackPar[20];

  SimPar_t             fSimPar;
                                        // histograms filled
  Hist_t               fHist;

  int                  fDirection;         // 1:downstream, -1:upstream  [direction of the particle]

  double               fEleE;           // electron energy

  int                  fCalorimeterType;

  int                  fNClusters;
  int                  fNCalPatRec;
  int                  fNMatchedTracks;
  int                  fNStrawHits;
  int                  fNCalHits;
  int                  fNGenp;          // N(generated particles)

                                        // fTrackNumber[i]: track number, 
                                        // corresponding to OBSP particle #i
                                        // or -1
  TStnArrayI           fTrackNumber;

  int                  fBestID;

  TStnTrack*           fTrack;
  TStnCluster*         fCluster;

  double               fLumWt;
  int                  fApplyCorrections;  // 0: do not apply momentum ant DT corrections

  int                  fEventPassedSelections; // 1: event passed analysis selections; 0: event didn't pass
  int                  fDnMax;
  int                  fN300;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TMomscaleAnaModule(const char* name="MomscaleAna", const char* title="MomscaleAna");
  ~TMomscaleAnaModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void               SetDirection    (int Dir  ) { fDirection     = Dir  ; }
  void               SetDnMax        (int N    ) { fDnMax         = N    ; }

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
  // void    BookCaloHistograms       (CaloHist_t*        Hist, const char* Folder);
  void    BookTimeClusterHistograms(TimeClusterHist_t* Hist, const char* Folder);
  void    BookMomscaleHistograms   (MomscaleHist_t*      Hist, const char* Folder);

  // void    FillCaloHistograms       (CaloHist_t*        Hist, TStnCrystal*     Crystal);
  void    FillTimeClusterHistograms(TimeClusterHist_t* Hist, TStnTimeCluster* Tc     , double Weight = 1.);

  void    FillMomscaleHistograms    (MomscaleHist_t*  Hist, 
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

  ClassDef(murat::TMomscaleAnaModule,0)
};
}
#endif
