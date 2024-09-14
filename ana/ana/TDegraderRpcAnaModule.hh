///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef murat_ana_TDegraderRpcAnaModule_hh
#define murat_ana_TDegraderRpcAnaModule_hh

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
class TDegraderRpcAnaModule: public murat::TAnaModule {
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

  struct DRpcHist_t {
    TH1F*    fNHitsVD09;
    TH1F*    fSMomVD09[2];                    //
    TH1F*    fNHitsVD10;
    TH1F*    fSMomVD10[2];                    //
    TH1F*    fNHitsVD13;
    TH1F*    fSMomVD13[2];                    //
    TH2F*    fSMomVD13VsSinTh;             //
    TH1F*    fCPath;
    TH2F*    fSMomVD13VsCPath;
  };

  struct T2Par_t {
    int    fItrk1;
    int    fItrk2;
    float  fSMom;
    float  fQ;
  };
  
  struct T2Hist_t {
    TH1F*    fSMom[2];                     //
    TH1F*    fQ;                        //
  };

//-----------------------------------------------------------------------------
  enum { kNEventHistSets         = 100 };
  enum { kNTrackHistSets         = 100 };
  enum { kNHelixHistSets         = 100 };
  enum { kNTrackStrawHitHistSets =  10 };
  enum { kNClusterHistSets       = 100 };
  enum { kNCaloHistSets          = 100 };
  enum { kNGenpHistSets          = 100 };
  enum { kNSimpHistSets          = 100 };
  enum { kNTimeClusterHistSets   =  10 };
  enum { kNDRpcHistSets          = 400 };
  enum { kNTrackIDHistSets       =  10 };
  enum { kNT2HistSets            =  10 };

  struct Hist_t {
    ClusterHist_t*         fCluster    [kNClusterHistSets];
    EventHist_t*           fEvent      [kNEventHistSets  ];
    GenpHist_t*            fGenp       [kNGenpHistSets   ];
    HelixHist_t*           fHelix      [kNHelixHistSets  ];
    SimpHist_t*            fSimp       [kNSimpHistSets   ];
    TimeClusterHist_t*     fTimeCluster[kNTimeClusterHistSets];
    TrackHist_t*           fTrack      [kNTrackHistSets  ];
    DRpcHist_t*            fDRpc       [kNDRpcHistSets   ];
    TStnTrackID::Hist_t*   fTrackID    [kNTrackIDHistSets];
    T2Hist_t*              fT2         [kNT2HistSets  ];
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
  T2Par_t              fT2Par  [100];

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
  int                  fNT2Par;

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

  TStepPointMC*        fHitVD09[100];
  int                  fNHitsVD09;
  float                fSMomVD09;
  int                  fQVD09;

  TStepPointMC*        fHitVD10[100];
  int                  fNHitsVD10;
  float                fSMomVD10;
  int                  fQVD10;

  TStepPointMC*        fHitVD13[100];
  int                  fNHitsVD13;
  float                fSMomVD13;
  int                  fQVD13;

  double               fCPath;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TDegraderRpcAnaModule(const char* name="murat_DrpcAna", const char* title="DegraderRpcAna");
  ~TDegraderRpcAnaModule();
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
  void    BookT2Histograms(T2Hist_t* Hist, const char* Folder);
  void    BookDRpcHistograms(DRpcHist_t*      Hist, const char* Folder);

  // void    FillCaloHistograms       (CaloHist_t*        Hist, TStnCrystal*     Crystal);
  void    FillTimeClusterHistograms(TimeClusterHist_t* Hist, TStnTimeCluster* Tc     , double Weight = 1.);

  void    FillDRpcHistograms    (DRpcHist_t*    Hist, 
                                 double         Weight = 1.);
  void    FillT2Histograms      (T2Hist_t*      Hist, T2Par_t* T2Par, double Weight = 1.);

  void    BookHistograms();
  void    FillHistograms();

  void    Debug();
//-----------------------------------------------------------------------------
// test
//-----------------------------------------------------------------------------
  void    Test001();

  ClassDef(murat::TDegraderRpcAnaModule,0)
};
}
#endif
