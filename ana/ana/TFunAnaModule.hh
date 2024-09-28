///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef murat_ana_TFunAnaModule_hh
#define murat_ana_TFunAnaModule_hh

#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TStepPointMCBlock.hh"
#include "Stntuple/obj/TSimpBlock.hh"

#include "Stntuple/base/TStnArrayI.hh"

//-----------------------------------------------------------------------------
// among other things, includes definitions of various histogram structures - 
// TrackHist_t etc
//-----------------------------------------------------------------------------
#include "murat/ana/TAnaModule.hh"
#include "murat/ana/SimpData_t.hh"

namespace murat {
class TFunAnaModule: public murat::TAnaModule {
public:
//-----------------------------------------------------------------------------
//  histograms
//-----------------------------------------------------------------------------
  struct FunHist_t {
    TH1F*    fTofVD13;
    TH2F*    fTofVD13VsCosth;
    TH2F*    fTofVD13VsZ0;
    TH2F*    fTofVD13VsPhi;
    TH2F*    fTofVD13VsRho;
  };
  
//-----------------------------------------------------------------------------
  enum { kNEventHistSets         = 100 };
  enum { kNTrackHistSets         = 400 };
  enum { kNSimpHistSets          = 100 };
  enum { kNFunHistSets           = 500 };

  struct Hist_t {
    EventHist_t*           fEvent  [kNEventHistSets  ];
    SimpHist_t*            fSimp   [kNSimpHistSets   ];
    TrackHist_t*           fTrack  [kNTrackHistSets  ];
    FunHist_t*             fFun    [kNFunHistSets    ];
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
  TString              fTrackBlockName;      // 
                                             // pointers to the data blocks used
  TStnTrackBlock*      fTrackBlock;
  TStepPointMCBlock*   fVDetBlock;
  TSimpBlock*          fSimpBlock;
                                        // additional parameters (assume nhelices/ntracks < 20)
  HelixPar_t           fHelixPar[20];
  TrackPar_t           fTrackPar[20];

  SimPar_t             fSimPar;
                                        // histograms filled
  Hist_t               fHist;

  TGenParticle*        fParticle;       // electron or muon
  int                  fDirection;         // 1:downstream, -1:upstream  [direction of the particle]

  TSimParticle*        fSimp;
  double               fEleE;           // electron energy

                                        // fTrackNumber[i]: track number, 
                                        // corresponding to OBSP particle #i
                                        // or -1
  TStnArrayI           fTrackNumber;

  int                  fBestID;

  TStnTrack*           fTrack;

  double               fMinETrig;       // minimal energy of the cluster to trigger on
                                        // Tcm - track-cluster matching
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TFunAnaModule(const char* name="murat_FunAna", const char* title="murat_FunAna");
  ~TFunAnaModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void               SetDirection    (int Dir  ) { fDirection     = Dir  ; }

  void               SetTrackBlockName         (const char* Name) { fTrackBlockName         = Name; }
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
  void    BookFunHistograms       (FunHist_t*     Hist, const char*  Folder);

  void    FillFunHistograms       (FunHist_t* Hist);

  void    BookHistograms();
  void    FillHistograms();

  void    Debug();
//-----------------------------------------------------------------------------
// test
//-----------------------------------------------------------------------------
  void    Test001();

  ClassDef(murat::TFunAnaModule,0)
};
}
#endif
