///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef murat_ana_TTrackCompModule_hh
#define murat_ana_TTrackCompModule_hh

#include "murat/ana/TAnaModule.hh"

#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TStnClusterBlock.hh"
#include "Stntuple/obj/TCalDataBlock.hh"
#include "Stntuple/obj/TStrawDataBlock.hh"
#include "Stntuple/obj/TGenpBlock.hh"
#include "Stntuple/obj/TSimpBlock.hh"
#include "Stntuple/obj/TStnTrackSeedBlock.hh"
#include "Stntuple/obj/TStnHelixBlock.hh"
#include "Stntuple/obj/TStepPointMCBlock.hh"

#include "Stntuple/base/TStnArrayI.hh"

#include "Stntuple/geom/TDiskCalorimeter.hh"

#include "Stntuple/alg/TStnTrackID.hh"
#include "Stntuple/alg/TEmuLogLH.hh"

namespace murat {
class TTrackCompModule: public murat::TAnaModule {
public:
					//  parameters
  struct DTrackHist_t : public HistBase_t {
    TH1F*    fDp;			// difference between the two track momenta
    TH1F*    fRMomErr10;
    TH1F*    fDT0;                      // shift in T0 
  };
//-----------------------------------------------------------------------------
//  fTrackHist[  0]: all tracks
//  fTrackHist[100]: Set C tracks
//-----------------------------------------------------------------------------
  enum { kPAR = 0, kDAR = 1 };

  enum { kNTrackPar          =    10 };

  enum { kNEventHistSets     =   100 };
  enum { kNTrackSeedHistSets =   100 };
  enum { kNTrackHistSets     = 10000 };
  enum { kNDTrackHistSets    =   100 };
  enum { kNSimpHistSets      =   100 };

  struct Hist_t {
    EventHist_t*     fEvent    [kNEventHistSets];
    TrackHist_t*     fTrack    [kNTrackHistSets];
    TrackSeedHist_t* fTrackSeed[kNTrackSeedHistSets];
    DTrackHist_t*    fDTrack   [kNDTrackHistSets];
  };


  struct Cut_t {
    double fXMin;
    double fXMax;
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used
  TStnTrackBlock*     fTrackBlock[2];	// [0]: TrkPatRec fit, [1]:CalPatRec fit
  TStnClusterBlock*   fClusterBlock;
  TGenpBlock*         fGenpBlock;
  TSimpBlock*         fSimpBlock;
  TStnTrackSeedBlock* fTrackSeedBlock;
  TStnHelixBlock*     fHelixBlock;
  TStepPointMCBlock*  fSpmcBlockVDet;

  TString             fTrackBlockName[2];
					
  TrackPar_t          fTrackPar[2][kNTrackPar];	// additional track parameters (assume ntracks < 10)
  SimPar_t            fSimPar;		        // additional parameters of the simulated MC particle
  Hist_t              fHist;		        // histograms filled

					        // cut values
  double              fPtMin;

  Cut_t               fDebugCut[100];

  int                fNHelices;         // 
  int                fNTrackSeeds;	// N reconstructed track seeds in the event
  int                fNGoodTrackSeeds;  // N track seeds passing quality cuts
  int                fNTracks    [2];	// 0:TrkPatRec 1:CalPatRec
  int                fNGoodTracks[2];
  int                fNSimp;		// N(simulated particles)

  TStnTrack*         fTrack;

  double             fMinETrig;

  double             fMinT0;
  double             fLumWt;

  TStnCluster*       fCluster;
  int                fNClusters;
  double             fEClMax;
  double             fTClMax;

  double             fKMaxRMC;             // RMC: closure approximation kMax

  int                fProcess;
  double             fPhotonE;

  double             fWeight;
//-----------------------------------------------------------------------------
// Track quality (TQ) MVA training ntuples
//-----------------------------------------------------------------------------
  int                 fWriteTrkQualMvaTree   ;       // =0:PAR   1:DAR  (default=-1)

  TFile*              fTrkQualMvaFile;
  TTree*              fTrkQualMvaTree;

  TmvaTrainData_t     fTrkQualMvaData;
  TmvaTrainBranches_t fTrkQualMvaBranch;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TTrackCompModule(const char* name="TrackComp", const char* title="TrackComp");
  ~TTrackCompModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Hist_t*            GetHist      () { return &fHist;      }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void    SetTrackBlockName(int I, const char* Name) { fTrackBlockName[I] = Name; }

  void    SetKMaxRMC      (double KMax) { fKMaxRMC = KMax ; }

  void    SetDebugCut(int I, double XMin, double XMax) {
    fDebugCut[I].fXMin = XMin;
    fDebugCut[I].fXMax = XMax;
  }

  void    SetWriteMvaTree(int Alg) { fWriteTrkQualMvaTree = Alg; }
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
  void    BookDTrackHistograms   (HistBase_t*   Hist, const char* Folder);

  void    FillDTrackHistograms   (HistBase_t*  Hist, TStnTrack* Trk1, TrackPar_t* Tp1, TStnTrack* Trk2, TrackPar_t* Tp2);

  void    BookHistograms();
  void    FillHistograms();

  void    Debug();
//-----------------------------------------------------------------------------
// test
//-----------------------------------------------------------------------------
  void    Test001();

  ClassDef(murat::TTrackCompModule,0)
};
}
#endif
