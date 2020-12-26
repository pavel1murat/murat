///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef __murat_ana_TAnaModule_hh__
#define __murat_ana_TAnaModule_hh__

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TTree.h"
#include "TBranch.h"


#include "Stntuple/loop/TStnModule.hh"

#include "Stntuple/obj/TStnTimeClusterBlock.hh"
#include "Stntuple/obj/TStnHelixBlock.hh"
#include "Stntuple/obj/TStnTrackSeedBlock.hh"
#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TStnClusterBlock.hh"
#include "Stntuple/obj/TCalDataBlock.hh"
#include "Stntuple/obj/TStepPointMCBlock.hh"
#include "Stntuple/obj/TStrawDataBlock.hh"
#include "Stntuple/obj/TGenpBlock.hh"
#include "Stntuple/obj/TSimpBlock.hh"
#include "Stntuple/obj/TCrvClusterBlock.hh"
#include "Stntuple/obj/TCrvPulseBlock.hh"

#include "Stntuple/base/TStnArrayI.hh"

#include "Stntuple/geom/TStnCrystal.hh"
#include "Stntuple/alg/TStnTrackID.hh"
#include "Stntuple/alg/TEmuLogLH.hh"
#include "Stntuple/alg/TStntuple.hh"

#include "murat/ana/TrackHist_t.hh"

#include "murat/ana/EventPar_t.hh"
#include "murat/ana/SimPar_t.hh"
#include "murat/ana/TrackPar_t.hh"

#include "murat/ana/CrvHist_t.hh"
#include "murat/ana/EventHist_t.hh"
#include "murat/ana/GenpHist_t.hh"
#include "murat/ana/SimpHist_t.hh"
#include "murat/ana/TrackSeedHist_t.hh"

#include "murat/ana/mva_data.hh"

namespace mu2e { 
  class MVATools;
};

namespace murat {

class TAnaModule: public TStnModule {
public:
					// make these definitions internal for the module to avoid 
					// dealing with namespaces
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:

  struct TmvaTrainData_t {
    float    fP;
    float    fPMC;
    float    fTanDip;
    float    fNActive;
    float    fNaFract;
    float    fNDoublets;                // total number of doublets
    float    fNDa;                      // number of doublets with all hits active
    float    fChi2Dof;
    float    fFitCons;
    float    fMomErr;
    float    fT0Err;
    float    fD0;
    float    fRMax;
    float    fNdaOverNa;
    float    fNzaOverNa;
    float    fNmaOverNm;
    float    fNdaOverNd;
    float    fZ1;			// Z-coordinate of the first hit
    float    fWeight;
  };

  struct TmvaTrainBranches_t {
    TBranch*  fP;
    TBranch*  fPMC;
    TBranch*  fTanDip;
    TBranch*  fNActive;
    TBranch*  fNaFract;
    TBranch*  fNDoublets;
    TBranch*  fNDa;
    TBranch*  fChi2Dof;
    TBranch*  fFitCons;
    TBranch*  fMomErr;
    TBranch*  fT0Err;
    TBranch*  fD0;
    TBranch*  fRMax;
    TBranch*  fNdaOverNa;
    TBranch*  fNzaOverNa;
    TBranch*  fNmaOverNm;
    TBranch*  fNdaOverNd;
    TBranch*  fZ1;			       // Z-coordinate of the first hit
    TBranch*  fWeight;			       // for background only
  };

  enum { kMaxNTrackID = 50 };

  TStnTrackID*        fTrackID_BOX;
  TStnTrackID*        fTrackID_MVA;

  int                 fNID;                    // number of different track ID's used, default=0
  TStnTrackID*        fTrackID[kMaxNTrackID];

  double              fMbTime;                 // microbunch time
  double              fMinT0;

  int                 fPdgCode;		       // determines which one
  int                 fGeneratorCode;      

  TEmuLogLH*          fLogLH;
  EventPar_t          fEvtPar;

  TStntuple*          fStnt;                   // STNTUPLE singleton
//-----------------------------------------------------------------------------
// TMVA training ntuples
//-----------------------------------------------------------------------------
  int                 fWriteTmvaTree   ;       // =0:PAR   1:DAR  (default=-1)

  TFile*              fTmvaFile;
  TTree*              fTmvaTree;

  TmvaTrainData_t     fTmvaData;
  TmvaTrainBranches_t fTmvaBranch;

  int                 fUseMVA;
  int                 fNMVA;	           // number of MVA classifiers used for tracks of the same type

  int                 fApplyCorr;          // default=1
//-----------------------------------------------------------------------------
// MVA-based classifiers
// so far, use up to 2. for PAR and DAR tracks separately, [0]:PAR, [1]:DAR
//-----------------------------------------------------------------------------
  mva_data*           fTrkQualMVA[20];
  int                 fBestID[2];          // best track ID for two ambig resolvers
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TAnaModule(const char* name="su2020_Ana", const char* title="Ana");
  ~TAnaModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  TEmuLogLH*         GetLogLH       () { return fLogLH; }
  TStnTrackID*       GetTrackID(int I) { return fTrackID[I]; }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void               SetMinT0        (double T0  ) { fMinT0         = T0   ; }
  void               SetPdgCode      (int    Code) { fPdgCode       = Code ; }
  void               SetGeneratorCode(int    Code) { fGeneratorCode = Code ; }
  void               SetApplyCorr    (int    Flag) { fApplyCorr     = Flag ; }
  void               SetWriteTmvaTree(int    Algo) { fWriteTmvaTree = Algo ; }
//-----------------------------------------------------------------------------
// MVA Training Codes: 
//
// 0060 : PAR dPf > 0.60
// 0070 : PAR dPf > 0.70
// 1060 : DAR dPf > 0.60
// 1070 : DAR dPf > 0.70
//-----------------------------------------------------------------------------
  void                SetMVA         (const char* Dataset, int MvaTrainingCode);
//-----------------------------------------------------------------------------
// overloaded methods of TStnModule
//-----------------------------------------------------------------------------
  virtual int     BeginJob();
  // virtual int     BeginRun();
  // virtual int     Event   (int ientry);
  virtual int     EndJob  ();
//-----------------------------------------------------------------------------
// other methods
//-----------------------------------------------------------------------------
  void    BookCrvClusterHistograms(CrvClusterHist_t* Hist, const char* Folder);
  void    BookCrvPulseHistograms  (CrvPulseHist_t*   Hist, const char* Folder);
  void    BookGenpHistograms      (GenpHist_t*       Hist, const char* Folder);
  void    BookEventHistograms     (EventHist_t*      Hist, const char* Folder);
  void    BookSimpHistograms      (SimpHist_t*       Hist, const char* Folder);
  void    BookTrackHistograms     (TrackHist_t*         Hist, const char* Folder);
  void    BookTrackIDHistograms   (TStnTrackID::Hist_t* Hist, const char* Folder);
  void    BookTrackSeedHistograms (HistBase_t*          HistR, const char* Folder);


  void    FillCrvClusterHistograms(CrvClusterHist_t*  Hist, TCrvCoincidenceCluster* CrvCl);
  void    FillCrvPulseHistograms  (CrvPulseHist_t*    Hist, TCrvRecoPulse* CrvCl);

  void    FillEventHistograms     (EventHist_t*  Hist, EventPar_t*  Evtpar  );

  void    FillGenpHistograms      (GenpHist_t*   Hist, TGenParticle* Genp   );
  void    FillSimpHistograms      (SimpHist_t*   Hist, TSimParticle* Simp   );

  void    FillTrackSeedHistograms (HistBase_t*  HistR, TStnTrackSeed* TrkSeed);

  void    FillTrackHistograms     (TrackHist_t* Hist, 
				   TStnTrack*          Trk, 
				   TrackPar_t*         Tp, 
				   SimPar_t*           SimPar,
				   double              Weight = 1.);

  int     InitTrackPar(TStnTrackBlock*    TrackBlock  , 
		       TStnClusterBlock*  ClusterBlock, 
		       TrackPar_t*        TrackPar    ,
		       SimPar_t*          SimPar      );
  
  void    PrintTrack(TStnTrack* Track, TrackPar_t* Tp, Option_t* Option) const ;

  //  ClassDef(TAnaModule,0)
};
}
#endif
