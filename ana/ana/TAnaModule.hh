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
#include "TMVA/Reader.h"

#include "Stntuple/loop/TStnModule.hh"

#include "Stntuple/obj/TStnTimeClusterBlock.hh"
#include "Stntuple/obj/TStnHelixBlock.hh"
#include "Stntuple/obj/TStnTrackSeedBlock.hh"
#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TStnClusterBlock.hh"
#include "Stntuple/obj/TCalDataBlock.hh"
#include "Stntuple/obj/TStepPointMCBlock.hh"
#include "Stntuple/obj/TStrawHitBlock.hh"
#include "Stntuple/obj/TGenpBlock.hh"
#include "Stntuple/obj/TSimpBlock.hh"
#include "Stntuple/obj/TCrvClusterBlock.hh"
#include "Stntuple/obj/TCrvPulseBlock.hh"

#include "Stntuple/base/TStnArrayI.hh"

#include "Stntuple/geom/TStnCrystal.hh"
#include "Stntuple/alg/TStnTrackID.hh"
#include "Stntuple/alg/TEmuLogLH.hh"
#include "Stntuple/alg/TStntuple.hh"

#include "murat/ana/EventPar_t.hh"
#include "murat/ana/SimPar_t.hh"
#include "murat/ana/SimpData_t.hh"
#include "murat/ana/HelixPar_t.hh"
#include "murat/ana/TrackPar_t.hh"

#include "murat/ana/ClusterHist_t.hh"
#include "murat/ana/CrvHist_t.hh"
#include "murat/ana/EventHist_t.hh"
#include "murat/ana/GenpHist_t.hh"
#include "murat/ana/SimpHist_t.hh"
#include "murat/ana/HelixHist_t.hh"
#include "murat/ana/TrackHist_t.hh"
#include "murat/ana/TrackSeedHist_t.hh"

#include "murat/ana/TTrkQualMva.hh"
#include "murat/ana/mva_data.hh"

namespace mu2e { 
  class MVATools;
}

namespace murat {

class TAnaModule: public TStnModule {
public:
					// make these definitions internal for the module to avoid 
					// dealing with namespaces
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:

  enum { kMaxNTrackID = 50 };

  TStnTrackID*        fTrackID_BOX;
  TStnTrackID*        fTrackID_MVA;

  int                 fNID;                    // number of different track ID's used, default=0
  TStnTrackID*        fTrackID[kMaxNTrackID];

  double              fMbTime;                 // microbunch time
  double              fMinT0;
  int                 fApplyCorr;              // default=1

  int                 fPDGCode;		       // determines which one
  int                 fMCProcessCode;      
  double              fEleE;                   // relevant event generated energy
  double              fEventWeight;            // weight applied when filling histograms
  int                 fBatchMode;

  EventPar_t          fEvtPar;

  TStntuple*          fStnt;                   // STNTUPLE singleton
//-----------------------------------------------------------------------------
// TrackQual MVA
//-----------------------------------------------------------------------------
  int                 fUseTrqMVA;
  int                 fNTrqMVA;	               // number of MVA classifiers used for tracks of the same type
  float               fTrqData[20];            // with a few spare words

  int                 fUsePidMVA;
  int                 fNPidMVA;	               // number of MVA classifiers used for tracks of the same type
  mva_data*           fPidMVA[20];

  mva_data*           fPidReader;              // ROOT-based implementation
  float               fPidData[20];            // with a few spare words
//-----------------------------------------------------------------------------
// MVA-based classifiers
// so far, use up to 2. for PAR and DAR tracks separately, [0]:PAR, [1]:DAR
//-----------------------------------------------------------------------------
  mva_data*           fTrqMVA[20];
  int                 fBestID[2];          // best track ID for two ambig resolvers
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TAnaModule(const char* name="MuratAna", const char* title="MuratAna");
  ~TAnaModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  TStnTrackID*       GetTrackID(int I) { return fTrackID[I]; }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void               SetMinT0         (double T0  ) { fMinT0         = T0   ; }
  void               SetPDGCode       (int    Code) { fPDGCode       = Code ; }
  void               SetMCProcessCode (int    Code) { fMCProcessCode = Code ; }

  void               SetSignalParticle(int PDGCode, int MCProcessCode) { 
    fPDGCode       = PDGCode; 
    fMCProcessCode = MCProcessCode; 
  }

  void               SetApplyCorr     (int    Flag) { fApplyCorr     = Flag ; }
//-----------------------------------------------------------------------------
// TRQ MVA Training Codes: 
//
// 0060 : PAR dPf > 0.60
// 0070 : PAR dPf > 0.70
// 1060 : DAR dPf > 0.60
// 1070 : DAR dPf > 0.70
//-----------------------------------------------------------------------------
  void                SetTrqMVA      (int Block, const char* Dataset, int MvaTrainingCode);
//-----------------------------------------------------------------------------
// PID MVA
//-----------------------------------------------------------------------------
  void                SetPidMVA      (const char* Dataset, int MvaTrainingCode);
//-----------------------------------------------------------------------------
// overloaded methods of TStnModule
//-----------------------------------------------------------------------------
  virtual int     BeginJob();
  virtual int     BeginRun();
  // virtual int     Event   (int ientry);
  virtual int     EndJob  ();
//-----------------------------------------------------------------------------
// other methods
//-----------------------------------------------------------------------------
  virtual void    BookClusterHistograms   (HistBase_t*          Hist, const char* Folder);
  virtual void    BookCrvClusterHistograms(CrvClusterHist_t*    Hist, const char* Folder);
  virtual void    BookCrvPulseHistograms  (CrvPulseHist_t*      Hist, const char* Folder);
  virtual void    BookGenpHistograms      (GenpHist_t*          Hist, const char* Folder);
  virtual void    BookEventHistograms     (HistBase_t*          Hist, const char* Folder);
  virtual void    BookHelixHistograms     (HistBase_t*          Hist, const char* Folder);
  virtual void    BookSimpHistograms      (HistBase_t*          Hist, const char* Folder);
  virtual void    BookTrackHistograms     (TrackHist_t*         Hist, const char* Folder);
  virtual void    BookTrackIDHistograms   (TStnTrackID::Hist_t* Hist, const char* Folder);
  virtual void    BookTrackSeedHistograms (HistBase_t*          Hist, const char* Folder);

  virtual void    FillClusterHistograms   (HistBase_t*          Hist, TStnCluster*          Cl   , double Weight = 1.);
  virtual void    FillCrvClusterHistograms(CrvClusterHist_t*  Hist, TCrvCoincidenceCluster* CrvCl);
  virtual void    FillCrvPulseHistograms  (CrvPulseHist_t*    Hist, TCrvRecoPulse*          Pulse);

  virtual void    FillEventHistograms     (HistBase_t*  Hist, EventPar_t*  Evtpar  );

  virtual void    FillGenpHistograms      (GenpHist_t*   Hist, TGenParticle* Genp   );

  virtual void    FillHelixHistograms     (HistBase_t*  Hist, TStnHelix* Hel, HelixPar_t* Help, double Weight = 1);

  virtual void    FillSimpHistograms      (HistBase_t*  Hist, TSimParticle*  Simp   , SimpData_t* Sd, double Weight = 1.);
  virtual void    FillTrackSeedHistograms (HistBase_t*  Hist, TStnTrackSeed* TrkSeed, double Weight = 1);

  virtual void    FillTrackHistograms     (TrackHist_t* Hist, 
                                           TStnTrack*   Trk, 
                                           TrackPar_t*  Tp, 
                                           SimPar_t*    SimPar,
                                           double       Weight = 1.);

  virtual int     InitTrackPar(TStnTrackBlock*    TrackBlock  , 
		       TStnClusterBlock*  ClusterBlock, 
		       TrackPar_t*        TrackPar    ,
		       SimPar_t*          SimPar      );
  
  virtual void    PrintTrack(TStnTrack* Track, TrackPar_t* Tp, Option_t* Option) const ;

  virtual double  BatchModeWeight(float lumi, int mode);

  ClassDef(murat::TAnaModule,0)
};
}
#endif
