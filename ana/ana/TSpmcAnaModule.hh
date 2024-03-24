///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef murat_ana_TSpmcAnaModule_hh
#define murat_ana_TSpmcAnaModule_hh

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#include "Stntuple/loop/TStnModule.hh"

#include "Stntuple/obj/TGenpBlock.hh"
#include "Stntuple/obj/TSimpBlock.hh"
#include "Stntuple/obj/TStepPointMCBlock.hh"

#include "Stntuple/base/TStnArrayI.hh"

#include "Stntuple/alg/TStntuple.hh"

#include "murat/ana/HistBase_t.h"
#include "murat/ana/SimPar_t.hh"
#include "murat/ana/VDetData_t.hh"

#include "murat/ana/AnaDefs.hh"

namespace murat {

class TSpmcAnaModule: public TStnModule {
public:
  enum { kMaxNSimp = 1000 };
//-----------------------------------------------------------------------------
//  histograms
//-----------------------------------------------------------------------------
  struct EventHist_t : public HistBase_t {
    TH1F*      fRunNumber;
    TH1F*      fEventNumber;
    TH1F*      fNSimp;
    TH1F*      fTMaxSimp[2];
    TH1F*      fTMaxSpmc;
    TH1F*      fLogWeight;
  };

  struct SimpHist_t : public HistBase_t {
    TH1F*      fVolumeID;		       //
    TH1F*      fStage;
    TH1F*      fGeneratorID;
    TH1F*      fTime;
    TH1F*      fStageDt;
    TH1F*      fParentPDG;
    TH1F*      fParentMom;
    TH1F*      fStartMom[2];
    TH1F*      fNStrawHits;

    TH2F*      fYVsX;
    TH2F*      fXEndVsZEnd;
    TH2F*      fYcVsZEnd;
    TH2F*      fPVD9VsZEnd;
    TH2F*      fYVsX_2480;
    TH2F*      fYVsX_2513;
    TH2F*      fCosThVsMom[2];
    TH2F*      fCosThVsMomPV;		// for antiprotons
  };

  struct StepPointMCHist_t : public HistBase_t {
    TH1F*      fVolumeID;		       //
    TH1F*      fGenIndex;		       //
    TH1F*      fSimID;
    TH1F*      fPDGCode[2];  // just different ranges
    TH1F*      fCreationCode;
    TH1F*      fParentSimID;
    TH1F*      fParentPDGCode;
    TH1F*      fEndProcessCode;

    TH1F*      fEDepTot;
    TH1F*      fEDepNio;
    TH1F*      fTime;			// in ns
    TH1F*      fTimeSec;		// in seconds, for R/A decays
    TH1F*      fStepLength;

    TH1F*      fMom[2];
    TH2F*      fCosThVsMom[2];
    TH1F*      fEKin;

    TH2F*      fYVsX;			// in local coord system
    TH2F*      fYcVsXc;			// in local coord system
    TH2F*      fCosThVsMomPV;		// for antiprotons
  };

  struct VDetHist_t : public HistBase_t {
    TH1F*      fIndex   ;
    TH1F*      fPDGCode ;		       //
    TH1F*      fGenCode ;		       // generator code
    TH1F*      fMom[2]  ;
    TH1F*      fTime    ;
    TH1F*      fPTime   ;                // proper time
    TH2F*      fYVsX    ;                // local coordinates
    TH2F*      fYcVsXc  ;                // trajectory axis
    TH2F*      fYcVsP   ;                // Yc vs P
    TH1F*      fPt      ;                // transverse mom
    TH1F*      fPp      ;                // momentum component parallel to the solenoid axis
    TH1F*      fTanTh   ;		       // tan (pitch angle)
    TH1F*      fEKin    ;

    TH1F*      fDt1508;

    TH2F*      fCosThVsMom[2];	       // cos (pitch angle) vs Mom
    TH2F*      fCosThVsMomPV;		// for antiprotons
    TH2F*      fTimeVsMom;
    TH2F*      fTimeVsMomW;
    TH2F*      fPTimeVsMom;
  };

  struct SimpData_t {
    int           fIndex;		// so far, not used
    TSimParticle* fParent;              // parent (in production vertex)
    int           fStage;               // for this particle
    TStepPointMC* fStepVD9;
    float         fPVD9;                // P(VD9), defined for stopped particles
    float         fY0;
    double        fWeight;
    double        fTau;                 // lifetime , in ns
    bool          fMuonDecay;
    TSimParticle* fEle;
  };

//-----------------------------------------------------------------------------
// assume Zlocal is normal to the virtual detector plane, YLocal points upwards
//-----------------------------------------------------------------------------
  struct SpmcData_t {
    TParticlePDG*  fParticle;		//
    float          fQ;                  // charge, e- : -1
    float          fM;
    float          fP;
    float          fPtLoc;
    float          fPxLoc;
    float          fPyLoc;
    float          fPzLoc;
    float          fEKin;
    float          fCosTh;              // cos(pitch angle wrt the beamline axis)
    float          fTanTh;              // tan(pitch angle wrt the beamline axis)
    float          fXLoc;               // X in the locall coord system of the virtual detector
    float          fYLoc;               // Y in the locall coord system of the virtual detector
    float          fR;                  // trajectory radius
    float          fX0;                 // X coordinate of the trajectory axis in the local coordinate system
    float          fY0;			// Y coordinate of the trajectory axis in the local coordinate system
    float          fSurvivalProb;       // survival prob, calculate once
    float          fDt1508;             // deltaT(15-8)
  };

//-----------------------------------------------------------------------------
  enum { kNEventHistSets        =   100 };
  enum { kNStepPointMCHistSets  = 10000 };
  enum { kNSimpHistSets         = 10000 };
  enum { kNVDetHistSets         = 10000 };

  struct Hist_t {
    EventHist_t*        fEvent      [kNEventHistSets      ];
    SimpHist_t*         fSimp       [kNSimpHistSets       ];
    VDetHist_t*         fVDet       [kNVDetHistSets       ];
    StepPointMCHist_t*  fStepPointMC[kNStepPointMCHistSets];
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used
  TGenpBlock*           fGenpBlock;  
  TSimpBlock*           fSimpBlock;  
  TStepPointMCBlock*    fSpmcBlock;
  TStepPointMCBlock*    fVDetBlock;
					// histograms filled
  Hist_t                fHist;

  TString               fSpmcBlockName;
  TString               fVDetBlockName;

  TSimParticle*         fMuon;		// pointer to stopped muon (pend=0)
  TSimParticle*         fParent;
  TSimParticle*         fProton;

  TParticlePDG*         fParticleCache[5000];

  TDatabasePDG*         fPdgDb;

  int                   fNVDetHits  ;
  int                   fNVDet;
  VDetData_t            fVDet[200];
  int                   fStageID;
  int                   fNSimp;
  int                   fStage;

  SimpData_t            fSimData[kMaxNSimp];

  TStntuple*            fStnt;
  double                fWeight;         // event weight, determined by the production cross section
  double                fTMaxSimp;	 // in seconds
  double                fTMaxSpmc;	 // in ns
					 // antiproton-specific : in the production vertex
  double                fPbarCosThPV;
  double                fPbarMomPV;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TSpmcAnaModule(const char* name="murat_SpmcAna", const char* title="murat SpmcAna");
  ~TSpmcAnaModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Hist_t*            GetHist        () { return &fHist;        }
//-----------------------------------------------------------------------------
// setters
//-----------------------------------------------------------------------------
  void SetSpmcBlockName(const char* Name) { fSpmcBlockName = Name; }
  void SetVDetBlockName(const char* Name) { fVDetBlockName = Name; }
  void SetStageID      (int ID) { fStageID = ID; }

  void          SetParticleCache(int PdgCode, TParticlePDG* P) { fParticleCache[2500+PdgCode] = P; }
  TParticlePDG* GetParticleCache(int PdgCode) { return fParticleCache[2500+PdgCode]; }
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
  void    BookEventHistograms        (HistBase_t* Hist, const char* Folder);
  void    BookSimpHistograms         (HistBase_t* Hist, const char* Folder);
  void    BookStepPointMCHistograms  (HistBase_t* Hist, const char* Folder);
  void    BookVDetHistograms         (HistBase_t* Hist, const char* Folder);

  void    FillEventHistograms        (HistBase_t* Hist);
  void    FillSimpHistograms         (HistBase_t* Hist, TSimParticle* Simp, SimpData_t* SimpData, double Weight = 1.);
  void    FillStepPointMCHistograms  (HistBase_t* Hist, TStepPointMC* Step, SpmcData_t* SpmcData, double Weight = 1.);
  void    FillVDetHistograms         (HistBase_t* Hist, TStepPointMC* Step, SpmcData_t* SpmcData, double Weight = 1.);

  void    InitSpmcData               (TStepPointMC* Step, SpmcData_t* SpmcData);

  void    BookHistograms();
  void    FillHistograms();


  void    Debug();
//-----------------------------------------------------------------------------
// test
//-----------------------------------------------------------------------------
  void    Test001();

  ClassDef(TSpmcAnaModule,0)
};

}
#endif
