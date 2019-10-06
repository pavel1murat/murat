///////////////////////////////////////////////////////////////////////////////
// analysis of the distributions for muons stopped in the ST
///////////////////////////////////////////////////////////////////////////////
#ifndef murat_ana_TMuonStopAnaModule_hh
#define murat_ana_TMuonStopAnaModule_hh

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

#include "Stntuple/loop/TStnModule.hh"
#include "Stntuple/obj/TSimpBlock.hh"
#include "Stntuple/obj/TStepPointMCBlock.hh"

#include "murat/ana/VDetData_t.hh"
#include "murat/ana/HistBase_t.h"

#include "murat/ana/AnaDefs.hh"

class TMuonStopAnaModule: public TStnModule {
public:
//-----------------------------------------------------------------------------
//  histograms
//-----------------------------------------------------------------------------
  struct SimpHist_t : public HistBase_t {
    TH1F*    fVolumeID;		       //
    TH1F*    fGeneratorID;
    TH1F*    fTime;
    TH1F*    fParentPDG;
    TH1F*    fParentMom;
    TH1F*    fStartMom;

    TH2F*    fYVsX;
    TH2F*    fXEndVsZEnd;
    TH2F*    fYVsX_2480;
    TH2F*    fYVsX_2513;
  };

  struct EventHist_t : public HistBase_t {
    TH1F*    fRunNumber  ;
    TH1F*    fEventNumber;
    TH1F*    fNVDetHits  ;
    TH1F*    fNVDetHits_9;
    TH1F*    fNVDetHits_13;
    TH1F*    fETot_13;
  };

  struct VDetHist_t : public HistBase_t {
    TH1F*    fIndex   ;
    TH1F*    fPDGCode ;		       //
    TH1F*    fGenCode ;		       // generator code
    TH1F*    fMomentum;
    TH1F*    fTime    ;
    TH2F*    fYVsX    ;                // different VD's have different orientation
    TH2F*    fYVsZ    ;                // fill both hist's
    TH1F*    fPt      ;                // transverse mom
    TH1F*    fPp      ;                // momentum component parallel to the solenoid axis
    TH1F*    fTanTh   ;		       // tan (pitch angle)
  };

  struct SimpData_t {
    int           fIndex;		// so far, not used
    TSimParticle* fParent;              // muon parent
  };
//-----------------------------------------------------------------------------
  enum { kNEventHistSets =  100 };
  enum { kNSimpHistSets  = 1000 };
  enum { kNVDetHistSets  = 1000 };

  struct Hist_t {
    EventHist_t* fEvent[kNEventHistSets];
    SimpHist_t*  fSimp [kNSimpHistSets ];
    VDetHist_t*  fVDet [kNVDetHistSets ];
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used
  //  TStepPointMCBlock*    fStepPointMCBlock;
  TSimpBlock*           fSimpBlock;  
  TStepPointMCBlock*    fVDetBlock;
					// transient data
  SimpData_t            fSimpData[100];

  TSimParticle*         fMuon;		// pointer to stopped muon (pend=0)
  TSimParticle*         fParent;
  TSimParticle*         fProton;

  int                   fNVDetHits  ;
  int                   fNVDetHits_9;
  int                   fNVDetHits_13;
  float                 fETot_13;
					// histograms filled
  Hist_t                fHist;

  int                   fNVDet;         // max number of the VD used
  VDetData_t            fVDet[200];     // helper , used in FillVDetHistograms

  TString               fVDetBlockName;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TMuonStopAnaModule(const char* name="MuonStopAna", const char* title="MuonStopAna");
  ~TMuonStopAnaModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Hist_t*            GetHist        () { return &fHist;        }
//-----------------------------------------------------------------------------
// setters
//-----------------------------------------------------------------------------
  void SetVDetBlockName(const char* Name) { fVDetBlockName = Name; }
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
  void    BookSimpHistograms  (HistBase_t* Hist, const char* Folder);
  void    BookVDetHistograms  (HistBase_t* Hist, const char* Folder);
  void    BookEventHistograms (HistBase_t* Hist, const char* Folder);

  void    FillEventHistograms (HistBase_t* Hist);
  void    FillSimpHistograms  (HistBase_t* Hist, TSimParticle* Simp, SimpData_t* SimpData);
  void    FillVDetHistograms  (HistBase_t* Hist, TStepPointMC* Step);

  void    BookHistograms();
  void    FillHistograms();

  VDetData_t*  GetVDetData(int VDetID);

  void    Debug();
//-----------------------------------------------------------------------------
// test
//-----------------------------------------------------------------------------
  void    Test001();

  ClassDef(TMuonStopAnaModule,0)
};

#endif
