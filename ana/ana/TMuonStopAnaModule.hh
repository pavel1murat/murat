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

    TH2F*    fYVsX;
    TH2F*    fYVsX_2480;
    TH2F*    fYVsX_2513;

  };

  struct EventHist_t : public HistBase_t {
    TH1F*    fRunNumber  ;
    TH1F*    fEventNumber;
    TH1F*    fNVdetHits  ;
    TH1F*    fNVdetHits_9;
    TH1F*    fNVdetHits_13;
    TH1F*    fETot_13;
  };

  struct VdetHist_t : public HistBase_t {
    TH1F*    fIndex   ;
    TH1F*    fPDGCode ;		       //
    TH1F*    fGenCode ;		       // generator code
    TH1F*    fMomentum;
    TH1F*    fTime    ;
    TH2F*    fYVsX    ;
  };

  struct SimpData_t {
    int           fIndex;		// so far, not used
    TSimParticle* fParent;              // muon parent
  };

//-----------------------------------------------------------------------------
  enum { kNEventHistSets =  100 };
  enum { kNSimpHistSets  = 1000 };
  enum { kNVdetHistSets  = 1000 };

  struct Hist_t {
    EventHist_t* fEvent[kNEventHistSets];
    SimpHist_t*  fSimp [kNSimpHistSets ];
    VdetHist_t*  fVdet [kNVdetHistSets ];
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used
  //  TStepPointMCBlock*    fStepPointMCBlock;
  TSimpBlock*           fSimpBlock;  
  TStepPointMCBlock*    fVdetBlock;
					// transient data
  SimpData_t            fSimpData[100];

  TSimParticle*         fMuon;
  TSimParticle*         fParent;
  TSimParticle*         fProton;

  int                   fNVdetHits  ;
  int                   fNVdetHits_9;
  int                   fNVdetHits_13;
  float                 fETot_13;
					// histograms filled
  Hist_t                fHist;
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
  void    BookVdetHistograms  (HistBase_t* Hist, const char* Folder);
  void    BookEventHistograms (HistBase_t* Hist, const char* Folder);

  void    FillEventHistograms (HistBase_t* Hist);
  void    FillSimpHistograms  (HistBase_t* Hist, TSimParticle* Simp, SimpData_t* SimpData);
  void    FillVdetHistograms  (HistBase_t* Hist, TStepPointMC* Step);

  void    BookHistograms();
  void    FillHistograms();

  void    Debug();
//-----------------------------------------------------------------------------
// test
//-----------------------------------------------------------------------------
  void    Test001();

  ClassDef(TMuonStopAnaModule,0)
};

#endif
