///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef murat_ana_TG4ValidationModule_hh
#define murat_ana_TG4ValidationModule_hh

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

#include "Stntuple/loop/TStnModule.hh"
#include "Stntuple/obj/TSimpBlock.hh"

#include "murat/ana/HistBase_t.h"
#include "murat/ana/AnaDefs.hh"

class TG4ValidationModule: public TStnModule {
public:
//-----------------------------------------------------------------------------
//  histograms
//-----------------------------------------------------------------------------
  struct SimpHist_t : public HistBase_t {
    TH1F*      fPDGCode[2];  // just different ranges
    TH1F*      fMomentum[2];
    TH1F*      fCosTheta;
    TH1F*      fTheta;
  };

  struct EventHist_t : public HistBase_t {
    TH1F*      fRunNumber;
    TH1F*      fEventNumber;
    TH1F*      fNPiMinus;
    TH1F*      fNPiPlus;
    TH1F*      fNPi0;
    TH1F*      fNPions;
    TH1F*      fNProtons;
    TH1F*      fNNeutrons;
  };


  struct SimpData_t {
    TSimParticle*  fParticle;
    TSimParticle*  fParent;
    TSimParticle*  fGParent;
  };

//-----------------------------------------------------------------------------
  enum { kNEventHistSets =  100 };
  enum { kNSimpHistSets  = 1000 };

  struct Hist_t {
    EventHist_t*  fEvent[kNEventHistSets];
    SimpHist_t*   fSimp [kNSimpHistSets];
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used
  TSimpBlock*           fSimpBlock;  
					// transient data
  SimpData_t            fSimpData[1000];

  TSimParticle*         fProton;
					// histograms filled
  Hist_t                fHist;

  int                   fNPiMinus;
  int                   fNPiPlus;
  int                   fNPi0;
  int                   fNNeutrons;
  int                   fNProtons;
  int                   fNPions;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TG4ValidationModule(const char* name="G4Validation", const char* title="G4Validation");
  ~TG4ValidationModule();
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
  void    BookSimpHistograms         (HistBase_t* Hist, const char* Folder);
  void    BookEventHistograms        (HistBase_t* Hist, const char* Folder);

  void    FillSimpHistograms         (HistBase_t* Hist, TSimParticle* Simp, SimpData_t* SimpData);
  void    FillEventHistograms        (HistBase_t* Hist);

  void    BookHistograms();
  void    FillHistograms();


  void    Debug();
//-----------------------------------------------------------------------------
// test
//-----------------------------------------------------------------------------
  void    Test001();

  ClassDef(TG4ValidationModule,0)
};

#endif
