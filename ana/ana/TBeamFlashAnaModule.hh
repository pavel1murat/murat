///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef murat_ana_TBeamFlashAnaModule_hh
#define murat_ana_TBeamFlashAnaModule_hh

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

#include "Stntuple/loop/TStnModule.hh"

#include "Stntuple/obj/TGenpBlock.hh"
#include "Stntuple/obj/TSimpBlock.hh"
#include "Stntuple/obj/TStepPointMCBlock.hh"

#include "Stntuple/base/TStnArrayI.hh"

#include "Stntuple/alg/TStnTrackID.hh"
#include "Stntuple/alg/TEmuLogLH.hh"

#include "murat/ana/HistBase_t.h"

#include "murat/ana/AnaDefs.hh"

class TBeamFlashAnaModule: public TStnModule {
public:
//-----------------------------------------------------------------------------
//  histograms
//-----------------------------------------------------------------------------
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
    TH1F*      fTime;
    TH1F*      fStepLength;

    TH1F*      fMomentum;

    TH2F*      fYVsZ;
    TH2F*      fYVsX;

    TH1F*      fGpPDGCode[2];
    TH2F*      fGpCosThVsMom;
  };

  struct EventHist_t : public HistBase_t {
    TH1F*      fRunNumber;
    TH1F*      fEventNumber;
  };


  struct SpmcData_t {
    TSimParticle*  fParticle;
    TSimParticle*  fParent;
    TSimParticle*  fGParent;
  };

//-----------------------------------------------------------------------------
  enum { kNEventHistSets        =  100 };
  enum { kNStepPointMCHistSets  = 1000 };

  struct Hist_t {
    EventHist_t*        fEvent       [kNEventHistSets];
    StepPointMCHist_t*  fStepPointMC [kNStepPointMCHistSets];
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used
  TStepPointMCBlock*    fStepPointMCBlock;
  TSimpBlock*           fSimpBlock;  
					// transient data
  SpmcData_t            fSpmcData[100];

  TSimParticle*         fProton;
					// histograms filled
  Hist_t                fHist;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TBeamFlashAnaModule(const char* name="BeamFlashAna", const char* title="BeamFlashAna");
  ~TBeamFlashAnaModule();
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
  void    BookStepPointMCHistograms  (HistBase_t* Hist, const char* Folder);
  void    BookEventHistograms        (HistBase_t* Hist, const char* Folder);

  void    FillStepPointMCHistograms  (HistBase_t* Hist, TStepPointMC* Step, SpmcData_t* Sd);
  void    FillEventHistograms        (HistBase_t* Hist);

  void    BookHistograms();
  void    FillHistograms();


  void    Debug();
//-----------------------------------------------------------------------------
// test
//-----------------------------------------------------------------------------
  void    Test001();

  ClassDef(TBeamFlashAnaModule,0)
};

#endif
