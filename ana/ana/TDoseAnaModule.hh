///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef murat_ana_TDoseAnaModule_hh
#define murat_ana_TDoseAnaModule_hh

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

class TDoseAnaModule: public TStnModule {
public:

  struct SpmcData_t {
    TSimParticle*  fParticle;
    TSimParticle*  fParent;
    TSimParticle*  fGParent;
    float          fGpTheta;
    float          fP;
    float          fCosTh;
    float          fTime;
    float          fX;			// local horizontal coord (X or Z)
    float          fY;			// local vertical   coord
  };

//-----------------------------------------------------------------------------
//  histograms
//-----------------------------------------------------------------------------
  struct SpmcHist_t : public HistBase_t {
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
    TH1F*      fCosTh;

    TH2F*      fYVsZ;
    TH2F*      fYVsX;
  };

  struct EventHist_t : public HistBase_t {
    TH1F*      fRunNumber;
    TH1F*      fEventNumber;
    TH1F*      fNSteps;
  };

//-----------------------------------------------------------------------------
  enum { kNEventHistSets = 100  };
  enum { kNSpmcHistSets  = 1000 };

  struct Hist_t {
    EventHist_t* fEvent[kNEventHistSets];
    SpmcHist_t*  fSpmc [kNSpmcHistSets ];
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used
  TStepPointMCBlock*    fSpmcBlock;
  TString               fSpmcBlockName;
					// histograms filled
  Hist_t                fHist;

  int                   fNSteps;
  int                   fPdgCode;
  int                   fGeneratorCode;

  SpmcData_t            fSpmcData[1000];
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TDoseAnaModule(const char* name="DoseAna", const char* title="DoseAna");
  ~TDoseAnaModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Hist_t*            GetHist        () { return &fHist;        }
//-----------------------------------------------------------------------------
// setters
//-----------------------------------------------------------------------------
  void    SetPdgCode      (int Code) { fPdgCode       = Code; }
  void    SetGeneratorCode(int Code) { fGeneratorCode = Code; }

  void    SetSpmcBlockName(const char* Name) { fSpmcBlockName = Name; }
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
  void    BookSpmcHistograms (HistBase_t* Hist, const char* Folder);
  void    BookEventHistograms(HistBase_t* Hist, const char* Folder);

  void    FillSpmcHistograms (HistBase_t* Hist, TStepPointMC* Step, SpmcData_t*  Sd);
  void    FillEventHistograms(HistBase_t* Hist);

  void    BookHistograms();
  void    FillHistograms();

  void    Debug();
//-----------------------------------------------------------------------------
// test
//-----------------------------------------------------------------------------
  void    Test001();

  ClassDef(TDoseAnaModule,0)
};

#endif
