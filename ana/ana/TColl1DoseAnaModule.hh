///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef murat_ana_TColl1DoseAnaModule_hh
#define murat_ana_TColl1DoseAnaModule_hh

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

class TColl1DoseAnaModule: public TStnModule {
public:
  enum {kMaxNSteps = 1000};
  enum {kNSlices   =   10};

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
    int            fBinX;               // depth bin
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

    TH2F*      fYVsX;
    TH2F*      fYVsXWtE;
    TH2F*      fYVsXDose;
  };

  struct EventHist_t : public HistBase_t {
    TH1F*      fRunNumber;
    TH1F*      fEventNumber;
    TH1F*      fNPbrAbsSteps;
  };

//-----------------------------------------------------------------------------
  enum { kNEventHistSets = 100  };
  enum { kNSpmcHistSets  = 1000 };

  struct Hist_t {
    EventHist_t* fEvent      [kNEventHistSets];
    SpmcHist_t*  fPbarAbsSpmc[kNSpmcHistSets ];
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used

  TStepPointMCBlock*    fPbarAbsSpmcBlock;

					// histograms filled
  Hist_t                fHist;

  int                   fNPbrAbsSteps;

  int                   fPdgCode;
  int                   fGeneratorCode;

  float                 fNPOT;
  int                   fNPerPOT;

  SpmcData_t            fPbarAbsSpmcData[1000];
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TColl1DoseAnaModule(const char* name="Coll1DoseAna", const char* title="Coll1DoseAna");
  ~TColl1DoseAnaModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Hist_t*            GetHist        () { return &fHist;        }
//-----------------------------------------------------------------------------
// setters
//-----------------------------------------------------------------------------
  void    SetPdgCode      (int Code) { fPdgCode       = Code; }
  void    SetGeneratorCode(int Code) { fGeneratorCode = Code; }
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
  void    BookPbrAbsSpmcHistograms(HistBase_t* Hist, const char* Folder);

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

  ClassDef(TColl1DoseAnaModule,0)
};

#endif
