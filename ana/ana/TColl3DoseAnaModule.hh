///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef murat_ana_TColl3DoseAnaModule_hh
#define murat_ana_TColl3DoseAnaModule_hh

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

class TColl3DoseAnaModule: public TStnModule {
public:
  enum {kMaxNSteps = 10000};
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

    TH2F*      fYVsZ;
    TH2F*      fYVsZWtE;

    TH2F*      fYVsZWtESlice[kNSlices];
    TH2F*      fYVsZDose    [kNSlices];

    TH2F*      fYVsX;
    TH2F*      fYVsXWtE;
  };

  struct EventHist_t : public HistBase_t {
    TH1F*      fRunNumber;
    TH1F*      fEventNumber;
    TH1F*      fNColl31Steps;
    TH1F*      fNColl32Steps;
    TH1F*      fNPbrAbsSteps;
  };

//-----------------------------------------------------------------------------
  enum { kNEventHistSets = 100  };
  enum { kNSpmcHistSets  = 1000 };

  struct Hist_t {
    EventHist_t* fEvent     [kNEventHistSets];
    SpmcHist_t*  fColl31Spmc[kNSpmcHistSets ];
    SpmcHist_t*  fColl32Spmc[kNSpmcHistSets ];
    SpmcHist_t*  fPbrAbsSpmc[kNSpmcHistSets ];
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used

  TStepPointMCBlock*    fColl31SpmcBlock;
  TStepPointMCBlock*    fColl32SpmcBlock;
  TStepPointMCBlock*    fPbrAbsSpmcBlock;

					// histograms filled
  Hist_t                fHist;

  int                   fNColl31Steps;
  int                   fNColl32Steps;
  int                   fNPbrAbsSteps;

  int                   fPdgCode;
  int                   fGeneratorCode;

  float                 fNPOT;
  int                   fNPerPOT;

  SpmcData_t            fColl31SpmcData[1000];
  SpmcData_t            fColl32SpmcData[1000];
  SpmcData_t            fPbrAbsSpmcData[1000];
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TColl3DoseAnaModule(const char* name="Coll3DoseAna", const char* title="Coll3DoseAna");
  ~TColl3DoseAnaModule();
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
  void    BookColl31SpmcHistograms(HistBase_t* Hist, const char* Folder);
  void    BookColl32SpmcHistograms(HistBase_t* Hist, const char* Folder);
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

  ClassDef(TColl3DoseAnaModule,0)
};

#endif
