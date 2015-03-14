//-----------------------------------------------------------------------------
// code for high-level analysis, making tables for the publication etc
// assume that each process has a histogram file 'HistFile' with all 
// the needed histograms defined in there
// represents data for a given run range
//-----------------------------------------------------------------------------
#ifndef murat_ana_aprocess_hh
#define murat_ana_aprocess_hh

#include "TNamed.h"
#include "TObjArray.h"
#include "TFile.h"
#include "TH1.h"

class TAnalysisDataset;
class analysis;

//-----------------------------------------------------------------------------
class aprocess : public TNamed {
public:
  enum {kNBins       = 20};
  enum {kNDebugFlags = 10};

  TString           fHistFileName;      // ! histogram file for this process
  Int_t             fMcFlag;		// ! 0=data, 1=MC BGR, 3=MC SIGNAL
  analysis*         fAnalysis;	        // ! backward reference to analysis channel
  TAnalysisDataset* fDsMetadata;        // ! dataset data
  Int_t             fRunRange;          // ! run range
  Int_t             fSrr;               // ! substitute run range...
  double            fAcc;		// ! 
  double            fTrEff;		// ! 
  double            fNTotal;            // ! total number of events (MC), for generated XSection
				        // ! defined only for MC datasets
//-----------------------------------------------------------------------------
// for the data (i.e. QCD) fXSec = -1, this means that fQEvents needs to be calculated 
// differently
//-----------------------------------------------------------------------------
  double     fXSec;			// ! cross section, nb
  double     fIntLumi;			// ! 
//   TObjArray* fListOfHistograms;         // ! for the data may need to 'assign' a histogram
//                                         //   to a background (fit result, for example)

  double     fTightIdSF;                // ! scale factor (tight ID efficiency etc)
  double     fLooseIdSF;                // ! scale factor (loose ID efficiency etc)
  double     fTrEffSF;		        // ! scale factor (trigger efficiency)
  double     fQEvents;                  // ! expected number of events in the histogram
  int        fNormalize;                // ! 1:do normalization, 0: don't do it
  TH1*       fNormHist;			// ! normalization histogram
  int        fTight[kNBins];            // ! tight-loose flags
  double     fFitQEvents[kNBins];       // ! normalization - per bin... assume only one fit at a time
				        // ! in principle, could ask analysis about it....
  int        fColor;                    // ! color to draw
  int        fStyle;                    // ! style to draw
  int        fDebug[100];               // ! different debug flags
//-----------------------------------------------------------------------------
// methods
//-----------------------------------------------------------------------------
  aprocess(const char* Name = "undefined");

  aprocess(const char*       Name    , 
	   int               McFlag  ,
	   int               RunRange,
	   analysis*         Analysis);

// 	   int         Srr,
// 	   const char* HistFileName, 
// 	   int         NTotal, 
// 	   double      XSec, 
// 	   double      IntLumi,
// 	   TH1*        NormHist);

  virtual ~aprocess();
//-----------------------------------------------------------------------------
// if Tight=1 then signal contribution will be scaled by tight ID efficiency,
// otherwise - by loose ID efficiency, tight/loose is defined by fTight[Bin]
// returns CLONE !!!
//-----------------------------------------------------------------------------
  int    get_h1(const char* Module  , 
		const char* HistSet , 
		int         Bin     , 
		const char* HistName, 
		int         Rebin   ,
		TH1*&       Hist    ) ;

  int    get_h1_unorm(const char* Module  , 
		      const char* HistSet , 
		      int         Bin     , 
		      const char* HistName, 
		      int         Rebin   ,
		      TH1*&       Hist    ,
		      double*     Anorm   ) ;

  const char*  GetDsID        () const { return GetTitle(); }
  int          GetMcFlag      () const { return fMcFlag; }

  const char*  GetHistFileName() const { return fHistFileName.Data() ; }

  double       GetIDEffSF(const char*  Module ,
			  const char*  HistSet,
			  int          Bin    ,
			  const char*  HistName) const ;

  double       GetTrEff  (const char*  Module ,
			  const char*  HistSet,
			  int          Bin    ,
			  const char*  HistName) const ;

  double       GetTrEffSF(const char*  Module ,
			  const char*  HistSet,
			  int          Bin    ,
			  const char*  HistName) const ;

  double       GetAcceptance(const char*  Module ,
			     const char*  HistSet,
			     int          Bin    ,
			     const char*  HistName) const ;

  double       GetFitQEvents(const char*  Module ,
			     const char*  HistSet,
			     int          Bin    ,
			     const char*  HistName) const ; 

  double       GetNTotal() const { return fNTotal; }
  int          GetStyle () const ; // { return fStyle ; }
  int          GetColor () const ; // { return fColor ; }
  const char*  GetLabel () const; // to be plotted on a TLegend...
//-----------------------------------------------------------------------------
// for Z->ll MC acceptance correction factor arises because of the M(ll) cut-off
// which is ~1.4 for M(ll) > 30 and ~1.96 for M(ll)>20
//-----------------------------------------------------------------------------
  double       GetAccCorrFactor() const; 
//-----------------------------------------------------------------------------
// this can be overloaded to deal with the missing MC datasets
// 'Option' = 'file' or 'job'
// for run range-dependent histograms use 'HistSet'="zee_rr_00" etc...
//-----------------------------------------------------------------------------
  virtual void GetHistogramName(const char* Module,
				const char* HistSet,
				int         Bin,
				const char* HistName,
				const char* Option,
				char*       FileHistName) const;
//-----------------------------------------------------------------------------
// setters
//-----------------------------------------------------------------------------
  void   SetHistFileName(const char* HistFileName) { fHistFileName = HistFileName; }
  void   SetNTotal      (int         N           ) { fNTotal       = N           ; }
  void   SetXSec        (double      XSec        ) { fXSec         = XSec        ; }
  void   SetIntLumi     (double      IntLumi     ) { fIntLumi      = IntLumi     ; }
  void   SetNormalize   (int         YesNo       ) { fNormalize    = YesNo       ; }
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  virtual void    Print(const char* Opt = "") const;

  virtual void    Print(const char* Module  ,
			const char* HistSet ,
			int         Bin     ,
			const char* HistName,
			const char* Opt     ) const;

  ClassDef(aprocess,0)

};


#endif
