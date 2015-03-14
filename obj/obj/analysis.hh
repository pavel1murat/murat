//-----------------------------------------------------------------------------
// code for high-level analysis, making tables for the publication etc
// assume that each process has a histogram file 'HistFile' with all 
// the needed histograms defined in there
// "signal" is the last background process
// analysis nas a name (i.e. "zee") and Title (i.e. "ze1sxd"), which is the 
// the same - ds_id 
//-----------------------------------------------------------------------------
#ifndef murat_ana_analysis_hh
#define murat_ana_analysis_hh

#include "TNamed.h"
#include "TClonesArray.h"

#include "aprocess.hh"
#include "murat/obj/TAnalysisDataset.hh"

class THStack;
class TLegend;

//-----------------------------------------------------------------------------
class analysis: public TNamed {
public:

  enum { kDebugFlags = 100 } ;

  struct process_dat_t {
    TString name;
    int     mcflag;
  };

  struct run_range_dat_t {
    aprocess*       fDat;	          // ! data
    aprocess*       fSig;                 // ! signal MC, one of the "backgrounds"
    TClonesArray*   fBgr;		  // ! list of background processes
    int             fNBgr;		  // ! number of background processes 
    double*         fQSig;                // ! 
    double*         fQBgr;                // ! 
    double*         fESig;                // ! 
    double*         fEBgr;                // ! 
    double*         fChi2;                // ! 
    double*         fQMcSig;              // ! # of signal MC events for different bins
    double          fIntLumi;		  // ! integrated luminosity
    		    
    TObjArray*      fListOfHistograms;
    TObjArray*      fListOfDmbHistograms;

    aprocess*  GetData    () { return fDat;     }
    aprocess*  GetSignalMC() { return fSig;     }
    double     GetIntLumi () { return fIntLumi; }
  };

  TString           fDataset;             // ! name of the analysis dataset
  TString           fVersion;             // ! may have multiple versions with different
				          //   histogram files
  TAnalysisDataset* fDsMetadata;	  // ! information about the datasets
  TString           fXSecModule;
  int               fNRunRanges;
  run_range_dat_t*  fRR;
  int               fLumiBin;             // ! defines which GRL is used

  int               fDebugFlag[kDebugFlags];  // !

  TH1*              fDatHist;
  TH1*              fBgrHist;
  THStack*          fBgrStack;
  TLegend*          fLegend;
					// legend position
  double            fLegendXMin;
  double            fLegendYMin;
  double            fLegendXMax;
  double            fLegendYMax;
//-----------------------------------------------------------------------------
// methods
//-----------------------------------------------------------------------------
  analysis(const char* Name, const char* Title, const char* Version, int LumiBin);
  virtual ~analysis();

  int              init();
  int              GetLumiBin()              { return fLumiBin            ; }
  int              NRunRanges ()             { return fNRunRanges         ; }
  run_range_dat_t* GetRunRange(int RunRange) { return (fRR+RunRange)      ; }
  aprocess*        GetData    (int RunRange) { return (fRR+RunRange)->fDat; }
  aprocess*        GetSignalMC(int RunRange) { return (fRR+RunRange)->fSig; }

  TH1*             GetDatHist () { return fDatHist; }
  TH1*             GetBgrHist () { return fBgrHist; }
  THStack*         GetBgrStack() { return fBgrStack;}

  aprocess*        GetProcess(const char* Name, int RunRange) const ;
  aprocess*        GetProcess(int  I, int RunRange) const ;

  int GetNProcesses(int RunRange) { return (fRR+RunRange)->fBgr->GetEntriesFast(); }

  TAnalysisDataset* GetDsMetadata()          { return fDsMetadata; }
//-----------------------------------------------------------------------------
// produce nice plots, 'nice' is the KEY word here
// if McFlag != 0 then plotted are only contributions of the background processes 
// with flag <= McFlag 
//-----------------------------------------------------------------------------
  void            plot        (const char* Module    , 
			       const char* HistSet   , 
			       int         RunRange  ,
			       int         Bin       , 
			       const char* HistName  , 
			       int         Rebin=1   ,
			       int         McFlag=0  ,
			       double      XMin=-999.,
			       double      XMax= 999.,
			       const char* AxisLabel=0,
			       int         Lumi=0    ,
			       int         YMax=0    ,
			       const char* YUnits=0,
                               float       AxisLabelSize=0,
			       int         Ndivisions=0,
			       int         CombineProcesses=0);

  void            plot_no_data(const char* Module    , 
			       const char* HistSet   , 
			       int         RunRange  ,
			       int         Bin       , 
			       const char* HistName  , 
			       int         Rebin=1   ,
			       int         McFlag=0  ,
			       double      XMin=-999.,
			       double      XMax= 999.,
			       const char* AxisLabel=0,
			       int         Lumi=0    );

  void            plot_mc_only(const char* Module  , 
			       const char* HistSet , 
			       int         RunRange,
			       int         Bin     , 
			       const char* HistName, 
			       int         Rebin=1 );
  
  void           plot_dmboverb(const char* Module  , 
			       const char* HistSet , 
			       int         RunRange,
			       int         Bin     ,
			       const char* HistName,
			       int         Rebin   ,
			       int         McFlag  ,
			       double      XMin    ,
			       double      XMax    ,
			       const char* Axislabel,
			       int         Lumi    ,
			       int         SubtractSignal=0);

//-----------------------------------------------------------------------------
// get 'dmb' (data-minus-background) histogram
// can either subtract all the MC contributions, or leave the signal in
//-----------------------------------------------------------------------------
  int            get_dmb (const char* Module        ,
			  const char* HistSet       ,
			  int         RunRange      ,
			  int         Bin           ,
			  const char* HistName      ,
			  int         SubtractSignal, 
			  int         Rebin         ,
			  TH1*&       Hist          );
//-----------------------------------------------------------------------------
// background model 10:
// get 'dmb_1' (data-minus-background) histogram, all the backgrounds except 
// the last, signal, one are subtracted, 
// signal is scaled to the data in the region [XMin,XMax] and also subtracted
//-----------------------------------------------------------------------------
  int            get_dmb_1 (const char* Module        ,
			    const char* HistSet       ,
			    int         RunRange      ,
			    int         Bin           ,
			    const char* HistName      ,
			    double      XMin          ,
			    double      XMax          ,
			    int         Rebin         ,
			    TH1*&       Hist          );
//-----------------------------------------------------------------------------
// get background histogram: sum of all the MC backgrounds minus MC signal...
//-----------------------------------------------------------------------------
  int            GetBackgroundHistogram (const char* Module        ,
					 const char* HistSet       ,
					 int         RunRange      ,
					 int         Bin           ,
					 const char* HistName      ,
					 int         Rebin         ,
					 TH1*&       Hist          );
//-----------------------------------------------------------------------------
// name of a histogram in a file... may be overloaded because of varying naming 
// conventions
//-----------------------------------------------------------------------------
  virtual void   GetHistogramName(const char* HistSet     ,
				  int         RunRange    ,
				  int         Bin         ,
				  const char* FitHistName ,
				  char*       FullHistName) = 0 ;
//-----------------------------------------------------------------------------
// provide for getting histograms produced for one of the processes 
// by a different module (Aidan's histograms for fakes come from 'ZAna'...
//-----------------------------------------------------------------------------
  virtual int   GetHistogram(aprocess*       Process ,
			     const char*     Module  ,
			     const char*     HistSet ,
			     int             RunRange,
			     int             Bin     ,
			     const char*     HistName,
			     int             Rebin   ,
			     TH1*&           Hist    ) = 0;
//-----------------------------------------------------------------------------
// name of a rebinned histogram to be 
// name of a histogram itself doesn't depend on the run range, however it may
// be overloaded!
//-----------------------------------------------------------------------------
  virtual void   get_hist_name(const char* Module        ,
			       const char* HistSet       ,
			       int         RunRange      ,
			       int         Bin           ,
			       const char* FitHistName   ,
			       int         Rebin         ,
			       char*       Name          ) = 0 ;

  TClonesArray* get_list_of_backgrounds(int RunRange) { return (fRR+RunRange)->fBgr; }

  aprocess*      get_background(int RunRange, const char* Name) {
    return (aprocess*) (fRR+RunRange)->fBgr->FindObject(Name);
  }

  aprocess*      get_background(int RunRange, int I) {
    return (aprocess*) (fRR+RunRange)->fBgr->UncheckedAt(I);
  }
//-----------------------------------------------------------------------------
// return normalization of the background histogram
// folder name is made out of HistName and Bin as follows: 
//
//                       sprintf(folder,"%s_%i",HistName,bin);
//
// so HistName coulde be "ele","tau","event" etc...
// 'aprocess' defines the RunRange...
//-----------------------------------------------------------------------------
  virtual double GetFitQEvents(const aprocess* Process, 
			       const char*     Module ,
			       const char*     HistSet,
			       int             RunRange,
			       int             Bin    ,
			       const char*     HistName) const ;

  virtual double GetIDEff     (const aprocess* Process, 
			       const char*     Module ,
			       const char*     HistSet,
			       int             RunRange,
			       int             Bin    ,
			       const char*     HistName) const ;
//-----------------------------------------------------------------------------
// 'HistSet' and 'Bin' account for tight/loose and also for 
// one lepton/several leptons, so returned is the ID efficiency SF per event
//-----------------------------------------------------------------------------
  virtual double GetIDEffSF   (const aprocess* Process, 
			       const char*     Module ,
			       const char*     HistSet,
			       int             RunRange,
			       int             Bin    ,
			       const char*     HistName) const ;

  virtual double GetTrEff     (const aprocess* Process, 
			       const char*     Module ,
			       const char*     HistSet,
			       int             RunRange,
			       int             Bin    ,
			       const char*     HistName) const ;

  virtual double GetTrEffSF   (const aprocess* Process, 
			       const char*     Module ,
			       const char*     HistSet,
			       int             RunRange,
			       int             Bin    ,
			       const char*     HistName) const ;

  virtual double GetAcceptance(const aprocess* Process, 
			       const char*     Module ,
			       const char*     HistSet,
			       int             RunRange,
			       int             Bin    ,
			       const char*     HistName) const ;

  virtual double GetCrossSection(const char*   Module  ,
				 const char*   HistSet ,
				 int           RunRange, 	
				 int           Bin     ,
				 const char*   HistName) const ;

  virtual void    Print(const char* Module  ,
			const char* HistSet ,
			int         Bin     ,
			const char* HistName,
			const char* Opt     ) const ;
//-----------------------------------------------------------------------------
// setters
//-----------------------------------------------------------------------------
  void SetLegendPosition(double XMin, double YMin, double XMax, double YMax) {
    fLegendXMin = XMin;
    fLegendYMin = YMin;
    fLegendXMax = XMax;
    fLegendYMax = YMax;
  }
//-----------------------------------------------------------------------------
// overloaded functions of TObject
//-----------------------------------------------------------------------------
  virtual void    Print(const char* Opt = "") const ;

  ClassDef(analysis,0)

};

#endif
