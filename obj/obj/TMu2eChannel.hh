//-----------------------------------------------------------------------------
// code for high-level analysis, making tables for the publication etc
// assume that each process has a histogram file 'HistFile' with all 
// the needed histograms defined in there
// name: tau_1 or tau_7
//-----------------------------------------------------------------------------
#ifndef murat_limits_TMu2eChannel_hh
#define murat_limits_TMu2eChannel_hh

#include "TNamed.h"

#include "murat/obj/analysis.hh"
//-----------------------------------------------------------------------------
class TMu2eChannel: public analysis {
public:
//-----------------------------------------------------------------------------
// variables used for 2-component fitting
//-----------------------------------------------------------------------------
  static TH1*       g_h1;		// ! signal-subtracted data distribution
  static TH1*       g_h2;		// ! MC W->tau nu distribution
  static TH1*       g_h3;		// ! QCD distribution

  TF1*              fBgrFitF;      // !

  double            fFitPar[10];   // !
  double            fMass;

  process_dat_t     fProc[10];
//-----------------------------------------------------------------------------
// Period = 1: 0D, 72/pb
//        = 2: 0D, 300/pb
//        = 3: 0H
//        = 4: 0I     (there were 2 different triggers, may need to split)
//-----------------------------------------------------------------------------
  TMu2eChannel(const char* Version, Int_t LumiBin, const char* Signal = "gpt000");
  ~TMu2eChannel();
//-----------------------------------------------------------------------------
// overloaded methods of the 'analysis' class
//-----------------------------------------------------------------------------
  virtual void GetHistogramName(const char* HistSet     ,
				int         RunRange    ,
				int         Bin         ,
				const char* FitHistName ,
				char*       FullHistName);
//-----------------------------------------------------------------------------
// provide for getting histograms produced for one of the processes 
// by a different module (Aidan's histograms for fakes come from 'ZAna'...
//-----------------------------------------------------------------------------
  virtual int    GetHistogram(aprocess*       Process ,
			      const char*     Module  ,
			      const char*     HistSet ,
			      int             RunRange,
			      int             Bin     ,
			      const char*     HistName,
			      int             Rebin   ,
                              TH1*&           Hist    );
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
			       char*       Name          );
//-----------------------------------------------------------------------------
// overloaded functions of 'analysis' class
//-----------------------------------------------------------------------------
  virtual double GetIDEffSF   (const aprocess*    Process, 
			       const char*        Module ,
			       const char*        HistSet,	
			       int                RunRange,
			       int                Bin    ,
			       const char*        HistName) const ;

  virtual double GetAcceptance(const aprocess*    Process, 
			       const char*        Module ,
			       const char*        HistSet,
			       int                RunRange,
			       int                Bin    ,
			       const char*        HistName) const ;

  virtual double GetFitQEvents(const aprocess*    Process, 
			       const char*        Module ,
			       const char*        HistSet,
			       int                RunRange,
			       int                Bin    ,
			       const char*        HistName) const ;

  virtual double GetTrEffSF    (const aprocess*    Process, 
				const char*        Module ,
				const char*        HistSet,
				int                RunRange,
				int                Bin     ,
				const char*        HistName) const;
//-----------------------------------------------------------------------------
// overloaded functions of TObject
//-----------------------------------------------------------------------------

  ClassDef(TMu2eChannel,0)

};

#endif
