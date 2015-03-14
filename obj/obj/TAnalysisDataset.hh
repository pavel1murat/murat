///////////////////////////////////////////////////////////////////////////////
// information related to ELECTRON_CENTRAL_18 trigger path
// and its time evolution
///////////////////////////////////////////////////////////////////////////////
#ifndef murat_obj_TAnalysisDataset_hh
#define murat_obj_TAnalysisDataset_hh

#include "TNamed.h"
#include "TObjArray.h"

class TDsMetadata;

class TAnalysisDataset: public TNamed {
public:

  struct run_range_t {
    int         fMinRun    ;
    int         fMaxRun    ;
    const char* fDsID      ;              // "0d" or ... "0m"
    int         fL3PathTag ;
    double      fIntLumi[4];    // [0]:e [1]:mu [2]:emu [3]:emusi
  };

				        // data for a run range
  struct dataset_t {
    int          fMinRun;
    int          fMaxRun;
    const char*  fDsID;
    const char*  fGoodRunList;
    int          fMcFlag;
    int          fNTotalMc[4]; // [0]:e [1]:mu [2]:emu [3]:emusi

    const char*  GetDsID  () const { return fDsID  ; };
    int          GetMcFlag() const { return fMcFlag; };
  };

  int          fNRunRanges;
  run_range_t* fRunRange;            // !
  run_range_t* fJPRunRange;          // ! Joint Physics run ranges (default)
  int          fLumiBin;             // 0:e-only, 1:mu-only, 2:emu-nosi,3:emusi

  TObjArray*   fListOfDatasets;      // ! list of TDsMetadata 

  int          fPlotMode;            // 1:note, 2:talk, 3:paper
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
public:
  TAnalysisDataset();

	                             // 'LumiBin': 0='eonly', 1:'muonly', 2='emu', 3='emusi'
 
  TAnalysisDataset(const char* Name, const char* Title, int LumiBin);
  ~TAnalysisDataset();

  const char*  Version        ()             { return GetTitle(); }
  int          NRunRanges     ()             { return fNRunRanges ; }
  run_range_t* GetRunRange    (int i       ) { return fRunRange+i ; }
  int          L3PathTag      (int RunRange) { return fRunRange[RunRange].fL3PathTag; }
  int          MinRun         (int RunRange) { return fRunRange[RunRange].fMinRun   ; }
  int          MaxRun         (int RunRange) { return fRunRange[RunRange].fMaxRun   ; }

				// LumiBin: 0='eonly', 1:'muonly', 2='emu', 3='emusi'

  virtual double       Luminosity     (int RunRange, int LumiBin);

  int          UsedRR         (int RunRange) { 
    return (fRunRange[RunRange].fMinRun > 0);
  }

  run_range_t* GetJPRunRange  (int RunRange) { return fJPRunRange+RunRange ; }
//-----------------------------------------------------------------------------
// setters
//-----------------------------------------------------------------------------
  void SetPlotMode(int Mode) { fPlotMode = Mode; }
//-----------------------------------------------------------------------------
// virtual functions
//-----------------------------------------------------------------------------
  virtual int  FindRunRange   (int RunNumber);
//-----------------------------------------------------------------------------
// different MC exist for different run ranges. 'Process'= 'zee' etc ... 
// translation "zee" -> 'ze0s01' happens inside
//-----------------------------------------------------------------------------
  virtual int  SubstituteRunRange(const char* Process , 
				  int         RunRange);

  TDsMetadata* GetDsMetadata(const char* Process) {
    return (TDsMetadata*) fListOfDatasets->FindObject(Process);
  };

  virtual const dataset_t*   GetDataset(const char* Process, int RunRange);

				// type = "data"/"mc"

  virtual int          GetDsID      (const char* Process , 
				     int         McFlag  ,
				     int         RunRange, 
				     char*       DsID    );

  virtual int          GetHistColor   (const char* Process, int McFlag) const ;
  virtual int          GetHistStyle   (const char* Process, int McFlag) const ;

  virtual int          GetHistFileName(const char* Process , 
				       int         McFlag  ,
				       int         RunRange, 
				       char*       HistFileName);

  virtual double       GetTightTauIDSF (const char* Process, int McFlag, int RunRange);
  virtual double       GetLooseTauIDSF (const char* Process, int McFlag, int RunRange);
  virtual double       GetTightEleIDSF (const char* Process, int McFlag, int RunRange);
  virtual double       GetLooseEleIDSF (const char* Process, int McFlag, int RunRange);
  virtual double       GetTriggerEffSF (const char* Process, int McFlag, int RunRange);
  virtual double       GetTrigEffErr   (const char* Process, int McFlag, int RunRange);
  virtual int          GetNTotalMc     (const char* Process, int McFlag, int RunRange);
  virtual double       GetXSec         (const char* Process, int McFlag);

//-----------------------------------------------------------------------------
// correction to the Z-->ll MC acceptance comes from the cutoff on M(ll)
// which results in more events, than defined by the Z-->ee cross section 
// being generated per unit of integrated luminosity
// in principle, can depend on the run range
//-----------------------------------------------------------------------------
  virtual double       GetAccCorrFactor(const char* Process, int McFlag, int RunRange);
  
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  virtual void         Print(const char* Opt="") const;

  ClassDef(TAnalysisDataset,0)
};

//-----------------------------------------------------------------------------
// a dataset needs to have a name!
//-----------------------------------------------------------------------------
class TDsMetadata: public TNamed {
public:
  const TAnalysisDataset::dataset_t* fData;

  TString fLabel;             // label to be printed on a histogram

  TDsMetadata(const char* Name, 
	      const char* Title, 
	      const TAnalysisDataset::dataset_t* Data, 
	      const char* Label=0);

				// Dataset doesnt' own anything
  ~TDsMetadata();

  const char* GetLabel() { return fLabel.Data(); }

  ClassDef(TDsMetadata,0)
};


#endif
