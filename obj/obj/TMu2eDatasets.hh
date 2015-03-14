///////////////////////////////////////////////////////////////////////////////
// book-keeping information related to W->enu strip from CPH1?? dataset
///////////////////////////////////////////////////////////////////////////////
#ifndef murat_ana_TMu2eDatasets_hh
#define murat_ana_TMu2eDatasets_hh

#include "murat/obj/TAnalysisDataset.hh"

class TMu2eDatasets: public TAnalysisDataset {
public:
  TDsMetadata*   fMu2e;                     // ! MC
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  TMu2eDatasets(const char* Name, const char* Title, int LumiBin);
  ~TMu2eDatasets();
//-----------------------------------------------------------------------------
// setters
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// other
//-----------------------------------------------------------------------------
  double  GetCentralEtSF(const char* Process, int McFlag, int RunRange);
  void    Help(int Mode);
//-----------------------------------------------------------------------------
// overloaded virtual functions of TAnalysisDataset
//-----------------------------------------------------------------------------
				        // Type = "data"/"mc"

  virtual int          GetHistColor   (const char* Process, int McFlag) const ;
  virtual int          GetHistStyle   (const char* Process, int McFlag) const ;

  virtual int          GetHistFileName(const char* Process, 
				       int         McFlag, 
				       int         RunRange, 
				       char*       HistFileName);

  virtual double  GetTightEleIDSF(const char* Process, int McFlag, int RunRange);
  virtual double  GetLooseEleIDSF(const char* Process, int McFlag, int RunRange);

  virtual double  GetTriggerEffSF (const char* Process, int McFlag, int RunRange);
  virtual double  GetTrigEffErr   (const char* Process, int McFlag, int RunRange);

  virtual int     SubstituteRunRange(const char* Process , 
				     int         RunRange);

  virtual double  GetAccCorrFactor(const char* Process, int McFlag, int RunRange);

  ClassDef(TMu2eDatasets,0)
};


#endif
