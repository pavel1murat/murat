///////////////////////////////////////////////////////////////////////////////
// book-keeping class for W->taunu analysis
// all datasets are assumed to be declared in the book=SNUEAU
//
// catalogs: see ncdf192:/cdf/ncdf192/datasets/catalog/
// ---------
// wtaunu :     we0s9t:OK        we0sat:OK           OK
// wenu   :     we0sfe:OK        we0sge:OK           OK
// wmunu  :     we0s7m:OK        we0s8m:OK           OK
// zee    :     ze1s6d:OK        ze1dad:OK           crashed
// zmumu  :     ze1s6m:OK        ze1s9m:OK           OK 
// ztautau:     ze0s8t:OK        ze0sat:OK           ..
//
// data on number of MC events for different trigger configurations come from:
// 
// cat $WORK_DIR/tau_met_doc/n2_mc.wtaunu.log | grep L3
//
// the log files have been produced with
//
// root [0] .L murat/ana/scripts/n_mc_events.C
// root [1] n2_mc("zee","TAU_MET") ; > n2_mc.zee.log
//
// root [2] good_run_list("DQM_V13:100",141544,156487)
// ....
// OFFLINE lumi in all  runs =  106.34802 pb^-1
// OFFLINE lumi in good runs =   89.32262 pb^-1
///////////////////////////////////////////////////////////////////////////////
#include "TSystem.h"
#include "obj/TAnalysisDataset.hh"

ClassImp(TAnalysisDataset)
ClassImp(TDsMetadata)

namespace {

  const int kNRunRanges = 1;

  TAnalysisDataset::run_range_t  run_range[kNRunRanges] = { 
//-----------------------------------------------------------------------------
// L3 path tag is a placeholder for different split into run ranges
//
// luminosities below correspond to DQM_V32:(300,010,310,311)
// for P0-P26 they are already corrected by a "Jaco-factor" of 1.019
// 
//    rmin    rmax  DSID  L3 path tag  int_lumi jaco-corr (1.019) 
//                         (settable)  e          mu         emu       emusi
//---------------- --------- ---------------------------------------------------
    1, 1000000, "0a",       1,      1.,           1.,        1.,        1.  //   0 
  };

};


//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
TAnalysisDataset::TAnalysisDataset(): TNamed() {
  //  fNRunRanges     = kNRunRanges;
  fNRunRanges     = 1;
  fRunRange       = run_range  ;
  fJPRunRange     = run_range  ;
  fLumiBin        = 0;
  fListOfDatasets = 0;
  fPlotMode       = 1;
}

//-----------------------------------------------------------------------------
TAnalysisDataset::TAnalysisDataset(const char* Name, const char* Title, int LumiBin): 
  TNamed(Name,Title) 
{
  fNRunRanges     = kNRunRanges;
  fRunRange       = run_range  ;
  fJPRunRange     = run_range  ;
  fListOfDatasets = new TObjArray();
  fLumiBin        = LumiBin;
  fPlotMode       = 1;
}


//-----------------------------------------------------------------------------
TAnalysisDataset::~TAnalysisDataset() {
  if (fListOfDatasets) delete fListOfDatasets;
}

//-----------------------------------------------------------------------------
int TAnalysisDataset::GetDsID(const char* Name    , 
			      int         McFlag  , 
			      int         RunRange, 
			      char*       DsID    ) {
  // data and MC datasets should have different names

  const dataset_t*  ds   = GetDataset(Name,RunRange);

  strcpy(DsID,ds->fDsID);

  return 0;
}

//-----------------------------------------------------------------------------
int TAnalysisDataset::FindRunRange(int RunNumber) {
  // 

  int rr = fNRunRanges;
  for (int i=0; i<fNRunRanges; i++) {
    if ((RunNumber >= fRunRange[i].fMinRun) && 
	(RunNumber <= fRunRange[i].fMaxRun)    ) {
      rr = i;
      break;
    }
  }

  if (rr == fNRunRanges) {
    Error("FindRunRange",Form("run %7i is outside the run range"));
  }

  return rr;
}

//-----------------------------------------------------------------------------
// virtual, to be overloaded
//-----------------------------------------------------------------------------
const TAnalysisDataset::dataset_t* TAnalysisDataset::GetDataset(const char* Process, 
								int         RunRange) 
{
  TDsMetadata*     dm;
  const dataset_t* ds;
  TString          dsid;
  int              mc_flag, srr;

  dm   = (TDsMetadata*) fListOfDatasets->FindObject(Process);
  if (dm) {
    srr = SubstituteRunRange(Process,RunRange);
    ds = dm->fData+srr;
  }
  else    {
    Error("GetDataset",Form("Process %s NOT FOUND",Process));
    ds = 0;
  }
  return ds;
}


//-----------------------------------------------------------------------------
int TAnalysisDataset::GetHistFileName(const char* Process , 
				      int         McFlag  , 
				      int         RunRange, 
				      char*       HistFileName) {
  Error("GetHistFileName","should not be called!");
  return 0;
}

//-----------------------------------------------------------------------------
int TAnalysisDataset::GetHistColor(const char* Process, int McFlag) const {
  Error("GetHistColor","should not be called!");
  return 1;
}

//-----------------------------------------------------------------------------
int TAnalysisDataset::GetHistStyle(const char* Process, int McFlag) const {
  Error("GetHistStyle","should not be called!");
  return 1;
}

//-----------------------------------------------------------------------------
double TAnalysisDataset::GetLooseTauIDSF(const char* Process, 
					 int         McFlag, 
					 int         RunRange) {
  Error("GetLooseTauIDSF","should not be called!");
  return 0.;
}

//-----------------------------------------------------------------------------
double TAnalysisDataset::GetTightTauIDSF(const char* Process , 
					 int         McFlag  , 
					 int         RunRange) {
  Error("GetTightTauIDSF","should not be called!");
  return 0.;
}

//-----------------------------------------------------------------------------
double TAnalysisDataset::GetLooseEleIDSF(const char* Process , 
					 int         McFlag  , 
					 int         RunRange) {
  Error("GetLooseEleIDSF","should not be called!");
  return 0.;
}

//-----------------------------------------------------------------------------
double TAnalysisDataset::GetTightEleIDSF(const char* Process , 
					 int         McFlag    , 
					 int         RunRange) {
  Error("GetTightEleIDSF","should not be called!");
  return 0.;
}

//-----------------------------------------------------------------------------
// TAU_MET trigger efficiency for because of the specifics of the MC generation 
// min. and max. run numbers 
// may be defined differently for the same run range for different datasets
//-----------------------------------------------------------------------------
double TAnalysisDataset::GetTriggerEffSF(const char* Process , 
					 int         McFlag    , 
					 int         RunRange) {
  Error("GetTriggerEffSF","should not be called!");
  return 0.;
}


//-----------------------------------------------------------------------------
// uncertainty on the trigger efficiency
//-----------------------------------------------------------------------------
double TAnalysisDataset::GetTrigEffErr(const char* Process , 
				       int         McFlag    , 
				       int         RunRange) {
  Error("GetTrigEffErr","should not be called!");
  return 0.;
}


//-----------------------------------------------------------------------------
// returns:
// MC  : total number of the MC events generated for the corresponding 
//       run range to be used as a denominator for acceptance calculation
// data: 0
// RunRange = -1: calculate sum over all used run ranges
//-----------------------------------------------------------------------------
int TAnalysisDataset::GetNTotalMc(const char* Process , 
				  int         McFlag    , 
				  int         RunRange) {
  int n(0);
  const dataset_t* ds;

  if (McFlag != 0) {
    if (RunRange >= 0) {

      ds   = GetDataset(Process,RunRange);

      if (ds != 0) n = ds->fNTotalMc[fLumiBin];
      else {
	Error("GetNTotalMc",Form("%s McFlag=%i run range %i : dataset is not defined",
				 Process,McFlag,RunRange));
      }
    }
    else {
      for (int i=0; i< fNRunRanges; i++) {
	ds = GetDataset(Process,i);
	if (ds) {
	  if (ds->fNTotalMc[fLumiBin] > 0) {
	    n += ds->fNTotalMc[fLumiBin];
	  }
	}
	else {
	  Error("GetNTotalMc",Form("%s McFlag=%i run range %i : dataset is not defined",
				   Process,McFlag,RunRange));
	}
      }
    }
  }

  return n;
}

//-----------------------------------------------------------------------------
int  TAnalysisDataset::SubstituteRunRange(const char* DsID, int RunRange) {
  Error("SubstituteRunRange",Form("DsID=%s, SHOULD NEVER BE CALLED !!!",DsID));
  return -1;
}


//-----------------------------------------------------------------------------
double  TAnalysisDataset::GetAccCorrFactor(const char* Pname   , 
					   int         McFlag   ,
					   int         RunRange) {
  Error("GetAccCorrFactor","SHOULD NEVER BE CALLED !!!");
  return -1;
}


//-----------------------------------------------------------------------------
double TAnalysisDataset::Luminosity(int RunRange, int LumiBin) { 
  return fRunRange[RunRange].fIntLumi[LumiBin] ; 
}


//-----------------------------------------------------------------------------
void   TAnalysisDataset::Print(const char* Opt) const {

  run_range_t* rr;

  for (int irr=0; irr<fNRunRanges; irr++) {
    rr = fRunRange+irr;

    printf ("---- run range %3i min_run=%6i maxrun=%6i DSID=%-10s %10.2f %10.2f %10.2f %10.2f\n",
	    irr,rr->fMinRun, rr->fMaxRun, rr->fDsID, 
	    rr->fIntLumi[0],rr->fIntLumi[1],rr->fIntLumi[2],rr->fIntLumi[3]);
  }

  fListOfDatasets->Print();
}




//-----------------------------------------------------------------------------
// returns:  cross section in pb
// MC  : cross section *BR, corresponding to the simulated process
// data: 0
//-----------------------------------------------------------------------------
double TAnalysisDataset::GetXSec(const char* Process, int McFlag) {

  double xs(0.);

  if (McFlag == 0) return 0;
//-----------------------------------------------------------------------------
// MC : need cross sections to normalize the backgrounds
//-----------------------------------------------------------------------------
  if      (strcmp(Process,"conv"   ) == 0) {
    xs = 1.;
  }
  else if (strcmp(Process,"dio"         ) == 0) {
    xs = 1.;
  }
  else if (strcmp(Process,"cosmics"     ) == 0) {
    xs = 1.;
  }
  else if (strcmp(Process,"pion_capture") == 0) {
    xs = 1.;
  }
  else {
    printf(" TAnalysisDataset::GetXSec ERROR: unknown process %s MCFLAG=%i\n",Process,McFlag);
    return 100.;
  }

  return xs;
}


//-----------------------------------------------------------------------------
TDsMetadata::TDsMetadata(const char* Name, 
                         const char* Title, 
			 const TAnalysisDataset::dataset_t* Data, 
			 const char* Label): 
  TNamed(Name,Title), fData (Data)
{
  // if label is not specified, use Name to plot 
  if (Label != 0) {
    fLabel = Label;
  }
  else {
    fLabel = Name;
  }
};


//-----------------------------------------------------------------------------
TDsMetadata::~TDsMetadata() {
}
