///////////////////////////////////////////////////////////////////////////////
// book-keeping class for MU2E analysis: datasets for Mu2e channel
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
// simple counting number of events in the MC
// murat/jobs/shell_scripts/parse_catalog_log /cdf/ncdf192/datasets/we0sge/catalog/ 210012 212133 | grep '#'
//
// to calculate numb of events for different good run list incarnations :
//
// root [4] .L murat/ana/scripts/list_of_runs.C 
// root [5] check_list_of_runs("/cdf/ncdf192/datasets/catalog/ze1s9m.catalog","DQM_V32")
///////////////////////////////////////////////////////////////////////////////
#include "TSystem.h"
#include "TEnv.h"

#include "obj/TMu2eDatasets.hh"
#include "obj/mu2e_datasets_one_rr.hh"

ClassImp(TMu2eDatasets)

//-----------------------------------------------------------------------------
// 'Title' is used to choose versioned datasets. 
//  So far the only choice is that of the Z-->ee MC: 
// "ze0s01" (GEN7) or "ze1sxd" (GEN6, not complete)
// it is obviously a very kludgy way to choose between different options....
//-----------------------------------------------------------------------------
namespace {
}

//-----------------------------------------------------------------------------
TMu2eDatasets::TMu2eDatasets(const char* Name, const char* Version, Int_t LumiBin): 
  TAnalysisDataset(Name,Version,LumiBin) 
{
//-----------------------------------------------------------------------------
// low-statistics analysis, have only one run range
//-----------------------------------------------------------------------------
  fNRunRanges = 1; 
  fRunRange   = run_range;
//-----------------------------------------------------------------------------
// data, use emu-nosi luminosity
// name and dataset name are the same
//-----------------------------------------------------------------------------
  fListOfDatasets->Add(new TDsMetadata("data"   ,"data01",data01,"data"));
  fListOfDatasets->Add(new TDsMetadata("conv"   ,"conv01",conv01,"CONV01"));
  fListOfDatasets->Add(new TDsMetadata("cosmics","cosm01",cosm01,"COSM01"));
  //  fListOfDatasets->Add(new TDsMetadata("dio"   ,"dio01" ,dio01 ,"DIO01" ));

  fListOfDatasets->Add(new TDsMetadata("dio"    ,"dio02" ,dio02 ,"DIO02" ));
}


//-----------------------------------------------------------------------------
TMu2eDatasets::~TMu2eDatasets() {
}

//-----------------------------------------------------------------------------
// naming conventions for the histogram files are defined here
// histograms for WZ3l analysis have been created by zzx_wz3l job
//-----------------------------------------------------------------------------
int TMu2eDatasets::GetHistFileName(const char* Process, 
				   int         McFlag, 
				   int         RunRange, 
				   char*       HistFileName) 
{
  char             dirname[200];
  //  const char       *nm;
  const dataset_t *ds;

  ds = GetDataset(Process,RunRange);
//-----------------------------------------------------------------------------
// default directory 
//----------------------------------------------------------------------------- 
  sprintf(dirname,"%s",gEnv->GetValue("mu2e.HistDir","UNDEFINED"));

  HistFileName[0] = 0;

  if (ds != 0) {
    sprintf(HistFileName,"%s/limits/%s.mu2e_limits.hist",dirname,ds->GetDsID());
  }

  return 0;
}


//-----------------------------------------------------------------------------
// because of the specifics of the MC generation min. and max. run numbers 
// may be defined differently for the same run range for different datasets
// 'taumet' : data only
//-----------------------------------------------------------------------------
int TMu2eDatasets::GetHistColor(const char* Process, int McFlag) const {
  int color(-1);

  if      (strcmp(Process,"data") == 0) {
    color = 1;
  }
  else if (strcmp(Process,"dio"   ) == 0) {
    color = 42; // 4;  // 7
  }
  else if (strcmp(Process,"conv"   ) == 0) {
    color = 4; // 42;  // 7
  }
  else {
    Error("GetHistColor",Form("Unknown process: %s, return 1",Process));
    color = 1;
  }

  return color;
}


//-----------------------------------------------------------------------------
// because of the specifics of the MC generation min. and max. run numbers 
// may be defined differently for the same run range for different datasets
// 'taumet' : data only
//-----------------------------------------------------------------------------
int TMu2eDatasets::GetHistStyle(const char* Process, int McFlag) const {
  int style(3001);

  if      (strcmp(Process,"data") == 0) {
    style = 1;
  }
  else if (strcmp(Process,"dio"   ) == 0) {
    style = 3001; // 3001;
  }
  else if (strcmp(Process,"dio"   ) == 0) {
    style = 3013; // 3001;
  }
  else {
    Error("GetHistStyle",Form("Unknown Process: %s, return 1",Process));
    style = 1;
  }

  return style;
}

//-----------------------------------------------------------------------------
double TMu2eDatasets::GetTightEleIDSF(const char* Process, int McFlag, int RunRange) {

  Error("GetTightEleIDSF", "Was called, returning 1");
  return 1;

}

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
double TMu2eDatasets::GetLooseEleIDSF(const char* Process, int McFlag, int RunRange) {
  double sf(1.);

  Error("GetLooseEleIDSF", "Was called, returning 1");

  return sf;
}

//-----------------------------------------------------------------------------
// TAU_MET trigegr efficiency for because of the specifics of the MC 
// generation min. and max. run numbers 
// may be defined differently for the same run range for different datasets
// !!! scale factor may be different for 1P and 3P taus....
// !!! add Bin to the parameterization

//   i, eff err =  0    0.86816    0.01168
//   i, eff err =  1    1.11432        nan
//   i, eff err =  2    1.10197        nan
//   i, eff err =  3    1.06244        nan
//   i, eff err =  4    1.09644        nan
//   i, eff err =  5    0.84106    0.02070
//   i, eff err =  6    0.89914    0.02080
//   i, eff err =  7    0.90295    0.03600
//   i, eff err =  8    0.90647    0.01502
//   i, eff err =  9    0.88848    0.01688
//   i, eff err = 10    0.90094    0.01317
//   i, eff err = 11    0.88623    0.01508
//   i, eff err = 12    0.84097    0.02103
//   i, eff err = 13    0.87019    0.01410
//   i, eff err = 14    0.83873    0.03597
//   i, eff err = 15    0.88838    0.01648
//   i, eff err = 16    0.87662    0.02258
//   i, eff err = 17    0.86447    0.01818
//   i, eff err = 18    0.87485    0.01143
//   i, eff err = 19    0.90534    0.01427
//   i, eff err = 20    0.85985    0.01662
//   i, eff err = 21    0.86524    0.01159
//-----------------------------------------------------------------------------
double TMu2eDatasets::GetTriggerEffSF(const char* Process, int McFlag, int RunRange) {

  Error("GetTriggerEffSF", "Was called, returning 1");
  return 1;

  double sf(0.);

  // double trig_eff_sf [] = { 
  //  0.86459  , 0.85781  , 0.86271  , 0.89026  , 0.89284  ,
  //  0.82183  , 0.86722  , 0.84179  , 0.87564  , 0.87456  ,
  //  0.89379  , 0.87990  , 0.85997  , 0.86215  , 0.82487  ,
  //  0.86157  , 0.90768  , 0.86197  , 0.86859  , 0.88281  ,
  //  0.86825  , 0.85534  ,-1.0      ,-1.       , -1       ,

//     0.86816, 0.86816, 0.86816, 0.86816, 0.86816,
//     0.84106, 0.89914, 0.90295, 0.90647, 0.88848,
//     0.90094, 0.88623, 0.84097, 0.87019, 0.83873,
//     0.88838, 0.87662, 0.86447, 0.87485, 0.90534,
//     0.85985, 0.86524, -1.    , -1.    , -1.0   ,
  //   -1
  // };

  if  (McFlag == 0) {
    sf = 1.;
  }
  else {
//-----------------------------------------------------------------------------
// MC : data-to-MC scale factors
//-----------------------------------------------------------------------------
				// efficiencies
//     if      (strcmp(Process,"wenu"   ) == 0) sf = trig_eff_sf[RunRange];
//     else if (strcmp(Process,"wmunu"  ) == 0) sf = trig_eff_sf[RunRange];
//     else if (strcmp(Process,"wtaunu" ) == 0) sf = trig_eff_sf[RunRange];
//     else if (strcmp(Process,"zee"    ) == 0) sf = trig_eff_sf[RunRange];
//     else if (strcmp(Process,"zmumu"  ) == 0) sf = trig_eff_sf[RunRange];
//     else if (strcmp(Process,"ztautau") == 0) sf = trig_eff_sf[RunRange];
//     else if (strcmp(Process,"qcd"    ) == 0) sf = trig_eff_sf[RunRange];
    sf = 1;
  }

  return sf;
}

//-----------------------------------------------------------------------------
double TMu2eDatasets::GetTrigEffErr(const char* Process, int McFlag, int RunRange) {

  double err(0.);

  Error("GetTriggerEffErr", "Was called, returning 0");
  return 0;

  double trig_eff_sf_err [] = { 
    // points 2-5 need to be updated

    0.01168,  0.01168,  0.01168,  0.01168,  0.01168,  
    0.02070,  0.02080,  0.03600,  0.01502,  0.01688,
    0.01317,  0.01508,  0.02103,  0.01410,  0.03597,
    0.01648,  0.02258,  0.01818,  0.01143,  0.01427,
    0.01662,  0.01159,
    -1
  };

 
  if (McFlag == 0) {
    err = 0.;
  }
  else {
//-----------------------------------------------------------------------------
// MC : data-to-MC scale factors
//-----------------------------------------------------------------------------
    if      (strcmp(Process,"wenu"    ) == 0) err = trig_eff_sf_err[RunRange];
    else if (strcmp(Process,"wmunu"   ) == 0) err = trig_eff_sf_err[RunRange];
    else if (strcmp(Process,"wtaunu"  ) == 0) err = trig_eff_sf_err[RunRange];
    else if (strcmp(Process,"zee_jets") == 0) err = trig_eff_sf_err[RunRange];
    else if (strcmp(Process,"zmm_jets") == 0) err = trig_eff_sf_err[RunRange];
    else if (strcmp(Process,"ztautau" ) == 0) err = trig_eff_sf_err[RunRange];
    else if (strcmp(Process,"qcd"     ) == 0) err = trig_eff_sf_err[RunRange];
  }

  return err;
}


//-----------------------------------------------------------------------------
// DsID: 'ze1sxd' , not 'ze1s6d' or 'ze1sad'...
//-----------------------------------------------------------------------------
int TMu2eDatasets::SubstituteRunRange(const char* Process, int RunRange) {
  int         rr;
				        // default: no substitution
  rr = RunRange;

  return rr;
}


//-----------------------------------------------------------------------------
// this is analysis specific
//-----------------------------------------------------------------------------
double TMu2eDatasets::GetAccCorrFactor(const char* Process , 
				       int         McFlag  ,
				       int         RunRange) {
  double corr_factor;
//-----------------------------------------------------------------------------
// MC: corr_factor = N(66<M(ll)<116)/N(tot)
//-----------------------------------------------------------------------------
  const char* dsid = GetDataset(Process,RunRange)->GetDsID();

  corr_factor = 1.;

  if      (strcmp(dsid,"ze1s01") == 0) {
					// ZE0S01 is generated with the M(ll) > 20
    corr_factor = 1.96;
  }
  else if (strcmp(dsid,"ze1s6d") == 0) {
    corr_factor = 1.96;
  }
  else if (strcmp(dsid,"ze1s9m") == 0) {
    corr_factor = 1.96;
  }
  else if (strcmp(dsid,"ze1sxm") == 0) {
    corr_factor = 1.96;
  }
  else if (strcmp(dsid,"ze1sad") == 0) {
    corr_factor = 1.96;
  }

  printf("[TMu2eDatasets::GetAccCorrFactor] dsid = %s, corr_factor = %12.5f\n",dsid,corr_factor);

  return corr_factor;
}



//-----------------------------------------------------------------------------
// results of the gaussian fit: TCdf10008_zee_xsec.cc: plot 51
//-----------------------------------------------------------------------------
double TMu2eDatasets::GetCentralEtSF(const char* Process, int McFlag, int RunRange) {

  double sf(0.);

  sf = 1.;

  printf(" ERROR: dummy TMu2eDatasets::GetCentralEtSF should not be called\n");
  return sf;


  double central_et_sf_bhelzl [] = { 
    1.0038, 1.0069, 1.0065, 1.0053, 1.0083, 1.0085, 1.0055, 1.0072, 1.0075, 1.0090,
    1.0088, 1.0103, 1.0124, 1.0116, 1.0082, 1.0112, 1.0114, 1.0113, 1.0039, 1.0039,
    1.0051, 1.0058, 1.0033, 1.0043, 1.0030, 1.0038, 1.0048,
    -1
  };

//   Error("GetTriggerEffErr", "Was called, returning 0");
//   return 0;

  if      (strcmp(Process,"bhelzl"  ) == 0) {
    sf = central_et_sf_bhelzl[RunRange];
  }
  
  return sf;
}

//-----------------------------------------------------------------------------
void TMu2eDatasets::Help(int Mode) {

  static char HelpText[3][100];
  static int  initialized(0);

  if (initialized == 0) {
    strcpy(HelpText[0],"known versions: gpt000: s-channel production of the RS G*, PYTHIA signal");
    strcpy(HelpText[1],"                gpt100: production of pp --> G* q/g. PT>100");
    strcpy(HelpText[2],"");
    initialized = 1;
  };

  for (int i=0; HelpText[i]!=0; i++) {
    printf("%s\n",HelpText[i]);
  }

}
