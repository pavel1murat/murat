///////////////////////////////////////////////////////////////////////////////
#ifndef murat_obj_datasets_one_rr_hh
#define murat_obj_datasets_one_rr_hh

#include "murat/obj/TAnalysisDataset.hh"

namespace {
  const int kNRunRanges = 1;

  TAnalysisDataset::run_range_t  run_range[1] = { 
//-----------------------------------------------------------------------------
//    rmin    rmax     DSID  L3 path tag  int_lumi           int_lumi
//                            (TAU_MET)   jaco_corr            rs_db
//                                         (1.019)
//---------------- --------- ---------------------------------------------------
         1, 100000,    "??",      -1,         1.,      1.,     1.,     1. //   0 
  };

//--------------------------------------------------------------------------------------------------------------------------
// datasets with names                                   rmin  rmax   name      grl     mcflag  NEVENTS 
//--------------------------------------------------------------------------------------------------------------------------
  TAnalysisDataset::dataset_t  data01[kNRunRanges] = {      1,500000,"data01","DQM_V34", 0,      1,      1,      1,      1 };
  TAnalysisDataset::dataset_t  conv01[kNRunRanges] = {      1,500000,"conv01","DQM_V34", 3,      1,      1,      1,      1 };
  TAnalysisDataset::dataset_t  cosm01[kNRunRanges] = {      1,500000,"cosm01","DQM_V34", 3,      1,      1,      1,      1 };
  //  TAnalysisDataset::dataset_t  dio01 [kNRunRanges] = {      1,500000,"dio01" ,"DQM_V34", 3,      1,      1,      1,      1 };
  TAnalysisDataset::dataset_t  dio02 [kNRunRanges] = {      1,500000,"dio02" ,"DQM_V34", 3,      1,      1,      1,      1 };
};

#endif
