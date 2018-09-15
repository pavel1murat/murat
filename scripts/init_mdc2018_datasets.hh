//
#ifndef __murat_scripts_init_mdc2018_datasets__
#define __murat_scripts_init_mdc2018_datasets__

#include "murat/scripts/dataset.hh"
//-----------------------------------------------------------------------------
// MDC2018 datasets (with corresponding histogram files)
//-----------------------------------------------------------------------------
dataset_t  d_mdc2018_0001;   //  CeEndpoint-mix.trackerMCCheck_read
//-----------------------------------------------------------------------------
void init_mdc2018_datasets() {
  
  const char* HistDir    = "/projects/hist/mu2e/v7_0_5";
  
  d_mdc2018_0001.fName = "d_mdc2018_0001";
  d_mdc2018_0001.fFn    = Form("%s/CeEndpoint-mix.MDC2018a.trackerMCCheck_read.hist",HistDir);
  d_mdc2018_0001.fLabel = "";
  d_mdc2018_0001.fNPOT  = 1.; // undefined

}

#endif
