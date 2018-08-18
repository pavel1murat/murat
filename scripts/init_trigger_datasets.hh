///////////////////////////////////////////////////////////////////////////////
// location of the datasets: catalogs/beamline_dose
///////////////////////////////////////////////////////////////////////////////
#ifndef __murat_scripts_init_trigger_datasets__
#define __murat_scripts_init_trigger_datasets__

#include "murat/scripts/dataset.hh"
//-----------------------------------------------------------------------------
// datasets
//-----------------------------------------------------------------------------
dataset_t  d_trig_645_0001;
dataset_t  d_fast_645_0001;
dataset_t  d_trig_645_0002;
dataset_t  d_fast_645_0002;
dataset_t  d_trig_645_0003;
dataset_t  d_fast_645_0003;

int        trigger_datasets_initialized(0);
//-----------------------------------------------------------------------------
void init_trigger_datasets() {
  const char* HistDir         = "/projects/hist/mu2e/v6_4_1";

  //  const char* DeltaFinderHistDir = "/projects/mu2e/results/deltaFinder";

  if (trigger_datasets_initialized == 1) return;
  
  d_trig_645_0001.fName  = "trig_645_0001";
  d_trig_645_0001.fFn    = Form("%s/validation_645_0001.calpatrec_trigger.trigger_ana.hist",HistDir);
  d_trig_645_0001.fLabel = "trig_645_0001";
  d_trig_645_0001.fNPOT  = 1.e6;

  d_fast_645_0001.fName  = "fast_645_0001";
  d_fast_645_0001.fFn    = Form("%s/validation_645_0001.calpatrec_trigger_fast.trigger_ana.hist",HistDir);
  d_fast_645_0001.fLabel = "fast_645_0001";
  d_fast_645_0001.fNPOT  = 1.e6;
  
  trigger_datasets_initialized  = 1;
}

#endif
