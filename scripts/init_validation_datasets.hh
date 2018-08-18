///////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////
#ifndef __murat_scripts_init_validation_datasets__
#define __murat_scripts_init_validation_datasets__

#include "murat/scripts/dataset.hh"
//-----------------------------------------------------------------------------
// datasets
//-----------------------------------------------------------------------------
dataset_t  d_val_641_track_comp_em;
dataset_t  d_val_705_track_comp_em;
dataset_t  d_val_705_track_comp_ep;

int        validation_datasets_initialized(0);
//-----------------------------------------------------------------------------
void init_validation_datasets() {
  
  const char* HistDir         = "/projects/hist/mu2e";

  if (validation_datasets_initialized == 1) return;
  
  d_val_641_track_comp_em.fName  = "val_641_track_comp_em";
  d_val_641_track_comp_em.fFn    = Form("%s/v6_4_1/validation_641_0001_track_comp_11_2_4.hist",HistDir);
  d_val_641_track_comp_em.fLabel = "val_641_track_comp_em";
  d_val_641_track_comp_em.fNPOT  = 1.e6;

  d_val_705_track_comp_em.fName  = "val_705_track_comp_em";
  d_val_705_track_comp_em.fFn    = Form("%s/v7_0_5/eminus_gun_stnmaker.track_comp.hist",HistDir);
  d_val_705_track_comp_em.fLabel = "val_705_track_comp_em";
  d_val_705_track_comp_em.fNPOT  = 1.e6;
  
  d_val_705_track_comp_ep.fName  = "val_705_track_comp_ep";
  d_val_705_track_comp_ep.fFn    = Form("%s/v7_0_5/eplus_gun_stnmaker.track_comp.hist",HistDir);
  d_val_705_track_comp_ep.fLabel = "val_705_track_comp_ep";
  d_val_705_track_comp_ep.fNPOT  = 1.e6;
  
  validation_datasets_initialized  = 1;
}

#endif
