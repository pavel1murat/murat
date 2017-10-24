///////////////////////////////////////////////////////////////////////////////
// location of the datasets: catalogs/beamline_dose
///////////////////////////////////////////////////////////////////////////////
#ifndef __murat_scripts_init_delta_finder_datasets__
#define __murat_scripts_init_delta_finder_datasets__

#include "murat/scripts/dataset.hh"
//-----------------------------------------------------------------------------
// datasets
//-----------------------------------------------------------------------------
dataset_t  d_delta_finder_arpan_0035;
dataset_t  d_delta_finder_arpan_0050;
dataset_t  d_delta_finder_arpan_0070;
dataset_t  d_delta_finder_arpan_0200;
dataset_t  d_delta_finder_dave;

int        delta_finder_datasets_initialized(0);
//-----------------------------------------------------------------------------
void init_delta_finder_datasets() {
  const char* HistDir    = "/projects/hist/mu2e/v6_1_4";

  if (delta_finder_datasets_initialized == 1) return;
  
  d_delta_finder_arpan_0035.fName  = "delta_finder_arpan_0035";
  d_delta_finder_arpan_0035.fFn    = Form("%s/delta_finder_arpan_0035.root",HistDir);
  d_delta_finder_arpan_0035.fLabel = "arpan_0035";
  d_delta_finder_arpan_0035.fNPOT  = 1.e6;

  d_delta_finder_arpan_0050.fName  = "delta_finder_arpan_0050";
  d_delta_finder_arpan_0050.fFn    = Form("%s/delta_finder_arpan_0050.root",HistDir);
  d_delta_finder_arpan_0050.fLabel = "arpan_0050";
  d_delta_finder_arpan_0050.fNPOT  = 1.e6;

  d_delta_finder_arpan_0070.fName  = "delta_finder_arpan_0070";
  d_delta_finder_arpan_0070.fFn    = Form("%s/delta_finder_arpan_0070.root",HistDir);
  d_delta_finder_arpan_0070.fLabel = "arpan_0070";
  d_delta_finder_arpan_0070.fNPOT  = 1.e6;

  d_delta_finder_arpan_0200.fName  = "delta_finder_arpan_0200";
  d_delta_finder_arpan_0200.fFn    = Form("%s/delta_finder_arpan_0200.root",HistDir);
  d_delta_finder_arpan_0200.fLabel = "arpan_0200";
  d_delta_finder_arpan_0200.fNPOT  = 1.e6;

  d_delta_finder_dave.fName        = "delta_finder_dave";
  d_delta_finder_dave.fFn          = Form("%s/delta_finder_dave.root",HistDir);
  d_delta_finder_dave.fLabel       = "dave_0035";
  d_delta_finder_dave.fNPOT        = 1.e6;

  delta_finder_datasets_initialized  = 1;
}

#endif
