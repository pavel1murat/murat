///////////////////////////////////////////////////////////////////////////////
// location of the datasets: catalogs/beamline_dose
///////////////////////////////////////////////////////////////////////////////
#ifndef __murat_scripts_init_beamline_dose_datasets__
#define __murat_scripts_init_beamline_dose_datasets__

#include "murat/scripts/dataset.hh"
//-----------------------------------------------------------------------------
// datasets
//-----------------------------------------------------------------------------
dataset_t  d_beamline_dose_622_0000_ts3;
dataset_t  d_beamline_dose_622_0001_ts1;
dataset_t  d_beamline_dose_622_0001_ts3;

int        beamline_dose_datasets_initialized(0);
//-----------------------------------------------------------------------------
void init_beamline_dose_datasets() {
  const char* HistDir    = "/projects/hist/mu2e/v6_1_4";

  if (beamline_dose_datasets_initialized == 1) return;
  
  d_beamline_dose_622_0000_ts3.fName  = "beamline_dose_622_0000_ts3";
  d_beamline_dose_622_0000_ts3.fFn    = Form("%s/beamline_dose.622_0000.g4s2.ts3coll_stn.coll3_dose_ana.hist",HistDir);
  d_beamline_dose_622_0000_ts3.fLabel = "622_0000_TS3";
  d_beamline_dose_622_0000_ts3.fNPOT  = 1.e6;

  d_beamline_dose_622_0001_ts1.fName  = "beamline_dose_622_0001_ts1";
  d_beamline_dose_622_0001_ts1.fFn    = Form("%s/beamline_dose.622_0001.g4s1.ts1pbarabs_stn.coll1_dose_ana.hist",HistDir);
  d_beamline_dose_622_0001_ts1.fLabel = "622_0001_TS1";
  d_beamline_dose_622_0001_ts1.fNPOT  = 1.e6;

  d_beamline_dose_622_0001_ts3.fName  = "beamline_dose_622_0001_ts3";
  d_beamline_dose_622_0001_ts3.fFn    = Form("%s/beamline_dose.622_0001.g4s2.ts3coll_stn.coll3_dose_ana.hist",HistDir);
  d_beamline_dose_622_0001_ts3.fLabel = "622_0001_TS3";
  d_beamline_dose_622_0001_ts3.fNPOT  = 1.e6;

  beamline_dose_datasets_initialized  = 1;
}

#endif
