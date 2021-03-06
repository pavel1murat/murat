//-----------------------------------------------------------------------------
// see catalogs for more information about the datasets
// catalogs on murat06: mu2e/app/users/murat/catalogs
//-----------------------------------------------------------------------------
#ifndef __murat_scripts_init_datasets__
#define __murat_scripts_init_datasets__

#include "murat/scripts/init_tracker_dose_datasets.hh"
#include "murat/scripts/init_harp_datasets.hh"
#include "murat/scripts/init_pion_yields_datasets.hh"
#include "murat/scripts/init_beamline_dose_datasets.hh"
#include "murat/scripts/init_delta_finder_datasets.hh"
#include "murat/scripts/init_trigger_datasets.hh"
#include "murat/scripts/init_mdc2018_datasets.hh"

//-----------------------------------------------------------------------------
void init_datasets() {
  init_tracker_dose_datasets();
  init_pion_yields_datasets();
  init_harp_datasets();
  init_beamline_dose_datasets();
  init_delta_finder_datasets();
  init_trigger_datasets();
  init_mdc2018_datasets();
}
#endif
