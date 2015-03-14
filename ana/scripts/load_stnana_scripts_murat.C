//
#include "TInterpreter.h"
#include "murat/ana/scripts/modules.hh"
//-----------------------------------------------------------------------------
int load_stnana_scripts_murat() {
  char        macro[200];
  const char* script[] = { 
#    "global_vars.cc",
#    "init_geometry.C",
    "cosmics.C",
    "dio_calib.C",
    "reco_eff_ana.C",
    "track.C",
    "vdet.C",
    0 
  };

  const char* work_dir = gSystem->Getenv("WORK_DIR");

  TInterpreter* cint = gROOT->GetInterpreter();
  
  for (int i=0; script[i] != 0; i++) {
    sprintf(macro,"%s/murat/ana/scripts/%s",work_dir,script[i]);
    if (! cint->IsLoaded(macro)) {
      cint->LoadMacro(macro);
    }
  }
  
  return 0;
}
