//
#include "TInterpreter.h"
#include "murat/ana/scripts/modules.hh"
//-----------------------------------------------------------------------------
int load_stnana_scripts_murat() {
  char        macro[200];
  const char* script[] = { 
    //    "global_vars.cc",
    //    "init_geometry.C",
    "calorimeter.C" ,  "PWD",
    "cosmics.C"     ,  "PWD",
    "dio_calib.C"   ,  "PWD",
    "genp.C"        ,  "PWD",
    "helix.C"       ,  "PWD",
    "hits.C"        ,  "PWD",
    "pid.C"         ,  "PWD",
    "reco_eff_ana.C",  "PWD",
    //    "rmc.C"  ,  "PWD",
    "simp.C"        ,  "PWD",
    "track.C"       ,  "PWD",
    "trigger.C"     ,  "PWD",
    "vdet.C"        ,  "PWD",
    0 
  };

  TString work_dir = gEnv->GetValue("Stnana.TestReleaseDir",gSystem->Getenv("PWD"));

  TInterpreter* cint = gROOT->GetInterpreter();
  
  for (int i=0; script[i] != 0; i+=2) {
    const char* dir = gSystem->Getenv(script[i+1]);
    if (dir) {
      sprintf(macro,"%s/murat/ana/scripts/%s",dir,script[i]);
      if (! cint->IsLoaded(macro)) cint->LoadMacro(macro);
    }
  }
  
  return 0;
}
