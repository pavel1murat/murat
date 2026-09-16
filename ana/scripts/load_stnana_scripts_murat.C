//
#include "TInterpreter.h"
#include "murat/ana/scripts/modules.hh"
//-----------------------------------------------------------------------------
int load_stnana_scripts_murat(int Debug = 0) {
  char        macro[200];
  const char* script[] = { 
    //    "global_vars.cc",
    //    "init_geometry.C",
    "bflash.C"      ,  "PWD",
    "calorimeter.C" ,  "PWD",
    "cosmics.C"     ,  "PWD",
    "det.C"         ,  "PWD",
    "dio_calib.C"   ,  "PWD",
    "dose.C"        ,  "PWD",
    "drpc.C"        ,  "PWD",
    "g4val.C"       ,  "PWD",
    "fun.C"         ,  "PWD",
    "genp.C"        ,  "PWD",
    "helix.C"       ,  "PWD",
    "hits.C"        ,  "PWD",
    "momscale.C"    ,  "PWD",
    "mustop.C"      ,  "PWD",
    "pid.C"         ,  "PWD",
    "pipenu.C"      ,  "PWD",
    "reco_eff_ana.C",  "PWD",
    //    "rmc.C"  ,  "PWD",
    "simp.C"        ,  "PWD",
    "spmc.C"        ,  "PWD",
    "track.C"       ,  "PWD",
    "trigger.C"     ,  "PWD",
    "vdet.C"        ,  "PWD",
    0 
  };

  TString work_dir = gEnv->GetValue("Stnana.TestReleaseDir",gSystem->Getenv("PWD"));

  TInterpreter* cint = gROOT->GetInterpreter();

  std::string   pkg("murat");
  
  for (int i=0; script[i] != 0; i+=2) {
    const char* dir = gSystem->Getenv(script[i+1]);
    if (dir) {
      if (strcmp(script[i+1],"PWD") == 0) sprintf(macro,"%s/%s/ana/scripts/%s",work_dir.Data(),pkg.data(),script[i]);
      else                                sprintf(macro,"%s/%s/ana/scripts/%s",dir            ,pkg.data(),script[i]);
      if (! cint->IsLoaded(macro)) {
        if (Debug != 0) {
          printf("[load_stnana_scripts_murat.C:%03i]  loading %s\n",__LINE__,macro);
        }
        cint->LoadMacro(macro);
      }
    }
  }
  
  return 0;
}
