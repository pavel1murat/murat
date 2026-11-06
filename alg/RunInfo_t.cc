//
#include <iostream>
#include <format>

#include "Stntuple/alg/TRunDb.hh"
#include "murat/alg/RunInfo_t.hh"

namespace murat {
//-----------------------------------------------------------------------------
int RunInfo_t::GetFromRunDb(int RunNumber,int PrintLevel) {
  int rc(0);
  
  auto* rundb = stntuple::TRunDb::Instance();

  fhicl::ParameterSet pset = rundb->GetTable("runinfo",RunNumber);

  // for now, assume success, TBF

  std::vector<float> trkp, trko, crvo;
  
  run_number     = pset.get<int>             ("run_number");
  max_ewtag      = pset.get<int>             ("max_ewtag");
  duration       = pset.get<float>           ("duration");
  run_type       = pset.get<int>             ("run_type");
  event_mode     = pset.get<int>             ("event_mode");
  cfo_rate       = pset.get<int>             ("cfo_rate");
  trigger_rate   = pset.get<int>             ("trigger_rate");
  ew_length      = pset.get<float>           ("ew_length");
  ref_channel    = pset.get<int>             ("ref_channel",-1);
  pulsed_channel = pset.get<std::vector<int>>("pulsed_channel",std::vector<int>());
  plane_flag     = pset.get<std::vector<int>>("plane_flag"    ,std::vector<int>());

  return rc;
}

}
