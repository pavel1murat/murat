//
#include <cstdio>
#include <iostream>
#include <format>
#include <string>

#include "Stntuple/alg/TRunDb.hh"
#include "murat/alg/RunInfo_t.hh"

namespace murat {
//-----------------------------------------------------------------------------
// convert "hh:mm:ss" (hours may exceed 24, minutes and seconds 0..59) to seconds
// returns -1 if the string is malformed
//-----------------------------------------------------------------------------
int RunInfo_t::DurationToSeconds(const std::string& Hms) {
  int hh(0), mm(0), ss(0);
  char extra(0);
  int n = std::sscanf(Hms.c_str(),"%d:%d:%d%c",&hh,&mm,&ss,&extra);
  if (n != 3)                       return -1;
  if (hh < 0 || mm < 0 || mm > 59 || ss < 0 || ss > 59) return -1;
  return hh*3600 + mm*60 + ss;
}
//-----------------------------------------------------------------------------
int RunInfo_t::GetFromRunDb(int RunNumber,int PrintLevel) {
  int rc(0);

  auto* rundb = stntuple::TRunDb::Instance();

  fhicl::ParameterSet pset = rundb->GetTable("runinfo",RunNumber);

  // for now, assume success, TBF

  std::vector<float> trkp, trko, crvo;

  run_number     = pset.get<int>             ("run_number");
  max_ewtag      = pset.get<int>             ("max_ewtag");

  std::string dur_str = pset.get<std::string> ("duration");   // "hh:mm:ss" in the file
  duration       = DurationToSeconds(dur_str);                // stored in seconds
  if (duration < 0) {
    if (PrintLevel > 0) {
      std::cout << std::format("RunInfo_t: bad duration \"{}\" for run {}\n",
                               dur_str,RunNumber);
    }
    rc = -1;
    duration = 0;
  }

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
