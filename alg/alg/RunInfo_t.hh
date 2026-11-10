#ifndef __runinfo_t_hh__
#define __runinfo_t_hh__

#include <string>
#include <vector>

namespace murat {
  class RunInfo_t {
  public:
                                        // run types
    enum {
      e_BEAM            = 1,
      e_COSMICS         = 2,
      e_INTERNAL_PULSER = 3,
      e_PULSE_INJECTION = 4,
      e_NOISE           = 5,
      e_MONTE_CARLO     = 100,
    };
  
    int              run_number;        //
    int              duration;          // run duration, seconds
    int              max_ewtag;         // max event window tag
    int              run_type;          // 1:beam 2:cosmics 3:internal_pulser
                                        // 4:pulse_injection 5:noise 100:MC
    int              event_mode;        // on/off spill
    int              cfo_rate;          // event rate, Hz
    int              trigger_rate;      // event rate to disk, Hz
    float            ew_length;         // event window length, sec
    int              ref_channel;       // for pulse injection runs
    std::vector<int> pulsed_channel;    // assume one per 8 [12 el]
    std::vector<int> plane_flag;        // should have 36 elements

    int    GetFromRunDb(int RunNumber, int PrintLevel);
    static int DurationToSeconds(const std::string& Hms);   // "hh:mm:ss" -> s, -1 on error
  };
}
#endif
