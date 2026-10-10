#ifndef __det_pos_hh__
#define __det_pos_hh__

namespace murat {
class DetPos_t {
public:
                                        // keep dimensions explicit
  float  trk_pos [3];                   // global position
  float  trk_off [3];                   // offset wrt above
  float  crv_off [3];                   // one global offset wrt nominal position
  float  calo_pos[2][3];                // two disks separately
  float  calo_off[2][3];

  int    GetFromRunDb(int RunNumber, int PrintLevel);
};
}
#endif
