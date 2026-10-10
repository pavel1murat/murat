//
#include <iostream>
#include <format>

#include "Stntuple/alg/TRunDb.hh"
#include "murat/alg/DetPos_t.hh"

namespace murat {
//-----------------------------------------------------------------------------
int DetPos_t::GetFromRunDb(int RunNumber,int PrintLevel) {
  int rc(0);
  
  auto* rundb = stntuple::TRunDb::Instance();

  fhicl::ParameterSet pset = rundb->GetTable("det_pos",RunNumber);

  // for now, assume success, TBF

  std::vector<float> trkp, trko, crvo;
  
  trkp = pset.get<std::vector<float>>("trk_pos");
  trko = pset.get<std::vector<float>>("trk_off");
  crvo = pset.get<std::vector<float>>("crv_off");
  for (int i=0; i<3; i++) {
    trk_pos[i] = trkp[i];
    trk_off[i] = trko[i];
    crv_off[i] = crvo[i];
  }

  std::vector<float> calp, calo;
  calp = pset.get<std::vector<float>>("calo_pos");
  calo = pset.get<std::vector<float>>("calo_off");

  for (int i=0; i<6; i++) {
    int idisk = i / 3;
    int ix    = i % 3;
    calo_pos[idisk][ix] = calp[i];
    calo_off[idisk][ix] = calo[i];
  }

  if (PrintLevel != 0) {
    std::cout << std::format("trk_pos    : {:10.3f} {:10.3f} {:10.3f}\n",
                             trk_pos[0],trk_pos[1],trk_pos[2]);
    std::cout << std::format("calo_pos[0]: {:10.3f} {:10.3f} {:10.3f}\n",
                             calo_pos[0][0],calo_pos[0][1],calo_pos[0][2]);
    std::cout << std::format("calo_pos[1]: {:10.3f} {:10.3f} {:10.3f}\n",
                             calo_pos[1][0],calo_pos[1][1],calo_pos[1][2]);
  }
  
  return rc;
}

}
