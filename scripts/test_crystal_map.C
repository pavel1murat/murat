/* /////////////////////////////////////////////////////////////////////////////

to run: 
-------
root.exe
.L lib/libmu2e_CalorimeterGeom.so
.L test_crystal_map.cc+
test_crystal_map()

///////////////////////////////////////////////////////////////////////////// */

#include "CalorimeterGeom/inc/HexMapper.hh"
#include "CalorimeterGeom/inc/SquareMapper.hh"

int test_crystal_map() {

  mu2e::HexMapper    hm;
  mu2e::SquareMapper sm;


  printf(" hm.nCrystalMax(1) = %5i sm.nCrystalMax(1) = %5i\n",
	 hm.nCrystalMax(1), sm.nCrystalMax(1));

  return 0;
}
