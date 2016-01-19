/* /////////////////////////////////////////////////////////////////////////////

to run: 
-------
root.exe
.L lib/libmu2e_CalorimeterGeom.so
.L test_crystal_map.C++
test1(5)
test2()

///////////////////////////////////////////////////////////////////////////// */

#include "CalorimeterGeom/inc/HexMapper.hh"
#include "CalorimeterGeom/inc/SquareMapper.hh"


//-----------------------------------------------------------------------------
// print out the crystal numerology
//-----------------------------------------------------------------------------
int test1(int N) {

  mu2e::SquareMapper map;
  
  for (int i=0; i<N; i++) {
     mu2e::SquLK x = map.lk(i);

    printf("i, x._i, x._y = %5i %5i %5i\n",i, x._l, x._k);
  }
  
  return 0;
}



//-----------------------------------------------------------------------------
// compare hex and square map definitions of nRings
//-----------------------------------------------------------------------------
int test2() {

  mu2e::HexMapper    hm;
  mu2e::SquareMapper sm;

  printf(" hm.nCrystalMax(1) = %5i sm.nCrystalMax(1) = %5i\n",
	 hm.nCrystalMax(1), sm.nCrystalMax(1));

  return 0;
}
