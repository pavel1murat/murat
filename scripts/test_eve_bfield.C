//



#include "murat/ana/TMu2eEveMagField.hh"


//-----------------------------------------------------------------------------
void test() {

  float x(0), y(0), z(0);
  TMu2eEveMagField* mf = new TMu2eEveMagField();

  for (double z = -1000; z<20000; z += 100) {

    TEveVector b = mf->GetField(x,y,z);

    printf("x,y,z,bx,by,bz = %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e \n",
	   x,y,z,b.fX,b.fY,b.fZ);
  }
}
