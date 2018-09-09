///////////////////////////////////////////////////////////////////////////////
// test syntax of complex structure initialization
///////////////////////////////////////////////////////////////////////////////
#include "cstdio"
#include "TVector3.h"

  struct Data_t {
    TVector3 tsul[2][2];
    TVector3 tsur[2][2];
    TVector3 tsdl[2][2];
    TVector3 tsdr[2][2];
  };

//  positions of the outer and inner points on the cylinder axis
//
  Data_t data = {
//            (x,y,z) outer                   (x,y,z) inner
    {
      {{  0.    ,   0.    ,  0.    }, {  0.0000,  65.5034,  0.0000}},  // nominal
      {{  0.4997,  -0.2817, -0.3827}, {  0.273 ,  65.5034, -0.2076}}   // measured
    },

    {
      {{405.0   , 417.5   ,  0.    }, {296.9065, 417.5   ,  0.0000}},  // nominal
      {{403.0654, 417.0060,  0.3523}, {296.9065, 417.3526,  0.1170}}   // measured
    },
    
    {
      {{  0.    ,   0.    ,  0.0000}, {  0.0000,  54.7010,  0.0000}},  // nominal
      {{ -0.0267,  -1.1981, -0.5653}, { -0.1483,  54.7010, -0.2530}}   // measured
    },
    
    { {{405.0000, 435.0000,  0.0000}, {324.3169, 435.0000,  0.0000}},  // nominal
      {{403.5686, 434.9822,  0.0944}, {324.3169, 435.0884,  0.1290}}   // measured
    }
  };

int test_vinit() {
  printf("data.tsdr[1][1] = { %12.5f %12.5f %12.5f}\n",
	 data.tsdr[1][1].Px(), data.tsdr[1][1].Py(), data.tsdr[1][1].Pz());
  
  return 0;
}
