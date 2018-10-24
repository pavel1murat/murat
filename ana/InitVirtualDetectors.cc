//

#include "murat/ana/InitVirtualDetectors.hh"

///////////////////////////////////////////////////////////////////////////////
// define X-offsets to histogram hits on virtual detectors conveniently
// NDet is not used so far, make sure dimension of VDet array is >= 10
///////////////////////////////////////////////////////////////////////////////
int InitVirtualDetectors(VDetData_t* VDet, int* NDet) {

  VDet[ 0] = { 0, 0, 0,     0.0, 0.0    };       // VD0  : unused
  VDet[ 1] = { 1, 1, 3,  3904.0, 2.6654 };       // VD1  : upstream   TS1
  VDet[ 2] = { 1, 1, 3,  3904.0, 2.3324 };       // VD2  : downstream TS1
  VDet[ 3] = { 3, 3, 1,     0.0, 2.4402 };       // VD3  : upstream   TS31
  VDet[ 4] = { 4, 3, 1,     0.0, 2.2162 };       // VD4  : downstream TS31
  VDet[ 5] = { 5, 3, 1,     0.0, 2.2063 };       // VD5  : upstream   TS32
  VDet[ 6] = { 6, 3, 1,     0.0, 2.0729 };       // VD6  : downstream TS32
  VDet[ 7] = { 7, 1, 3, -3904.0, 2.2019 };       // VD7  : upstream   TS5
  VDet[ 8] = { 8, 1, 3, -3904.0, 1.9333 };       // VD8  : downstream TS5
  VDet[ 9] = { 9, 1, 3, -3904.0, 1.5884 };       // VD9  : upstream   ST
  VDet[10] = {10, 1, 3, -3904.0, 1.3799 };       // VD10 : downstream ST

  *NDet            = 11;

  return 0;
}
