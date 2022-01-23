//

#include "murat/ana/InitVirtualDetectors.hh"

///////////////////////////////////////////////////////////////////////////////
// define X-offsets to histogram hits on virtual detectors conveniently
// NDet is not used so far, make sure dimension of VDet array is >= 10
// the magnetic field value is the absolute value of the vector
// use it to convert pT into radius
///////////////////////////////////////////////////////////////////////////////
namespace murat {

int InitVirtualDetectors(murat::VDetData_t* VDet, int* NDet) {
//            id ix iz       x         z       bfield      phi
//-----------------------------------------------------------------------------
  VDet[ 0] = { 0, 0, 0,     0.00,     0.00,   0.0      ,   0.0 };       // VD0  : unused
  VDet[ 1] = { 1, 1, 3,  3904.00, -4003.99,   2.6627833,   0.0 };       // VD1  : upstream   TS1
  VDet[ 2] = { 1, 1, 3,  3904.00, -3004.01,   2.3327647,   0.0 };       // VD2  : downstream TS1
  VDet[ 3] = { 3, 3, 1,   824.99,     0.00,   2.4380166, -90.0 };       // VD3  : upstream   TS31
  VDet[ 4] = { 4, 3, 1,    25.01,     0.00,   2.2344449, -90.0 };       // VD4  : downstream TS31
  VDet[ 5] = { 5, 3, 1,   -25.01,     0.00,   2.2245229, -90.0 };       // VD5  : upstream   TS32
  VDet[ 6] = { 6, 3, 1,  -824.99,     0.00,   2.0705379, -90.0 };       // VD6  : downstream TS32
  VDet[ 7] = { 7, 1, 3, -3904.00,  2930.01,   2.2008908,   0.0 };       // VD7  : upstream   TS5
  VDet[ 8] = { 8, 1, 3, -3904.00,  3929.99,   1.9731202,   0.0 };       // VD8  : downstream TS5
  VDet[ 9] = { 9, 1, 3, -3904.00,  5465.99,   1.5943823,   0.0 };       // VD9  : ST_In
  VDet[10] = {10, 1, 3, -3904.00,  6176.01,   1.3833433,   0.0 };       // VD10 : ST_Out
  VDet[11] = {11, 1, 3, -3904.00, 10181.10,   1.0116359,   0.0 };       // VD11 : TT_Mid
  VDet[12] = {12, 1, 3, -3904.00, 10181.10,   1.0116359,   0.0 };       // VD12 : TT_MidInner
  VDet[13] = {13, 1, 3, -3904.00,  8539.89,   1.0371811,   0.0 };       // VD13 : TT_FrontHollow

					// dummy initialization (until needed)
  for (int i=14; i<98; i++) {
    VDet[i] = {i, 0, 0,     0.0,      0.0,    1.,       0.0 };
  }

  VDet[91] = {91, 1, 3,  3904.00, -4144.28 ,  2.7786501,   0.0 };       // VD91 : in front of the first PBAR window
  VDet[92] = {92, 1, 3,  3904.00, -4144.01 ,  2.7769741,   0.0 };       // VD92 : right after the first PBAR window

  VDet[98] = {98, 1, 3,  3046.12,  -857.884,  2.2879650, -45.0 };       // VD98 : mid-bend TSu 
  VDet[99] = {99, 1, 3, -3046.12,   857.884,  2.1390811, -45.0 };       // VD99 : mid-bend TSd

  *NDet            = 100;

  return 0;
}
}
