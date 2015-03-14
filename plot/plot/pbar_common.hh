// -*- Mode: C++ -*- 
//------------------------------------------------------------------------------
// June 01 1996 P.Murat CDF (INFN-Pisa & ITEP-Moscow)
// hepevt.hh : interface to Fortran /HEPEVT/ common-block
//  revision history:
//  ----------------
// *0001 Aug 26 1997 P.Murat: rename and move to generator/
//------------------------------------------------------------------------------

#ifndef __PBAR_COMMON_HH__
#define __PBAR_COMMON_HH__

struct PbarCommon_t {
  int      IFIRST_TIME;        // event number
  float    P2MAX;              // p2max
  float    P2C;                // p2c
};
				// this is what comes from STDHEP

extern PbarCommon_t*    gPbarCommon;

void InitPbarCommon();
//------------------------------------------------------------------------------
//  end of pbar_common.hh
//------------------------------------------------------------------------------
#endif
