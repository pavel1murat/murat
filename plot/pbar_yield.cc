//
#include "stdio.h"

#include "murat/plot/pbar_yield.hh"


extern "C" {
  float pbar_yield_(float* P0, float* P, float* Theta);
}


//-----------------------------------------------------------------------------
float get_pbar_yield(float P0, float P, float Theta) {
  printf("emoe\n");
  float xsec = pbar_yield_(&P0,&P,&Theta);
  return xsec;
}
