//
#include "math.h"
#include "TMath.h"
#include "TGraph.h"
#include "plot/pbar_common.hh"

//-----------------------------------------------------------------------------
// p -in GeV
//-----------------------------------------------------------------------------
void plot_pbar_kinematics(float P0) {

  InitPbarCommon();

  double thcm[200], px[200], py[200] , plab[200], thlab[200];

  double MP = 0.938;


  double e0 = sqrt(P0*P0+MP*MP);
				// center of mass velocity

  double v0 = P0/(e0+MP);

  double s = (e0+MP)*(e0+MP)-P0*P0;

  double sqrts = sqrt(s);

  double q = sqrts-2*MP; // 2 times pbar energy

  // max antiproton momentum in the CM system

  double e = q/2.;
  double p = sqrt(e*e - MP*MP);

  double vcm = sqrt(1-s/(e0*e0));

  printf(" p0, e0, sqrts, vcm, e, p: %10.3f %10.3f %10.3f %10.5f %10.3f %10.3f\n",
	 P0, e0, sqrts, vcm, e, p);

  int np = 181;
  for (int i=0; i<np; i++) {

    thcm[i] = i*TMath::Pi()/180;
				// px,py - in the lab system

    px[i] = (p*cos(thcm[i])+vcm*e)/sqrt(1-vcm*vcm);
    py[i] =  p*sin(thcm[i]);

    thlab[i] = atan2(py[i],px[i])*180/TMath::Pi();
    plab [i] = sqrt(px[i]*px[i]+py[i]*py[i]);

    printf("i, thcm, px, py, thlab, plab = %3i %10.3f %10.3f %10.3f %10.3f %10.3f\n",
	   i,thcm[i],px[i],py[i],thlab[i],plab[i]);
  }

  TGraph* g = new TGraph(np, thlab,plab);

  g->Draw("ALP");
}
