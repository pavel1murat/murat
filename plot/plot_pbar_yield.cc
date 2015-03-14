//

//
#include "stdio.h"

#include "TMath.h"
#include "TAxis.h"
#include "TPad.h"
#include "TGraph.h"

#include "murat/plot/murat_plot_functions.hh"
#include "murat/plot/pbar_common.hh"

// extern "C" {
//   float pbar_yield_(float* P0, float* P, float* Theta);
// }

// //-----------------------------------------------------------------------------
// float get_pbar_yield(float P0, float P, float Theta) {
//   printf("emoe\n");
//   float xsec = pbar_yield_(&P0,&P,&Theta);
//   return xsec;
// }

//-----------------------------------------------------------------------------
// plot MARS pbar production cross sections
//-----------------------------------------------------------------------------
extern "C" void set_p2max_(float*);

void plot_pbar_yield(float P0, float Theta) {

  //  InitPbarCommon();

  float  p[1000], xs[1000];
  double e, cf;

  int np = 160;

  double const MP = 0.93825;	// in GeV
				// this is 2011-2012 default
  float p2max[2] = {1.3, 2};

  set_p2max_(&p2max[0]);
  //   gPbarCommon->P2MAX = 1.3;

  for (int i=0; i<np; i++) {
    p[i]  = i*0.05+1.e-12; // to avoid zeroes

    // calculate correction factor

    e     = sqrt(p[i]*p[i]+MP*MP);
    cf    = 1.539e6*e/(p[i]*p[i])/(2*TMath::Pi()*TMath::Sin(Theta));
    xs[i] = get_pbar_yield(P0,p[i],Theta)*cf;
    //    printf(" i, p, xs : %3i   %10.3f  %12.5e\n",i,p[i],xs[i]);
  }


  TGraph* g = new TGraph(np,p,xs);

  g->SetName(Form("gr_%.0f_%.0f",P0,Theta*180/TMath::Pi()));
  g->SetTitle(Form("MARS #bar{p} yield, P0 =%5.1f GeV/c, #theta_{lab} =%6.1f deg",
		   P0,Theta*180/TMath::Pi()));

  g->GetXaxis()->SetRangeUser(0,7.999);
  //  g->GetYaxis()->SetRangeUser(1.e-15,5.e2);
  g->GetYaxis()->SetRangeUser(1.e-5,5.e2);
  g->GetYaxis()->SetTitle("Ed^{3}#sigma/dp^{3}, #mub/(GeV/c)^{2}");
  g->GetYaxis()->SetTitleOffset(1.3);

  gPad->SetLogy(1);
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  g->Draw("ALP");
//-----------------------------------------------------------------------------
// now redefine P2MAX and do the same
//-----------------------------------------------------------------------------
				// this is 2013 default
  //  gPbarCommon->P2MAX = 2.0;

  set_p2max_(&p2max[1]);

  for (int i=0; i<np; i++) {
    p[i]  = i*0.05+1.e-12; // to avoid zeroes

    // calculate correction factor

    e     = sqrt(p[i]*p[i]+MP*MP);
    cf    = 1.539e6*e/(p[i]*p[i])/(2*TMath::Pi()*TMath::Sin(Theta));
    xs[i] = get_pbar_yield(P0,p[i],Theta)*cf;
    //    printf(" i, p, xs : %3i   %10.3f  %12.5e\n",i,p[i],xs[i]);
  }

  TGraph* g2 = new TGraph(np,p,xs);

  g2->SetName(Form("gr2_%.0f_%.0f",P0,Theta*180/TMath::Pi()));
  g2->SetTitle(Form("MARS #bar{p} yield, P0 =%5.1f GeV/c, #theta_{lab} =%6.1f deg",
		   P0,Theta*180/TMath::Pi()));

  g2->GetXaxis()->SetRangeUser(0,7.999);
  g2->GetYaxis()->SetRangeUser(1.e-15,5.e2);
  g2->GetYaxis()->SetTitle("Ed^{3}#sigma/dp^{3}, #mub/(GeV/c)^{2}");
  g2->GetYaxis()->SetTitleOffset(1.3);

  g2->SetLineColor(2);
  g2->Draw("LP,same");
}
