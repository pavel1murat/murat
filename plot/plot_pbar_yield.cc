//

//
#include "stdio.h"

#include "TMath.h"
#include "TAxis.h"
#include "TPad.h"
#include "TGraph.h"
#include "TCanvas.h"

#include "Stntuple/alg/TStntuple.hh"

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

void plot_pbar_yield(float PBeam, float Theta) {

  //  InitPbarCommon();

  float  p[1000], xs[1000];
  double cf;  // , e ;

  int np = 160;

  //  double const MP = 0.93825;	// in GeV
				// this is 2011-2012 default
  float p2max[2] = {1.3, 2};

  TStntuple::PBar_Striganov_SetP2Max(p2max[0]);

  for (int i=0; i<np; i++) {
    p[i]  = i*0.05+1.e-12; // to avoid zeroes

    // calculate correction factor, 1.539e6 is the inelastic cross section on Ta in mubarn

    //    e     = sqrt(p[i]*p[i]+MP*MP);
    cf    = 1; // 1.539e6*e/(p[i]*p[i])/(2*TMath::Pi()*TMath::Sin(Theta));

    xs[i] = TStntuple::PBar_Striganov_Ed3SigdP3(PBeam,p[i],Theta)*cf;
    //    printf(" i, p, xs : %3i   %10.3f  %12.5e\n",i,p[i],xs[i]);
  }

  TCanvas* c = new TCanvas("c_plot_pbar_yield","Pbar yield",1200,800);
  c->Modified();
  c->Update();

  TGraph* g = new TGraph(np,p,xs);

  g->SetName(Form("gr_%.0f_%.0f",PBeam,Theta*180/TMath::Pi()));
  g->SetTitle(Form("MARS #bar{p} yield, P0 =%5.1f GeV/c, #theta_{lab} =%6.1f deg",
		   PBeam,Theta*180/TMath::Pi()));

  g->GetXaxis()->SetRangeUser(0,7.999);
  //  g->GetYaxis()->SetRangeUser(1.e-15,5.e2);
  g->GetYaxis()->SetRangeUser(1.e-5,5.e2);
  g->GetYaxis()->SetTitle("Ed^{3}#sigma/dp^{3}, #mub/(GeV/c)^{2}");
  g->GetYaxis()->SetTitleOffset(1.3);

  c->SetLogy(1);
  c->SetGridx(1);
  c->SetGridy(1);

  g->Draw("ALP");
//-----------------------------------------------------------------------------
// now redefine P2MAX and do the same
//-----------------------------------------------------------------------------
				// P2MAX = 2.0 is the 2013 default

  TStntuple::PBar_Striganov_SetP2Max(p2max[1]);

  for (int i=0; i<np; i++) {
    p[i]  = i*0.05+1.e-12; // to avoid zeroes

    // calculate correction factor

    //    e     = sqrt(p[i]*p[i]+MP*MP);
    cf    = 1. ;  // 1.539e6*e/(p[i]*p[i])/(2*TMath::Pi()*TMath::Sin(Theta));

    //    xs[i] = get_pbar_yield(P0,p[i],Theta)*cf;
    xs[i] = TStntuple::PBar_Striganov_Ed3SigdP3(PBeam,p[i],Theta)*cf;
    //    printf(" i, p, xs : %3i   %10.3f  %12.5e\n",i,p[i],xs[i]);
  }

  TGraph* g2 = new TGraph(np,p,xs);

  g2->SetName(Form("gr2_%.0f_%.0f",PBeam,Theta*180/TMath::Pi()));
  g2->SetTitle(Form("MARS #bar{p} yield, P0 =%5.1f GeV/c, #theta_{lab} =%6.1f deg",
		   PBeam,Theta*180/TMath::Pi()));

  g2->GetXaxis()->SetRangeUser(0,7.999);
  g2->GetYaxis()->SetRangeUser(1.e-15,5.e2);
  g2->GetYaxis()->SetTitle("Ed^{3}#sigma/dp^{3}, #mub/(GeV/c)^{2}");
  g2->GetYaxis()->SetTitleOffset(1.3);

  g2->SetLineColor(2);
  g2->Draw("LP,same");
}
