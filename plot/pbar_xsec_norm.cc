///////////////////////////////////////////////////////////////////////////////
// contains two functions: pbar_xsec_norm and pbar_xsec_norm_2 filling the same
// histogram : p:costh for antiprotons with Striganov's cross section :
//
// a) using function returning d3sigma/dp3
// b) using function returning d2sigma/dp/dcosth
//
// the purpose is to check the normalization by calculating the total antiproton 
// production cross-section
// TStntuple constructor sets Striganov's p2max to 2 GeV/c - our current default
///////////////////////////////////////////////////////////////////////////////
#include "murat/plot/murat_plot_functions.hh"
#include <math.h>
#include "TH2D.h"
#include "TRandom3.h"
#include "TMath.h"
#include "Stntuple/alg/TStntuple.hh"

//------------------------------------------------
void pbar_xsec_norm(double PBeam, long int NEvents) {

  TStntuple* stnt = TStntuple::Instance();

  TRandom3 rn;

  double x[2];

  double pbeam = PBeam;
  double pmax  (5.);

  //  TH2D* h2 = new TH2D("h2","h2",200,0,5,200,-1,1);
  TH2D* h2 = new TH2D("h2_pbar_xsec_norm","pbar_xsec_norm",1000,0,5,200,-1,1);

  for (int i=0; i<NEvents; i++) {

    rn.RndmArray(2,x);

    double plab  = pmax*x[0];
    double costh = 2*x[1]-1;
    double th    = TMath::ACos(costh);
    double w     = stnt->PBar_Striganov_d2N(pbeam,plab,th);

//    printf("pbeam: %12.5f plab: %12.5f thlab: %12.5f w: %12.5e\n",pbeam,plab,th,w);

    h2->Fill(plab,costh,w);
  }

  h2->Draw();
}


//-----------------------------------------------------------------------------
// plot 2D d^2sigma/dp/dcosth for pA-->pbarX (Striganov)
// such that one could integrate the cross section directly 
// for 1e8 events, the integral, divided by the number of events, is 1.51327e-04 mubarn
//-----------------------------------------------------------------------------
void pbar_xsec_norm_2(double PBeam, long int NEvents) {

  double pbeam   (PBeam);              // [GeV/c ], pbeam = 8.9 GeV corresponds to Ekin = 8 GeV
  //  double density = 19.3 ;	       // [g/cm^3], density of tungsten
  //  double length  = 16.;                // [cm    ], length of the production target
  //  double rho     = 0.3;                // [cm    ], radius of the production target
  //  double npot    = 3.6e20;             // [      ], Mu2e N(protons on target in 3 years)
  //  double A       = 183.84;             // [      ], atomic weight
  //  double xsec_Ta = 1.539e6;            // [mubarn], total inelastic pA xsection on Ta
  double mP      = 0.938272;           // [ GeV  ], proton mass
  double twopi   = 2*M_PI;

  //  float p2max = 2.0;

  TH2D* h_pbar_bgr = new TH2D("h_pbar_bgr","h_pbar_bgr: cross sec", 1000,0,5,200,-1,1);

  // double binx = h_pbar_bgr->GetXaxis()->GetBinWidth(1);
  // double biny = h_pbar_bgr->GetYaxis()->GetBinWidth(1);

  double   rn[2];
  TRandom3 trn3;
  
  TStntuple* stnt = TStntuple::Instance();

  for (int i=0; i<NEvents; i++) {
    trn3.RndmArray(2,rn);
    float  plab  = 5*rn[0];
    double costh = 2*rn[1]-1;
    float  e     = sqrt(plab*plab+mP*mP);
    float  thlab = TMath::ACos(costh);
//-----------------------------------------------------------------------------
// pa2pbarx_ returns E*d^3sigma/d^3p, integrate over dp*dcosth:
// w    = (E*d^3sigma/dp^3)*(1/E*dp^3)
// dp^3 = 2pi*p^2*dp*dcosth (integrated over phi)
//-----------------------------------------------------------------------------
    double cf = twopi*plab*plab/e;
    double w  = stnt->PBar_Striganov_Ed3SigdP3(pbeam,plab,thlab)*cf; 

//    printf("pbeam: %12.5f plab: %12.5f thlab: %12.5f cf: %12.5f w: %12.5e\n",pbeam,plab,thlab,cf,w);

    h_pbar_bgr->Fill(plab,costh,w);
  }

  h_pbar_bgr->Draw();
}
