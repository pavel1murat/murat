//
#include "murat/plot/murat_plot_functions.hh"

extern "C" void set_p2max_(float* p2max);
//-----------------------------------------------------------------------------
// plot 2D d^2sigma/dp/dcosth for pA-->pbarX (Striganov)
// such that one could integrate the cross section directly 
// for 1e8 events, the integral, divided by the number of events, is 1.51327e-04 mubarn
//-----------------------------------------------------------------------------
void pbar_bgr(long int NEvents = 1000000) {
  
  float  pbeam   = 8.9 ;	       // [GeV/c ], corresponds to Ekin = 8 GeV
  double density = 19.3 ;	       // [g/cm^3], density of tungsten
  double length  = 16.;                // [cm    ], length of the production target
  double rho     = 0.3;                // [cm    ], radius of the production target
  double npot    = 3.6e20;             // [      ], Mu2e N(protons on target in 3 years)
  double A       = 183.84;             // [      ], atomic weight
  double xsec_Ta = 1.539e6;            // [mubarn], total inelastic pA xsection on Ta
  double mP      = 0.93825;            // [ GeV  ], proton mass
  double twopi   = 2*M_PI;

  float p2max = 2.0;

  TH2D* h2 = new TH2D("h2","cross sec", 1000,0,5,200,-1,1);

  double binx = h2->GetXaxis()->GetBinWidth(1);
  double biny = h2->GetYaxis()->GetBinWidth(1);

  double   rn[2];
  TRandom3 trn3;
  
  set_p2max_(&p2max);

  for (int i=0; i<NEvents; i++) {
    trn3.RndmArray(2,rn);
    float  plab  = 5*rn[0];
    double costh = 2*rn[1]-1;
    float  e     = sqrt(plab*plab+mP*mP);
    float  thlab = TMath::ACos(costh);
    double sinth = sqrt(1-costh*costh);
//-----------------------------------------------------------------------------
// pa2pbarx_ returns E*d^3sigma/d^3p, integrate over dp*dcosth:
// w    = (E*d^3sigma/dp^3)*(1/E*dp^3)
// dp^3 = 2pi*p^2*dp*dcosth (integrated over phi)
//-----------------------------------------------------------------------------
    double cf = twopi*plab*plab/e;
    double w  = TStntuple::PBar_Striganov_Ed3SigdP3(pbeam,plab,thlab)*cf; 

    h2->Fill(plab,costh,w);
  }

  h2->Draw();

  double xs_integral = h2->Integral()*binx*biny/NEvents;
  
  printf(" >>> Cross section integral: %12.5e\n", xs_integral);
}
