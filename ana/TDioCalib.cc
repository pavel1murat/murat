//-----------------------------------------------------------------------------
// analysis code to implement detector momentum/energy calibrations 
// using DIO electrons at half field
//-----------------------------------------------------------------------------

#include "ana/TDioCalib.hh"

ClassImp(TDioCalib)

//-----------------------------------------------------------------------------
TDioCalib::TDioCalib(): TObject() {

  InitDetectorResponseFunction();
  fDioFileName = "hist/dio_half_field.dio_calib.hist";

  fHist.fDioMom  = gh1(fDioFileName.Data(),"DioCalib","gen_0/p");
  fHist.fRecoMom = gh1(fDioFileName.Data(),"DioCalib","trk_0/p2");

}

//-----------------------------------------------------------------------------
TDioCalib::~TDioCalib() {
}
//-----------------------------------------------------------------------------
// A1*exp((x-x0)^2/2sig^2 for (x-x0)/sig > -alp
// A1*A2*(B-(x-x0)/sig)^{-n}
// A2 = 
// P[0] = N
// P[1] = x[0]
// P[2] = sigma
// P[3] = alpha  (> 0)
// P[4] = n
//-----------------------------------------------------------------------------
// 2014-07-16: so far  - a placeholder for the response function
//-----------------------------------------------------------------------------
double TDioCalib::ReconstructedMomentum(double P0) {
  double p, dp;

  p  = P0;
  dp = fHist.fDetResponse->GetRandom();
  p  = P0+dp;

  return p;
}

//-----------------------------------------------------------------------------
// need to have it normalized to the unit area
//-----------------------------------------------------------------------------
double TDioCalib::f_crystal_ball(double* X, double* P) {
  double f;

  double mean, sigp, alpha, n;

  mean  = P[1];  // data come from the 50 MeV/c efficiency file
  sigp  = P[2];
  alpha = P[3];
  n     = P[4];

  double a, b, abs_alpha;

  double dx = (X[0]-mean)/sigp;
  abs_alpha = fabs(alpha);

  if (dx > -abs_alpha) {
    if (dx < 0){
      f = P[0]*TMath::Exp(-dx*dx/2.);
    }
    else {
      f = P[0]*((1-P[5])*TMath::Exp(-dx*dx/2.)+P[5]*TMath::Exp(-dx/P[6]));
    }
  }
  else {
    a = TMath::Power((n/abs_alpha),n)*TMath::Exp(-(alpha*alpha)/2);
    b = n/abs_alpha-abs_alpha;
    f = P[0]*a*TMath::Power((b-dx),-n);
  }

  return f;
}

//-----------------------------------------------------------------------------
// 1  anorm        1.37278e-01   2.48487e-01   2.21904e-04  -1.06923e-04
// 2  mean         4.93037e+01   1.31601e+00   1.44907e-03  -8.96999e-04
// 3  sigma        3.52361e-01   5.38599e-01  -3.31518e-03   8.09890e-05
// 4  alpha        5.88541e-01   3.05837e+00  -2.82805e-03  -1.21048e-04
// 5  N            2.93476e+00   1.68556e+01   3.27537e-03   2.07664e-05
// 6  tail_frac    1.93462e-02   3.50928e-01  -1.76909e-02   4.93193e-04
// 7  tail_lamb    4.53700e+00   2.30049e+02  -4.10784e-01  -2.05113e-05
//-----------------------------------------------------------------------------
void TDioCalib::InitDetectorResponseFunction() {
  double par[] = {
    1.37278e-01, 
    -0.7,   // -1.4,            // 4.93037e+01,
    3.52361e-01,
    5.88541e-01,
    2.93476e+00,
    1.93462e-02,
    4.53700e+00
  };

  fHist.fDetResponse = new TH1F("det_resp","response",200,-10,10);

  double x[2], wt;

  for (int i=1; i<200; i++) {
    x[0] = fHist.fDetResponse->GetBinCenter(i);
    wt   = f_crystal_ball(x,par);
    fHist.fDetResponse->SetBinContent(i,wt);
  }
}


//-----------------------------------------------------------------------------
// 2014-07-16: initial version (35-65 MeV/c);
// so far, use linear interpolation
//-----------------------------------------------------------------------------
double TDioCalib::TrkRecoEff(double P0) {

  int const np = 6;

  double pi[np] = { 35.,   40.,   45.,   50.,  55.,  60.};
  double ei[np] = { 24.,  493., 1386., 1678., 827., 366.};

  double eff (0);

  int    loc = 0;

  if ((P0 > pi[0]) && (P0 < pi[np-1])) {
    for (int i=0; i<np-1; i++) {
      if (P0 < pi[i+1]) {
	loc = i;
	break;
      }
    }

    eff = ei[loc]+(ei[loc+1]-ei[loc])/(pi[loc+1]-pi[loc])*(P0-pi[loc]);
  }

  return eff;
}


//-----------------------------------------------------------------------------
// estimate of the number of muons decayed in orbit 
// per year: 1.2*10e20*1.6*e-3*0.391
// 
// spectrum parameterization from Czarnecki et al, arXiv:1106.4756
//-----------------------------------------------------------------------------
double TDioCalib::f_dio_spectrum(double* X, double* P) {

  double a5(8.6434), a6(1.16874), a7(-1.87828e-2), a8(9.16327e-3);
  double emu(105.194), mAl(25133.);
  
  double de, de5, p;

  de  = emu-X[0]-X[0]*X[0]/(2*mAl);

  de5 = de*de*de*de*de;

  p   = P[0]*1e-17*de5*(a5 + de*(a6+de*(a7+a8*de)));

  if (p < 0) p = 0;

  return p;

}


//-----------------------------------------------------------------------------
// most of the numbers obrained with 'dev2' release, which is effectively 
// 'cvs co -D 2013-09-04 Offline'
//-----------------------------------------------------------------------------
void TDioCalib::HalfField(double Scale) {
  
  double  p, p_rec, wt, eff;

  fHist.fDioExpected = (TH1F*) fHist.fDioMom->Clone("h_dio_expected");

  fHist.fDioExpected->Reset();
//-----------------------------------------------------------------------------
// model detector response
//-----------------------------------------------------------------------------
  int nbins = fHist.fDioMom->GetNbinsX();

  for (int ib=1; ib<=nbins; ib++) {
					// p in MeV/c
    p = fHist.fDioMom->GetBinCenter(ib);
					// have a response function

    wt  = fHist.fDioMom->GetBinContent(ib);
					// reconstruction efficiency
    eff = TrkRecoEff(p);

    if ((p > 35.) && (p < 75.)) {
//-----------------------------------------------------------------------------
// consider response,
// h_dio_expected - the DIO momentum spectrum, corrected for the detector response
// it makes a template,  will need multiple templates for different momentum scales
// 1000 is just a number, fix it later
//-----------------------------------------------------------------------------
      for (int i=0; i<100000; i++) {
					// sample the "detector response" function,
					// involves the random number generation
	p_rec = ReconstructedMomentum(p)*Scale;
	fHist.fDioExpected->Fill(p_rec,wt*eff);
      }
    }
  }

  fHist.fRecoMom->SetMarkerStyle(24);
  fHist.fRecoMom->SetMarkerSize (0.7);
  fHist.fRecoMom->GetXaxis()->SetRangeUser(20,69.9999);
  fHist.fRecoMom->Draw("ep");

  double qn = fHist.fRecoMom->Integral()/fHist.fDioExpected->Integral();

  fHist.fDioExpected->Scale(qn);

  fHist.fDioExpected->SetFillColor(1);
  fHist.fDioExpected->SetFillStyle(3002);
  fHist.fDioExpected->Draw("same");

//-----------------------------------------------------------------------------
// now calculate chi2 in the region 40-55 MeV - bins 201-275
//-----------------------------------------------------------------------------
  double chi2(0), ydat, yf, err;
  int    n1   = 0;

  for (int ib=211; ib<=260; ib++) {
    ydat = fHist.fRecoMom->GetBinContent(ib);
    err  = fHist.fRecoMom->GetBinError(ib);
    yf   = fHist.fDioExpected->GetBinContent(ib);

    chi2 += (yf-ydat)*(yf-ydat)/(err*err+1.e-12);
    n1   += 1;
  }

  printf(" Scale = %10.5f, chi2/N = %12.5f/%5i \n",Scale,chi2, n1);

}
