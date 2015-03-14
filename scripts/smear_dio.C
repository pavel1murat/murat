//-----------------------------------------------------------------------------
// smeared DIO spectrum in Al
// this file includes 2 routines:
// - smear_dio_spectrum()
// - save_dio_hist()
//
// root [17] smear_dio_spectrum()
// n_stopped_muons =   5.7600e+17
// n_dio =    1.194e+17
// n_conv =      6.363
//  crystal_ball integral =   6.3012e+00
// ymin, N(conv), N(DIO), N(DIO smeared) =   103.0000      5.418      0.966      1.156
// ymin, N(conv), N(DIO), N(DIO smeared) =   103.1000      5.307      0.529      0.704
// ymin, N(conv), N(DIO), N(DIO smeared) =   103.2000      5.178      0.271      0.414
// ymin, N(conv), N(DIO), N(DIO smeared) =   103.3000      5.026      0.128      0.234
// ymin, N(conv), N(DIO), N(DIO smeared) =   103.4000      4.846      0.055      0.127
// ymin, N(conv), N(DIO), N(DIO smeared) =   103.5000      4.631      0.020      0.065
// ymin, N(conv), N(DIO), N(DIO smeared) =   103.6000      4.370      0.006      0.032
// ymin, N(conv), N(DIO), N(DIO smeared) =   103.7000      4.051      0.001      0.015
// ymin, N(conv), N(DIO), N(DIO smeared) =   103.8000      3.653      0.000      0.006
// ymin, N(conv), N(DIO), N(DIO smeared) =   103.9000      3.149      0.000      0.002
//-----------------------------------------------------------------------------
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TMath.h"
#include "include/Stntuple/val/stntuple_val_functions.hh"


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
namespace {
  TGraphErrors* gr;
  TH1D*         h_conv;
  TH1D*         h_dio;
  TH1D*         h_dio_2;
  TH1D*         h_dio_smeared;

  int    kNBins  = 500 ;
  double kEhmin  = 101.;
  double kEhmax  = 106.;



};

//-----------------------------------------------------------------------------
double f_eff(double* X, double* P) {
  double f, dx;

  dx = (X[0]-P[1])/P[2];

  if (dx < 0) {
    f = 0;
  }
  else {
    f = P[0]*(1-TMath::Exp(-dx));
  }

  return f;
}

//-----------------------------------------------------------------------------
// the efficiency numbers used in the fit are the numbers of the reconstructed 
// "Set C" tracks, no momentum window cuts
// 1.4/2. factor accounts for generation in cos(theta) in [-0.6,0.8] window
//
// also, for conversion electrons get 4889 "Set C" tracks, use this number to 
// introduce a scale factor (sf = 4889/5591) to renormalize the DIO efficiency
//
// root [7] .L smear_dio.C
// root [35] eff_vs_mom()
//  FCN=0.599931 FROM MIGRAD    STATUS=CONVERGED     106 CALLS         107 TOTAL
//                      EDM=2.54207e-09    STRATEGY= 1      ERROR MATRIX ACCURATE 
//   EXT PARAMETER                                   STEP         FIRST   
//   NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
//    1  p0           3.78838e-01   3.84800e-02   8.02894e-06   5.18293e-03
//    2  p1           8.51550e+01   7.90913e-01   2.50166e-04  -7.00714e-05
//    3  p2           7.70675e+00   2.53132e+00   4.00162e-04  -8.11730e-05
//
// the fit result: eff = 0.378838*(1-exp(-(p-85.1550)/7.70675));
// 
// 
//------------------------------------------------------------------------------
void  eff_vs_mom() {

  int const np = 6;

  double dio_sf = 4889./5591.;
					// momenta
  double p[np]  = {  80., 87.,     90.,   95.,  100.,  105.};
					// numbers of Set C tracks 
  double qn[np] = {   7., 1368., 2784., 4459., 5473., 5591.};

  double ex[np], y[np], ey[np];
					// total number of generated events  - 10000,
					// -0.6 < cos_the < 0.8 gives a factor of 1.4/2.0 = 0.7 
  for (int i=0; i<np; i++) {
    y [i] = qn[i]/(10000*2./1.4)*dio_sf;
    ex[i] = 0;
    ey[i] = 0.02;
  }

  gr = new TGraphErrors(np,p,y,ex,ey);

  TF1* fit_eff = new TF1("fit_eff",f_eff,94,106,3);

  fit_eff->SetParameter(0,0.40);
  fit_eff->SetParameter(1,80.);
  fit_eff->SetParameter(2,10.);

  gr->Fit("fit_eff");
  gr->Draw("ALP");
}


//-----------------------------------------------------------------------------
// need to have it normalized to the unit area
//-----------------------------------------------------------------------------
double f_crystal_ball(double* X, double* P) {
  double f;
  double dp(0);
  //  double sigp(0.288), alpha(0.698), n(2.694);

  // choose mean between 100 and 105...

  // previous value
  //  double sigp(0.230), alpha(0.600), n(3.0);

  // v2: parameterize v4_1_1 trk_1/p

  double sigp(0.2345), alpha(0.3267), n(9.449);

  double a, b, abs_alpha;

  double dx = (X[0]+dp-P[1])/sigp;
  abs_alpha = fabs(alpha);

  if (dx > -abs_alpha) {
    f = P[0]*1e-17*TMath::Exp(-dx*dx/2.);
  }
  else {
    a = TMath::Power((n/abs_alpha),n)*TMath::Exp(-(alpha*alpha)/2);
    b = n/abs_alpha-abs_alpha;
    f = P[0]*1.e-17*a*TMath::Power((b-dx),-n);
  }

  return f;
}

//-----------------------------------------------------------------------------
// track reconstruction efficiency parameterization: fit results - see above
// the fit result: eff = 0.378838*(1-exp(-(p-85.1550)/7.70675));
//-----------------------------------------------------------------------------
double trk_reco_eff(double P) {

  double dx = (P-85.1550)/7.70675;
  if (dx < 0) dx = 0;

  double e = 0.378838*(1-exp(-dx));

  return e;
}

//-----------------------------------------------------------------------------
// estimate of the number of muons decayed in orbit 
// per year: 1.2*10e20*1.6*e-3*0.391
// 
// spectrum parameterization from Czarnecki et al, arXiv:1106.4756
//-----------------------------------------------------------------------------
double f_dio_spectrum(double* X, double* P) {

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
TF1  *f_conv(0), *f_dio(0);

void smear_dio_spectrum() {

  double mu_capture_prob(0.609);
					// roughly speaking, fraction of captures within the 
					// active window [700,1700] ns  (0.53 = 2594./4889.)
  double live_fraction(0.53);
					// efficiency 

  double pmax(104.973), pconv(104.973), eff_conv(0.4889*14./20.);

					// total number of 8 GeV protons on target in 3 years

  double n_prot      (1.2e20 * 3.);

					// assumed, for illustration purposes, B(mu->e) conversion 
  double br_mu2e_conv(1.e-16);
					// muon (transport+stopping) efficiency : 1.6e-3

  double n_stopped_muons = n_prot*1.6e-3;

  printf("n_stopped_muons = %12.4e\n",n_stopped_muons);
					// track reco eff for DIO's is accounted for lated 
  double n_dio = n_stopped_muons*(1-mu_capture_prob)*live_fraction;

  printf("n_dio = %12.3e\n",n_dio);

  double n_conv = n_stopped_muons*mu_capture_prob*br_mu2e_conv*live_fraction*eff_conv;

  printf("n_conv = %10.3f\n",n_conv);

  double emin =  90.;
  double emax = 110.;

  if (f_dio) delete f_dio;
  f_dio = new TF1("dio_0",f_dio_spectrum,emin,emax,1);
  f_dio->SetParameter(0,n_dio);


  if (f_conv) delete f_conv;
  f_conv = new TF1("conv_0",f_crystal_ball,emin,emax,2);

  f_conv->SetParameter(0,1.);                // normalization
  f_conv->SetParameter(1,pconv);

  double anorm = n_conv/f_conv->Integral(emin,emax);

  f_conv->SetParameter(0,anorm);

  TH1::SetDefaultSumw2(1);

  h_conv = new TH1D("h_conv","conv"         ,kNBins,kEhmin,kEhmax);
  h_dio  = new TH1D("h_dio" ,"DIO electrons",kNBins,kEhmin,kEhmax);

  double kEhmin2 = kEhmin-2;
  int kNBins2    = kNBins+200;

  h_dio_2 = new TH1D("h_dio_2" ,"DIO electrons 2",kNBins2,kEhmin2,kEhmax);

  double y1, y2, bin, p;

  bin = (kEhmax-kEhmin)/kNBins;
//-----------------------------------------------------------------------------
// acceptance*reconstruction efficiency for conversion electrons is already accounted for
// but not efficiency for the DIO's
// account for mean particle energy losses
//-----------------------------------------------------------------------------
  double mean_eloss(0.913);

  for (int i=0; i<kNBins; i++) {
    p = kEhmin+(i+0.5)*bin;
					// here we account for the mean energy losses, assumed to be constant 

    y1 = f_conv->Eval(p+mean_eloss)*bin;
    h_conv->SetBinContent(i+1,y1);

    y2 = f_dio->Eval(p+mean_eloss)*bin*trk_reco_eff(p);
    h_dio->SetBinContent(i+1,y2);
  }
//-----------------------------------------------------------------------------
// fill h_dio_2 histogram - the one with EMIN = emin-1 - use it to generate the smeared
// histogram
//-----------------------------------------------------------------------------
  for (int i=0; i<kNBins2; i++) {
    p = kEhmin2+(i+0.5)*bin;
					// here we account for the mean energy losses, assumed to be constant 

    y2 = f_dio->Eval(p+mean_eloss)*bin*trk_reco_eff(p);
    h_dio_2->SetBinContent(i+1,y2);
  }

  printf(" crystal_ball integral = %12.4e\n",f_conv->Integral(kEhmin,kEhmax));

  // f_conv->Draw();

  // f_dio->Draw("same");


  TCanvas* c = new_slide("c_smear_dio","results",2,2,1300,800);
  TPad* p1 = (TPad*) c->GetPrimitive("p1");
//-----------------------------------------------------------------------------
// draw conversion histogram and unsmeared DIO
//-----------------------------------------------------------------------------
  p1->cd(1);
  h_conv->SetMaximum(0.1);
  h_conv->Draw();

  h_dio->Draw("same");
//-----------------------------------------------------------------------------
// now - the smearing function
//-----------------------------------------------------------------------------
  p1->cd(2);

  TF1* f_smear = new TF1("f_smear",f_crystal_ball,-10,10,2);
  
  f_smear->SetParameter(0,1.);                // normalization
  f_smear->SetParameter(1,0);

  f_smear->Draw();


  double bin_sm = (10.-(-10.))/2000.;
  TH1D* h_smear = new TH1D("h_smear"," smearing", 2000,-10,10);

  for (int i=0; i<2000; i++) {
    p = -10+bin_sm*(i+0.5);
    y1 = f_smear->Eval(p);
    h_smear->SetBinContent(i+1,y1);
  }

  h_smear->Draw("same");
//-----------------------------------------------------------------------------
// check random numbers 
//-----------------------------------------------------------------------------
  p1->cd(3);
  TH1D* h_rnd_smear = new TH1D("h_rnd_smear"," RND smearing", 2000,-10,10);

  for (int i=0; i<100000; i++) {
    p = f_smear->GetRandom();
    h_rnd_smear->Fill(p);
  }

  h_rnd_smear->Draw();
//-----------------------------------------------------------------------------
// original and smeared DIO histograms
//-----------------------------------------------------------------------------
  p1->cd(4);

  h_dio_smeared = new TH1D("h_dio_smeared","DIO smeared",kNBins,kEhmin,kEhmax);

  double  p0, weight, dp;
  int nev_smear(100000);

  for (int i=0; i<kNBins2; i++) {
					// here we account for the mean energy losses, assumed to be constant 
    p0     = h_dio_2->GetBinCenter(i+1);
    weight = h_dio_2->GetBinContent(i+1)/nev_smear;

					// smearing - produce weighted histogram 
    for (int j=0;  j<nev_smear; j++) {
      dp = f_smear->GetRandom();

      p = p0+dp;

      h_dio_smeared->Fill(p,weight);
    }
  }

  h_dio->Draw();
  h_dio_smeared->SetLineColor(2);
  h_dio_smeared->Draw("same");
  h_conv->Draw("same");

  for (int i=201; i<300; i+=10) {
    double ymin     = h_dio_smeared->GetBinLowEdge(i);
    double qn_conv  = h_conv->Integral(i,kNBins);
    double qn_dio_1 = h_dio->Integral(i,kNBins);
    double qn_dio_2 = h_dio_smeared->Integral(i,kNBins);

    printf("ymin, N(conv), N(DIO), N(DIO smeared) = %10.4f %10.3f %10.3f %10.3f\n",
	   ymin, qn_conv, qn_dio_1, qn_dio_2);
  }
}

//-----------------------------------------------------------------------------
void save_dio_hist() {

  TFile* f;

//-----------------------------------------------------------------------------
// conversions
//-----------------------------------------------------------------------------
  f = TFile::Open("conv01.mu2e_limits.hist","recreate");

  TDirectory* d_ana = f->mkdir("Ana");
  f->cd("Ana");

  TDirectory* d_mu2e = d_ana->mkdir("Mu2eLimits");
  d_mu2e->cd();

  TDirectory* d_hist = d_mu2e->mkdir("Hist");
  d_hist->cd();

  TDirectory* d_trk = d_hist->mkdir("trk_1");
  d_trk->cd();
  TH1D* h_conv_p = (TH1D*) h_conv->Clone("p");
  //  h_conv_p->Write();

  //  d_trk->Write();
  f->Write();
  f->Close();
//-----------------------------------------------------------------------------
// now - DIO histograms
//-----------------------------------------------------------------------------
  f = TFile::Open("dio01.mu2e_limits.hist","recreate");

  d_ana = f->mkdir("Ana");
  f->cd("Ana");

  d_mu2e = d_ana->mkdir("Mu2eLimits");
  d_mu2e->cd();

  d_hist = d_mu2e->mkdir("Hist");
  d_hist->cd();

  d_trk = d_hist->mkdir("trk_1");
  d_trk->cd();

  TH1D* h_dio_p = (TH1D*) h_dio->Clone("p");
  //  h_dio_p->Write();

  //  d_trk->Write();
  f->Write();
  f->Close();
//-----------------------------------------------------------------------------
// now - DIO smeared histograms
//-----------------------------------------------------------------------------
  f = TFile::Open("dio02.mu2e_limits.hist","recreate");

  d_ana = f->mkdir("Ana");
  f->cd("Ana");

  d_mu2e = d_ana->mkdir("Mu2eLimits");
  d_mu2e->cd();

  d_hist = d_mu2e->mkdir("Hist");
  d_hist->cd();

  d_trk = d_hist->mkdir("trk_1");
  d_trk->cd();

  TH1D* h_dio_smeared_p = (TH1D*) h_dio_smeared->Clone("p");
  //  h_dio_p->Write();

  //  d_trk->Write();
  f->Write();
  f->Close();

//-----------------------------------------------------------------------------
// cosmics + radiative pion capture + antiprotons - 0.2 events/ 1.5 MeV ;
// set uncertainty to 10% of the value
//-----------------------------------------------------------------------------
  f = TFile::Open("cosm01.mu2e_limits.hist","recreate");

  d_ana = f->mkdir("Ana");
  f->cd("Ana");

  d_mu2e = d_ana->mkdir("Mu2eLimits");
  d_mu2e->cd();

  d_hist = d_mu2e->mkdir("Hist");
  d_hist->cd();

  d_trk = d_hist->mkdir("trk_1");
  d_trk->cd();

  TH1D* h_cosm = (TH1D*) h_dio_smeared->Clone("p");
  h_cosm->Reset();

  double cosmic_bgr_level = 0.2/(kNBins*1.5/(kEhmax-kEhmin));

  int nb = h_cosm->GetNbinsX();
  for (int i=0; i<nb; i++) {
    h_cosm->SetBinContent(i+1,cosmic_bgr_level);
    h_cosm->SetBinError  (i+1,cosmic_bgr_level/10.);
  }

  f->Write();
  f->Close();
//-----------------------------------------------------------------------------
// finally, empty mock data histogram
//-----------------------------------------------------------------------------
  f = TFile::Open("data01.mu2e_limits.hist","recreate");

  d_ana = f->mkdir("Ana");
  f->cd("Ana");

  d_mu2e = d_ana->mkdir("Mu2eLimits");
  d_mu2e->cd();

  d_hist = d_mu2e->mkdir("Hist");
  d_hist->cd();

  d_trk = d_hist->mkdir("trk_1");
  d_trk->cd();

  TH1D* h_data = (TH1D*) h_dio_smeared->Clone("p");
  h_data->Reset();

  f->Write();
  f->Close();

  delete f;
}
