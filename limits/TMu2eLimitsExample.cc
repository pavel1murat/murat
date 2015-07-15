///////////////////////////////////////////////////////////////////////////////
// learn to use Tom's MCLIMIT
// needs libmclimit.so to be loaded, this is why it is a separate subpackage
// 
///////////////////////////////////////////////////////////////////////////////
#include "TH1F.h"
#include "math.h"
#include "string.h"
#include "limits/TMu2eLimitsExample.hh"


ClassImp(TMu2eLimitsExample)

//-----------------------------------------------------------------------------
TMu2eLimitsExample::TMu2eLimitsExample() {
					// 1. create a new mclimit_csm object
  fMcLimit   = new mclimit_csm();

  fNullHyp   = 0;
  fNullHypPe = 0;

  fTestHyp   = 0;
  fTestHypPe = 0;

  fBgrHist  = 0;
  fSigHist  = 0;

  fBgrLevel = 0.1;

  fNPseudoExp = 10000;


//-----------------------------------------------------------------------------
// 2. define a histogram to fit
//-----------------------------------------------------------------------------
  fXMin  =   0;
  fXMax  = 500;
  fNBins = 100;

  fDataHist = new TH1F("hist","M(ZZ)",fNBins,fXMin,fXMax);
  fBgrHist  = new TH1F("bgr" ,"bgr"  ,fNBins,fXMin,fXMax);
  fSigHist  = new TH1F("sig" ,"sig"  ,fNBins,fXMin,fXMax);

  Init();
}



//-----------------------------------------------------------------------------
TMu2eLimitsExample::~TMu2eLimitsExample() {

  if (fBgrHist ) delete fBgrHist;
  if (fSigHist ) delete fSigHist;

  if (fNullHyp  ) delete fNullHyp;
  if (fNullHypPe) delete fNullHypPe;

  if (fTestHyp  ) delete fTestHyp;
  if (fTestHypPe) delete fTestHypPe;

  delete fMcLimit;
}


//-----------------------------------------------------------------------------
int TMu2eLimitsExample::Init() {

  for (int i=1; i<=fNBins; i++) {
    fDataHist->SetBinContent(i,fBgrLevel);
  }

  fMcLimit->set_datahist(fDataHist,"data_hist_name");
//-----------------------------------------------------------------------------
// null hypothesis
//-----------------------------------------------------------------------------
  fNullHyp = new csm_model();

  for (int i=1; i<=fNBins; i++) {
    fBgrHist->SetBinContent(i,fBgrLevel);
    fBgrHist->SetBinError(i,fBgrLevel/10.);
  }

  const char* pname[2] = {"b1","b2"};

  //  const int np_bgr = 2;
  double   nps_lo[2]  = { -0.01, -0.01};
  double   nps_hi[2]  = {  0.01,  0.01};

  TH1*   bgr_lo_shape[2];
  TH1*   bgr_hi_shape[2];

  double  bgr_lo_sigma[2], bgr_hi_sigma[2];

  for (int i=0; i<2; i++) {
    bgr_lo_shape[i] = NULL;
    bgr_hi_shape[i] = NULL;
    bgr_lo_sigma[i] = 1.;
    bgr_hi_sigma[1] = 1.;
  }

  fNullHyp->add_template(fBgrHist,
			 1.0,
			 0,
			 pname,
			 nps_lo,
			 nps_hi,
			 bgr_lo_shape,
			 bgr_lo_sigma,
			 bgr_hi_shape,
			 bgr_hi_sigma,
			 0,                   // no bin-by-bin errors 
			 0,
			 "mass");

  fNullHypPe = (csm_model*) fNullHyp->Clone();

  fMcLimit->set_null_hypothesis(fNullHyp);
  fMcLimit->set_null_hypothesis_pe(fNullHypPe);
//-----------------------------------------------------------------------------
// 3. test hypothesis: background + signal in one bin
//-----------------------------------------------------------------------------
  fTestHyp = new csm_model();

				        // background template
  fTestHyp->add_template(fBgrHist,
			 1.0,           // normalization scale factor
			 0,             // no nuissance parameters
			 pname,
			 nps_lo,
			 nps_hi,
			 bgr_lo_shape,
			 bgr_lo_sigma,
			 bgr_hi_shape,
			 bgr_hi_sigma,
			 0,                 // no bin-by bin errors ...
			 0,                 // background is fixed
			 "mass");

					// signal template
  double x, w;
  double x0 = 301;

  for (int i=1; i<=fNBins; i++) {
    x = fSigHist->GetBinCenter(i);
    w = fSigHist->GetBinWidth(i);
    //    printf("i,x,w,x0,fabs(x-x0) = %3i %10.3f %10.3f %10.3f %10.3f\n",i,x,w,x0,fabs(x-x0));
    if (fabs(x-x0) <= w/2.) {
      printf("------------ i,x,x0 = %3i %10.3f %10.3f\n",i,x,x0);
      fSigHist->SetBinContent(i,1.);
      fSigHist->SetBinError(i,0.001);
    }
  }

  const char*   names[2] = {"P0", "P1"};
  double lo_side_err[2] = {-0.01, -0.01};
  double hi_side_err[2] = { 0.01,  0.01};

  double scale_factor = 1.;
  int n_var_par = 0;

  TH1* lo_shape_hist[100];
  TH1* hi_shape_hist[100];

  for (int i=0; i<100; i++) {
    lo_shape_hist[i] = NULL;
    hi_shape_hist[i] = NULL;
  }
  
  double sig_lo_shape[100], sig_hi_shape[100];

  fTestHyp->add_template(fSigHist,
			 scale_factor,
			 n_var_par,
			 names,
			 lo_side_err,
			 hi_side_err,
			 lo_shape_hist,
			 sig_lo_shape,
			 hi_shape_hist,
			 sig_hi_shape,
			 0,                     // no errors , just scaling
			 1,                     // signal has to be scaled
			 "mass");
  
  fTestHyp->print();

  fMcLimit->set_test_hypothesis(fTestHyp);

  fTestHypPe = (csm_model*) fTestHyp->Clone();

  fMcLimit->set_test_hypothesis_pe(fTestHypPe);

  return 0;
}


//-----------------------------------------------------------------------------
int TMu2eLimitsExample::FindLimit(int Npe) {

  if (Npe > 0) fMcLimit->set_npe(Npe);
  else         fMcLimit->set_npe(fNPseudoExp);

  double chi2 = 1.; // fMcLimit->chisquared();
  double s95  = fMcLimit->s95();

  printf("chi2 = %10.3f, s95=%10.5f\n",chi2,s95);
	
  return 0;
}
