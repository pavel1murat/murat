///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "TH1.h"
#include "TH2.h"
#include "TObjString.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TMath.h"

#include "math.h"
#include "string.h"
#include "limits/TLHChannel_ZZ_4l.hh"
#include "Stntuple/val/stntuple_val_functions.hh"
// #include "murat/ana/TAnaUtils.hh"

ClassImp(TLHChannel_ZZ_4l)

namespace {

  //                M(ZZ)  PT(ZZ)
  double x[8][2] = { {196.,  35.},
		     {190.,  30.},
		     {234.,  10.},
		     {192.,  27.},
					// high mass
		     {321.,  47.},
		     {315., 114.},
		     {325., 127.},
		     {334.,  44.}
  };

};

//-----------------------------------------------------------------------------
TLHChannel_ZZ_4l::TLHChannel_ZZ_4l(int HighMassOnly): TLHChannel("ZZ_4l","ZZ_4l") {

  fLhHist = new TH1F("ZZ_4l_Likelihood","ZZ 4l likelihood",400,-200,0);

  for (int i=0; i<8; i++) {
    fData[i][0] = x[i][0];
    fData[i][1] = x[i][1];
  }

  fHighMassOnly = HighMassOnly;
}

int TLHChannel_ZZ_4l::Init() {
//-----------------------------------------------------------------------------
// need to generate random numbers according to this histogram
//                    M     PT
//-----------------------------------------------------------------------------
  double  xbin_width, ybin_width, xmin, ymin, p;
  //  double  mass[100], prob[100];
  //  double  tot_prob;

  int     ix, iy;

  TLHChannel::Init();			// defines fProbHist

  xbin_width = fProbHist->GetXaxis()->GetBinWidth(1);  // all bins are the same
  ybin_width = fProbHist->GetYaxis()->GetBinWidth(1);  // all bins are the same

  xmin = fProbHist->GetXaxis()->GetXmin();
  ymin = fProbHist->GetYaxis()->GetXmin();

  fLhData = 0;

  int    i0;
  if (fHighMassOnly == 0) i0 = 0;
  else                    i0 = 4;

  for (int i=i0; i<8; i++) {
    ix       = (int) ((fData[i][0]-xmin)/xbin_width + 1);
    iy       = (int) ((fData[i][1]-ymin)/ybin_width + 1);
    p        = fProbHist->GetBinContent(ix,iy);
    fLhData += 2*TMath::Log(p);
  }

  printf(" >>> TLHChannel_ZZ_4l::fLhrData = %10.5f\n",fLhData);
  return 0;
}

//-----------------------------------------------------------------------------
TLHChannel_ZZ_4l::~TLHChannel_ZZ_4l() {
}


//-----------------------------------------------------------------------------
// return code = 1 : pseudoexperiment succeeded
//             = 0 : pseudoexperiment failed
//-----------------------------------------------------------------------------
int TLHChannel_ZZ_4l::PseudoExperiment(double& Likelihood) {

  TRandom3  r3;
  int       success, nhm, nlm, ix, iy, nn(-1), nevents, high_mass;
  double    mass[100], prob[100], msave[100], ptsave[100], psave[100];
  double    lh(-1.e10), m, pt, p, xbin_width, xmin, ymin, /*prob_hm,*/ tot_prob, ybin_width;
  
  nhm       = 0;
  nlm       = 0;
  tot_prob  = 1;
  //  prob_hm   = 1;

  xbin_width = fProbHist->GetXaxis()->GetBinWidth(1);  // all bins are the same
  ybin_width = fProbHist->GetYaxis()->GetBinWidth(1);  // all bins are the same

  xmin = fProbHist->GetXaxis()->GetXmin();
  ymin = fProbHist->GetYaxis()->GetXmin();

  nevents   = r3.Poisson(fNExpected);  // 5.8 is the mean.... for SM ZZ->4l in 6/fb

  Likelihood = 0;

  for (int ievent=0; ievent<nevents; ievent++) {

    ((TH2F*) fProbHist)->GetRandom2(m,pt);
    
    ix    = (int) ((m -xmin)/xbin_width + 1);
    iy    = (int) ((pt-ymin)/ybin_width + 1);
    p     = fProbHist->GetBinContent(ix,iy);

    if (m > 300) {
      mass[nhm] = m;
      prob[nhm] = p;
      nhm       = nhm+1;
      high_mass = 1;
    }
    else {
      nlm       = nlm+1;
      high_mass = 0;
    }

    if (fHighMassOnly == 0) {
					// consider all events
      lh = lh + 2*TMath::Log(p);
    }
    else {
					// consider only high-mass events
      if (high_mass) {
	lh = lh + 2*TMath::Log(p);
      }
    }
  }
//-----------------------------------------------------------------------------
// pseudoexperiment is generated
// there are 4 high mass events within 20 GeV from each other 
// and also 4 low mass events. Count low mass events because the total is 
// also higher than the expected total...
//-----------------------------------------------------------------------------
  if (nhm < 4) {
    lh = 2*TMath::Log(0.5);
                                                            goto FILL_HIST;
  }
//-----------------------------------------------------------------------------
// number of high-mass events is greater than 4
// if low-mass events are considered, check their number as well
//-----------------------------------------------------------------------------
  if (fHighMassOnly == 0) {
    if (nlm < 4) {
      lh = 2*TMath::Log(0.5);
                                                            goto FILL_HIST;
    }
  }
//-----------------------------------------------------------------------------
// check if masses of the high mass events are within 20 GeV from each other
// first order masses
//-----------------------------------------------------------------------------
  for (int i1=0; i1<nhm-1; i1++) {
    for (int i2=i1+1; i2<nhm; i2++) {
      if (mass[i1] > mass[i2]) {
	m        = mass[i1];
	mass[i1] = mass[i2];
	mass[i2] = m;
	
	p        = prob[i1];
	prob[i1] = prob[i2];
	prob[i2] = p;
      }
    }
  }
//-----------------------------------------------------------------------------
// nn : number of events within 20 GeV window
//-----------------------------------------------------------------------------
  nn = 0;
  for (int i1=0; i1<nhm-3; i1++) {
    nn      = 1;
    for (int i2=i1+1; i2<nhm; i2++) {
      if (mass[i2]-mass[i1] < 20) {
	nn++;
      }
      else {
	break;
      }
    }
  }

  if (nn < 4) {
//-----------------------------------------------------------------------------
// total number of high-mass events within 20 GeV is less than 4
//-----------------------------------------------------------------------------
      lh = 2*TMath::Log(0.05);
  }
					// fill the likelihood histogram
 FILL_HIST:;
  fLhHist->Fill(lh);
  Likelihood = lh;
//-----------------------------------------------------------------------------
// search mode: pseudoexperiment is considered a success if the calculated 
// likelihood is less then the likelihood of the data
//-----------------------------------------------------------------------------
  if (Likelihood > fLhData) success = 0;
  else {
    success = 1;
    printf(" TLHChannel_ZZ_4l: tot_prob = %12.5e lh  = %12.5e fLHData = %12.5e success = %i\n",
	   tot_prob, lh, fLhData, success);
    printf(" nevents, nlm, nhm, nn = %3i %3i %3i %3i\n",nevents,nlm,nhm,nn);

    printf("msave : ");
    for (int i=0; i<nevents; i++) printf("%12.5e ",msave[i]);
    printf("\n");
    printf("ptsave: ");
    for (int i=0; i<nevents; i++) printf("%12.5e ",ptsave[i]);
    printf("\n");
    printf("psave : ");
    for (int i=0; i<nevents; i++) printf("%12.5e ",psave [i]);
    printf("\n");
  }

  return success;
}



