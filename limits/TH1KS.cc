//

#include "murat/limits/TH1KS.hh"
#include "TMath.h"

ClassImp(TH1KS)

//-----------------------------------------------------------------------------
TH1KS::TH1KS():TH1F() {
}

//-----------------------------------------------------------------------------
TH1KS::TH1KS(const char* Name, const char* Title, int NBins, float XLow, float XHigh)
  :TH1F(Name,Title,NBins,XLow,XHigh) {
}

//-----------------------------------------------------------------------------
TH1KS::~TH1KS() {
}

//______________________________________________________________________________
Double_t TH1KS::KolmogorovTest(const TH1 *h2, int NPExp, Option_t *option) const
{
   //  Statistical test of compatibility in shape between
   //  THIS histogram and h2, using Kolmogorov test.
   //
   //     Default: Ignore under- and overflow bins in comparison
   //
   //     option is a character string to specify options
   //         "U" include Underflows in test  (also for 2-dim)
   //         "O" include Overflows     (also valid for 2-dim)
   //         "N" include comparison of normalizations
   //         "D" Put out a line of "Debug" printout
   //         "M" Return the Maximum Kolmogorov distance instead of prob
   //         "X" Run the pseudo experiments post-processor with the following procedure:
   //             make pseudoexperiments based on random values from the parent
   //             distribution and compare the KS distance of the pseudoexperiment
   //             to the parent distribution. Bin the KS distances in a histogram,
   //             and then take the integral of all the KS values above the value
   //             obtained from the original data to Monte Carlo distribution.
   //             The number of pseudo-experiments nEXPT is currently fixed at 1000.
   //             The function returns the integral.
   //             (thanks to Ben Kilminster to submit this procedure). Note that
   //             this option "X" is much slower.
   //
   //   The returned function value is the probability of test
   //       (much less than one means NOT compatible)
   //
   //  Code adapted by Rene Brun from original HBOOK routine HDIFF
   //
   //  NOTE1
   //  A good description of the Kolmogorov test can be seen at:
   //    http://www.itl.nist.gov/div898/handbook/eda/section3/eda35g.htm
   //
   //  NOTE2
   //  see also alternative function TH1::Chi2Test
   //  The Kolmogorov test is assumed to give better results than Chi2Test
   //  in case of histograms with low statistics.
   //
   //  NOTE3 (Jan Conrad, Fred James)
   //  "The returned value PROB is calculated such that it will be
   //  uniformly distributed between zero and one for compatible histograms,
   //  provided the data are not binned (or the number of bins is very large
   //  compared with the number of events). Users who have access to unbinned
   //  data and wish exact confidence levels should therefore not put their data
   //  into histograms, but should call directly TMath::KolmogorovTest. On
   //  the other hand, since TH1 is a convenient way of collecting data and
   //  saving space, this function has been provided. However, the values of
   //  PROB for binned data will be shifted slightly higher than expected,
   //  depending on the effects of the binning. For example, when comparing two
   //  uniform distributions of 500 events in 100 bins, the values of PROB,
   //  instead of being exactly uniformly distributed between zero and one, have
   //  a mean value of about 0.56. We can apply a useful
   //  rule: As long as the bin width is small compared with any significant
   //  physical effect (for example the experimental resolution) then the binning
   //  cannot have an important effect. Therefore, we believe that for all
   //  practical purposes, the probability value PROB is calculated correctly
   //  provided the user is aware that:
   //     1. The value of PROB should not be expected to have exactly the correct
   //  distribution for binned data.
   //     2. The user is responsible for seeing to it that the bin widths are
   //  small compared with any physical phenomena of interest.
   //     3. The effect of binning (if any) is always to make the value of PROB
   //  slightly too big. That is, setting an acceptance criterion of (PROB>0.05
   //  will assure that at most 5% of truly compatible histograms are rejected,
   //  and usually somewhat less."

   TString opt = option;
   opt.ToUpper();

   Double_t prob = 0;
   TH1 *h1 = (TH1*)this;
   if (h2 == 0) return 0;
   TAxis *axis1 = h1->GetXaxis();
   TAxis *axis2 = h2->GetXaxis();
   Int_t ncx1   = axis1->GetNbins();
   Int_t ncx2   = axis2->GetNbins();

   // Check consistency of dimensions
   if (h1->GetDimension() != 1 || h2->GetDimension() != 1) {
      Error("KolmogorovTest","Histograms must be 1-D\n");
      return 0;
   }

   // Check consistency in number of channels
   if (ncx1 != ncx2) {
      Error("KolmogorovTest","Number of channels is different, %d and %d\n",ncx1,ncx2);
      return 0;
   }

   // Check consistency in channel edges
   Double_t difprec = 1e-5;
   Double_t diff1 = TMath::Abs(axis1->GetXmin() - axis2->GetXmin());
   Double_t diff2 = TMath::Abs(axis1->GetXmax() - axis2->GetXmax());
   if (diff1 > difprec || diff2 > difprec) {
      Error("KolmogorovTest","histograms with different binning");
      return 0;
   }

   Bool_t afunc1 = kFALSE;
   Bool_t afunc2 = kFALSE;
   Double_t sum1 = 0, sum2 = 0;
   Double_t ew1, ew2, w1 = 0, w2 = 0;
   Int_t bin;
   Int_t ifirst = 1;
   Int_t ilast = ncx1;
   // integral of all bins (use underflow/overflow if option)
   if (opt.Contains("U")) ifirst = 0;
   if (opt.Contains("O")) ilast = ncx1 +1;
   for (bin = ifirst; bin <= ilast; bin++) {
      sum1 += h1->GetBinContent(bin);
      sum2 += h2->GetBinContent(bin);
      ew1   = h1->GetBinError(bin);
      ew2   = h2->GetBinError(bin);
      w1   += ew1*ew1;
      w2   += ew2*ew2;
   }
   if (sum1 == 0) {
      Error("KolmogorovTest","Histogram1 %s integral is zero\n",h1->GetName());
      return 0;
   }
   if (sum2 == 0) {
      Error("KolmogorovTest","Histogram2 %s integral is zero\n",h2->GetName());
      return 0;
   }

   // calculate the effective entries.
   // the case when errors are zero (w1 == 0 or w2 ==0) are equivalent to
   // compare to a function. In that case the rescaling is done only on sqrt(esum2) or sqrt(esum1)
   Double_t esum1 = 0, esum2 = 0;
   if (w1 > 0)
      esum1 = sum1 * sum1 / w1;
   else
      afunc1 = kTRUE;    // use later for calculating z

   if (w2 > 0)
      esum2 = sum2 * sum2 / w2;
   else
      afunc2 = kTRUE;    // use later for calculating z

   if (afunc2 && afunc1) {
      Error("KolmogorovTest","Errors are zero for both histograms\n");
      return 0;
   }


   Double_t s1 = 1/sum1;
   Double_t s2 = 1/sum2;

   // Find largest difference for Kolmogorov Test
   Double_t dfmax =0, rsum1 = 0, rsum2 = 0;

   for (bin=ifirst;bin<=ilast;bin++) {
      rsum1 += s1*h1->GetBinContent(bin);
      rsum2 += s2*h2->GetBinContent(bin);
      dfmax = TMath::Max(dfmax,TMath::Abs(rsum1-rsum2));
   }

   // Get Kolmogorov probability
   Double_t z, prb1=0, prb2=0, prb3=0;

   // case h1 is exact (has zero errors)
  if  (afunc1)
      z = dfmax*TMath::Sqrt(esum2);
  // case h2 has zero errors
  else if (afunc2)
      z = dfmax*TMath::Sqrt(esum1);
  else
     // for comparison between two data sets
     z = dfmax*TMath::Sqrt(esum1*esum2/(esum1+esum2));

   prob = TMath::KolmogorovProb(z);

   // option N to combine normalization makes sense if both afunc1 and afunc2 are false
   if (opt.Contains("N") && !(afunc1 || afunc2 ) ) {
      // Combine probabilities for shape and normalization,
      prb1 = prob;
      Double_t d12    = esum1-esum2;
      Double_t chi2   = d12*d12/(esum1+esum2);
      prb2 = TMath::Prob(chi2,1);
      // see Eadie et al., section 11.6.2
      if (prob > 0 && prb2 > 0) prob *= prb2*(1-TMath::Log(prob*prb2));
      else                      prob = 0;
   }
   // X option. Pseudo-experiments post-processor to determine KS probability
   //   const Int_t nEXPT = 1000;
   if (opt.Contains("X") && !(afunc1 || afunc2 ) ) {
      Double_t dSEXPT;
      //      TH1 *hExpt = (TH1*)(gDirectory ? gDirectory->CloneObject(this,kFALSE) : gROOT->CloneObject(this,kFALSE));

      TH1 *hExpt =  (TH1*) this->Clone();
      // make nEXPT experiments (this should be a parameter)
      prb3 = 0;
      for (Int_t i=0; i < NPExp; i++) {
         hExpt->Reset();
         hExpt->FillRandom(h1,(Int_t)esum2);
         dSEXPT = KolmogorovTest(hExpt,NPExp,"M");
         if (dSEXPT>dfmax) prb3 += 1.0;
      }
      prb3 /= (Double_t)NPExp;
      delete hExpt;
   }

   // debug printout
   if (opt.Contains("D")) {
      printf(" Kolmo Prob  h1 = %s, sum bin content =%g  effective entries =%g\n",h1->GetName(),sum1,esum1);
      printf(" Kolmo Prob  h2 = %s, sum bin content =%g  effective entries =%g\n",h2->GetName(),sum2,esum2);
      printf(" Kolmo Prob     = %g, Max Dist = %g\n",prob,dfmax);
      if (opt.Contains("N"))
         printf(" Kolmo Prob     = %f for shape alone, =%f for normalisation alone\n",prb1,prb2);
      if (opt.Contains("X"))
         printf(" Kolmo Prob     = %f with %d pseudo-experiments\n",prb3,NPExp);
   }
   // This numerical error condition should never occur:
   if (TMath::Abs(rsum1-1) > 0.002) Warning("KolmogorovTest","Numerical problems with h1=%s\n",h1->GetName());
   if (TMath::Abs(rsum2-1) > 0.002) Warning("KolmogorovTest","Numerical problems with h2=%s\n",h2->GetName());

   if(opt.Contains("M"))      return dfmax;
   else if(opt.Contains("X")) return prb3;
   else                       return prob;
}
