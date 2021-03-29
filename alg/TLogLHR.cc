// my implementation of the Feldman-Cousins algorithm
// 
#include "murat/alg/TLogLHR.hh"
#include "TCanvas.h"

ClassImp(TLogLHR)
//-----------------------------------------------------------------------------
// CL > 0: CL has the meaning of a probability - 0.9, 0.95 .. .0.99 etc
//    < 0: CL is the "discovery probability", probability corresponding
//         to 5 gaussian sigma level, 
//-----------------------------------------------------------------------------
TLogLHR::TLogLHR(const char* Name, double CL, int DebugLevel):
  TNamed(Name,Name),
  fRn()
{
  if (CL > 0) fCL = CL;
  else        fCL = 1-TMath::Erf(5./sqrt(2)); // always two-sided: 5.7330314e-07

  fDebugLevel = DebugLevel;
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  fHist.fLHb    = new TH1D(Form("h_bg_prob_%s"   ,GetName()),"LHb"    ,MaxNx,-0.5,MaxNx-0.5);
  fHist.fLHs    = new TH1D(Form("h_bs_prob_%s"   ,GetName()),"LHs"    ,MaxNx,-0.5,MaxNx-0.5);

  fHist.fLogLHb = new TH1D(Form("h_log_lhb_%s"   ,GetName()),"LogLHb" ,NxLogLH,-180,120);
  fHist.fLogLHs = new TH1D(Form("h_log_lhs_%s"   ,GetName()),"LogLHs" ,NxLogLH,-180,120);
  fHist.fLogLHr = new TH1D(Form("h_log_lhr_%s"   ,GetName()),"LogLHr" ,NxLogLH,-180,120);
  fHist.fLogLHw = new TH1D(Form("h_log_lhw_%s"   ,GetName()),"LogLHw" ,NxLogLH,-180,120);

  fHist.fPTail  = new TH1D(Form("h_ptail_%s"     ,GetName()),"PTail"  ,NxLogLH,-180,120);
  
  for (int i=0; i<MaxNx; i++) {
    fHist.fLHPoi[i] = new TH1D(Form("h_lh_poi_%03i_%s",i,GetName()),Form("poisson prob for <m> = %i",i),
			       MaxNx,-0.5,MaxNx-0.5);
  }
  
  fNExp         = 10000000;
}

void TLogLHR::InitPoissonDist(double Mean, double* Prob, int N) {
  // array 'Prob' should have at least N elements
  
  Prob[0] = TMath::Exp(-Mean);
  for (int i=1; i<N; i++) {
    Prob[i] = Prob[i-1]*Mean/i;
  }
}


void TLogLHR::InitLogLHr(double* Num, double* Denom) {
  // weight is given by the numerator histogram

  fHist.fLogLHr->Reset();
  fHist.fLogLHw->Reset();

  for (int i=0; i<NxLogLH; i++)  fHist.fLogLHw->SetBinContent(i+1,1.e12);

  double x0 = fHist.fLogLHr->GetXaxis()->GetXmin();
  double bx = fHist.fLogLHr->GetBinWidth(1);
  
  for (int i=0; i<MaxNx; i++) {
    double llh_num = Num[i];
    double llh_den = Denom[i];
    double llhr = log(llh_num/llh_den);
   
    fHist.fLogLHr->Fill(llhr,llh_num);

    int bin = floor((llhr-x0)/bx)+1;
    fHist.fLogLHw->SetBinContent(bin,llhr);
  }
//-----------------------------------------------------------------------------
// integrate the tail
//-----------------------------------------------------------------------------
  fHist.fPTail->Reset();
  for (int i=0; i<MaxNx; i++) {
    double sum = 0;
    for (int i1=0; i1<i; i1++) {
      sum += fHist.fLogLHr->GetBinContent(i1+1);
    }
    fHist.fPTail->SetBinContent(i+1,sum);
  }
}



void TLogLHR::Init(double Bgr, double Sig) {
  //  fSigMean = SigMean;

  fMeanBgr = Bgr;
  fMeanSig = Sig;

  double poi[MaxNx];
  
  InitPoissonDist(fMeanBgr         , fLHb, MaxNx);
  InitPoissonDist(fMeanBgr+fMeanSig, fLHs, MaxNx);
//-----------------------------------------------------------------------------
// [re]-initialize 1D histograms with the probabilities and integral probabilities
//-----------------------------------------------------------------------------
  for (int i=0; i<MaxNx; i++) {
    fHist.fLHb->SetBinContent(i+1,fLHb[i]);
    fHist.fLHs->SetBinContent(i+1,fLHs[i]);

    double log_lhb = log(fLHb[i]);
    double log_lhs = log(fLHs[i]);

    fHist.fLogLHb->Fill(log_lhb,fLHs[i]);
    fHist.fLogLHs->Fill(log_lhs,fLHs[i]);

    InitPoissonDist(double(i), poi, MaxNx);
    for (int i1=0; i1<MaxNx; i1++) {
      fHist.fLHPoi[i]->SetBinContent(i1+1,poi[i1]);
    }
  }

  InitLogLHr(fLHs,fLHb);
}

void TLogLHR::PrintData(const char* Title, char DataType, void* Data, int MaxInd) {

  int k      = 0;
  int findex = 0;
  int n_per_line(20);

  double*   pd = (double*) Data;
  int*      pi = (int*   ) Data;
  
  printf("    ---------------------------- %s:",Title);
  for (int i=0; i<MaxInd; i++) {
    if (k % n_per_line == 0) {
      printf("\n%03i",n_per_line*findex);
      k = 0;
    }
    if      (DataType == 'd') printf(" %9.2e",pd[i]);
    else if (DataType == 'i') printf(" %9i"  ,pi[i]);
    k++;
  }

  if (k > 0) printf("\n");
}

//-----------------------------------------------------------------------------
// print results
//-----------------------------------------------------------------------------
void TLogLHR::PrintProbs(int N) {

  printf(" <bgr> = %10.4f <sig> = %10.4f\n",fMeanBgr,fMeanSig);
  printf(" IMin, IMax, Prob = %3i %3i %12.5f\n",fIMin,fIMax,fProb);

  PrintData("LHb"     ,'d',fLHb      ,N);
  PrintData("BestSig" ,'d',fBestSig  ,N);
  PrintData("LHs"     ,'d',fLHs      ,N);
  PrintData("BestProb",'d',fBestProb ,N);
  PrintData("LhRatio" ,'d',fLhRatio  ,N);
}


double TLogLHR::PTail(double LogLHr) {
  double prob(0);
  
  for (int i=0; i<NxLogLH; i++) {
    double log_lhr  = fHist.fLogLHw->GetBinContent(i+1);
					// skip uninitialized bins
    if (log_lhr < 1.e12) {
      if (log_lhr < LogLHr+1.e-6) {
	//      prob += fHist.fPTail->GetBinContent(i+1);
	prob += fHist.fLogLHr->GetBinContent(i+1);
      }
      else {
	break;
      }
    }
  }
  return prob;
}


void TLogLHR::ConfInterval(double Bgr, int N, double SMin, double SMax, int NSteps, double* Prob) {
  // 'Bgr' : background expectation
  // 'N'   : number of measured events
  // assume that the "measured" signal s = N-Bgr;

  double step(0);
  int    ns(1);

  printf(">>> TLogLHR::ConfInterval : Bgr = %12.5e N = %5i, SMin=%12.5e\n",Bgr,N,SMin);
    
  if (NSteps > 1) {
    step = (SMax-SMin)/NSteps;
    ns   = NSteps+1;
  }
//-----------------------------------------------------------------------------
// measurement result, computation-wise, this place can be optimized 
//-----------------------------------------------------------------------------
  double sbest  = N-Bgr;
  if (sbest < 0) sbest = 0;
  
  Init(Bgr,sbest);

  for (int is=0; is<ns; is++) {
    double sig = SMin+is*step;
    double tot = Bgr+sig;           // 'tested' mean for given signal assumption
//-----------------------------------------------------------------------------
// likelihood corresponding to the measurement of N events
//-----------------------------------------------------------------------------
    InitPoissonDist(tot,fLHb,MaxNx);
					// for debugging, dont' need otherwise
    for (int ib=0; ib<MaxNx; ib++) {
      fHist.fLHb->SetBinContent(ib+1,fLHb[ib]);
    }
    
    InitLogLHr     (fLHb,fLHs);
    
    double lh      = fLHb[N];
    double lh_best = fHist.fLHPoi[N]->GetBinContent(N);
    double log_lhr = log(lh) - log(lh_best);
    
    Prob[is]       = PTail(log_lhr);
    
    printf("is, lh, lh_best, log_lhr, prob: %3i %12.5e %12.5e %12.5e %12.5e\n",
	   is, lh, lh_best, log_lhr, Prob[is]);
  }
}
