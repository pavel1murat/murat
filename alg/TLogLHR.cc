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

  double alpha = 1-TMath::Erf(5./sqrt(2));
  
  if (CL > 0) fCL = CL;
  else        fCL = 1-alpha/2; // always two-sided: 1-5.7330314e-07

  fDebugLevel = DebugLevel;
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  fHist.fLHb    = new TH1D(Form("h_lhb_%s"     ,GetName()),"LHb"    ,MaxNx,-0.5,MaxNx-0.5);
  fHist.fLHs    = new TH1D(Form("h_lhs_%s"     ,GetName()),"LHs"    ,MaxNx,-0.5,MaxNx-0.5);
  fHist.fLHt    = new TH1D(Form("h_lht_%s"     ,GetName()),"LHt"    ,MaxNx,-0.5,MaxNx-0.5);

  fHist.fLogLHb = new TH1D(Form("h_log_lhb_%s" ,GetName()),"LogLHb" ,NxLogLH,-180,120);
  fHist.fLogLHb->SetMarkerStyle(20);
  fHist.fLogLHb->SetMarkerSize(1);
  fHist.fLogLHb->GetYaxis()->SetRangeUser(5.e-12,1);

  fHist.fLogLHs = new TH1D(Form("h_log_lhs_%s" ,GetName()),"LogLHs" ,NxLogLH,-180,120);
  fHist.fLogLHs->SetMarkerStyle(20);
  fHist.fLogLHs->SetMarkerSize(1);
  fHist.fLogLHs->GetYaxis()->SetRangeUser(5.e-12,1);
  
  fHist.fLogLHr = new TH1D(Form("h_log_lhr_%s" ,GetName()),"LogLHr" ,NxLogLH,-180,120);
  fHist.fLogLHw = new TH1D(Form("h_log_lhw_%s" ,GetName()),"LogLHw" ,NxLogLH,-180,120);
  fHist.fLogLHn = new TH1D(Form("h_log_lhn_%s" ,GetName()),"LogLHn" ,NxLogLH,-180,120);

  fHist.fPTail  = new TH1D(Form("h_ptail_%s"   ,GetName()),"PTail"  ,NxLogLH,-180,120);
  
  // calculate factorials, do that only once, assume MaxNx to be larre enough,
  // so having up calculation done up to MaxNx-1 included is enough
  fFactorial[0] = 1;
  for (int i=1; i<MaxNx; i++) fFactorial[i] = fFactorial[i-1]*i;

  for (int i=0; i<MaxNx; i++) fData[i] = new MData_t();
  
  fNExp         = 10000000;
}

void TLogLHR::InitPoissonDist(double Bgr, double Sig, int N, double* Prob, MData_t** Data, int NMax) {
  // N<0 means no prior knowledge 
  // array 'Prob' should have at least NMax elements
  // 'N' is the number of measured events - it constrains the background fluctuations
  // and truncates the distribution

  printf(">>> InitPoissonDist: Bgr=%12.5e Sig=%12.5e N=%3i,NMax=%5i\n",Bgr,Sig,N,NMax);

  double mean = Bgr+Sig;
  Prob[0] = TMath::Exp(-mean);
  if (N < 0) {
    for (int i=1; i<NMax; i++) {
      Prob[i] = Prob[i-1]*mean/i;
      if (Prob[i] < 1.e-12) Prob[i] = 1.e-12;
    }
  }
  else {
    // if N >= 0, we know that in the performed measurement, background is not higher than N
    // so the measurement constrains the probability calculation
    for (int i=1; i<NMax; i++) {
      Prob[i] = 0;
      if (i <= N) Prob[i] = Prob[0]*pow(mean,i)/fFactorial[i];
    }
  }
}


void TLogLHR::Init(double MuB, double MuS, int N, MData_t** Data) {
  // 'N' is the 'measured number of events'

  fMuB = MuB;
  fMuS = MuS;
  printf("INit: MuB=%12.5e MuS: %12.5e N:%3i\n",MuB,MuS,N);
  
  InitPoissonDist(MuB,  0, N, fLHb, fData, MaxNx);
  InitPoissonDist(MuB,MuS, N, fLHs, fData, MaxNx);
//-----------------------------------------------------------------------------
// [re]-initialize 1D histograms with the probabilities and integral probabilities
//-----------------------------------------------------------------------------
  for (int i=0; i<MaxNx; i++) {
    fHist.fLHb->SetBinContent(i+1,fLHb[i]);
    fHist.fLHs->SetBinContent(i+1,fLHs[i]);

    double log_lhb(1.e-12); if (fLHb[i] > 0) log_lhb = log(fLHb[i]);
    double log_lhs(1.e-12); if (fLHs[i] > 0) log_lhs = log(fLHs[i]);

    fHist.fLogLHb->Fill(log_lhb,fLHs[i]);
    fHist.fLogLHs->Fill(log_lhs,fLHs[i]);

    fData[i]->N       = i;
    fData[i]->mub     = MuB;
    fData[i]->lhb     = fLHb[i];
    fData[i]->log_lhb = log_lhb;

    fData[i]->mus     = MuS;
    fData[i]->lhs     = fLHs[i];
    fData[i]->log_lhs = log_lhs;
  }

//-----------------------------------------------------------------------------
// order in log_lhs
//-----------------------------------------------------------------------------
  for (int k1=0; k1<MaxNx-1; k1++) {
    double x1 = fData[k1]->log_lhr();
    for (int k2=k1+1; k2<MaxNx; k2++) {
      double x2 = fData[k2]->log_lhr();
      if (x2 < x1) {
	MData_t* x = fData[k1];
	fData[k1]  = fData[k2];
	fData[k2]  = x;
      }
    }
  }

  PrintData(fData);
}


void TLogLHR::PrintData(MData_t** Data, int MaxInd) {

  int maxind = MaxInd;
  if (maxind < 0) maxind = MaxNx;

  printf("----------------------------------------------------------------------------------- \n");
  printf("  i   N      mub         lhb         log_lhb        mus          lhs        log_lhs\n");
  printf("----------------------------------------------------------------------------------- \n");
  for (int i=0; i<maxind; i++) {
    printf("%3i %3i %12.5e %12.5e  %12.5e %12.5e %12.5e %12.5e\n",
	   i,Data[i]->N,
	   Data[i]->mub,Data[i]->lhb,Data[i]->log_lhb,
	   Data[i]->mus,Data[i]->lhs,Data[i]->log_lhs);
  }
}


double TLogLHR::PTail(MData_t** Data, double LogLHr) {
  // Data are ordered in increasing log_lhs
  double prob(0);
  
  for (int i=0; i<MaxNx; i++) {
    printf(" PTail i:%4i lhs[i]: %12.5e log_lhr[i]: %12.5e LogLHr: %12.5e prob:%12.5e\n",
	   i,Data[i]->lhs,Data[i]->log_lhr(),LogLHr, prob);

    if (Data[i]->log_lhr() > LogLHr) break;
    prob += Data[i]->lhs;
  }

  printf(" PTail prob = %12.5e\n",prob);

  return prob;
}


void TLogLHR::ConfInterval(double Bgr, int N, double SMin, double SMax, int NSteps, double* Prob) {
  // 'Bgr' : background expectation
  // 'N'   : number of measured events
  // assume that the "measured" signal s = N-Bgr;

  double step(0);

  if (NSteps > 1) step = (SMax-SMin)/(NSteps-1);
//-----------------------------------------------------------------------------
// measurement result, computation-wise, this place can be optimized 
//-----------------------------------------------------------------------------
  printf(">>> ConfInterval   : Bgr=%12.5e N = %5i SMin=%12.5e\n",Bgr,N,SMin);
 
  for (int is=0; is<NSteps; is++) {
    double sig = SMin+is*step;
//-----------------------------------------------------------------------------
// likelihood corresponding to the measurement of N events 
// we have an estimate of B, and know that background didn't fluctuate above N
//-----------------------------------------------------------------------------
    printf("is = %3i, sig=%12.5e Bgr = %12.5e N = %3i\n",is,sig,Bgr,N);
    
    Init(Bgr,sig,-1,fData);

    //    double log_lhs = fData[N]->log_lhs;          // - log(lh_best);
    double log_lhr = log(1-fCL);
//-----------------------------------------------------------------------------
// ptail - ratio of probabilites integrated over the region x < log_lhr
//-----------------------------------------------------------------------------
    Prob[is]       = PTail(fData,log_lhr);
    
    printf("is, log_lhr, prob: %3i %12.5e %12.5e\n",
	   is, log_lhr, Prob[is]);
  }

  printf(" >>> TLogLHR::ConfInterval : END\n");
}
