// my implementation of the Feldman-Cousins algorithm
// 
#include "murat/alg/TLogLHR.hh"
#include "TCanvas.h"

ClassImp(TLogLHR)
//-----------------------------------------------------------------------------q
// CL > 0: CL has the meaning of a probability - 0.9, 0.95 .. .0.99 etc
//    < 0: CL is the "discovery probability", probability corresponding
//         to 5 gaussian sigma level, 
//-----------------------------------------------------------------------------
TLogLHR::TLogLHR(const char* Name, double CL, int DebugLevel):
  TNamed(Name,Name),
  fRn()
{
  SetCL(CL);
  
  fDebugLevel = DebugLevel;
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  fHist.fLHb    = new TH1D(Form("h_lhb_%s"     ,GetName()),"LHb"    ,MaxNx,-0.5,MaxNx-0.5);
  fHist.fLHb->SetMarkerStyle(20);
  fHist.fLHb->SetMarkerSize(1);
  fHist.fLHb->SetLineColor(kBlue);
  fHist.fLHb->SetMarkerColor(kBlue);
  
  fHist.fLHs    = new TH1D(Form("h_lhs_%s"     ,GetName()),"LHs"    ,MaxNx,-0.5,MaxNx-0.5);
  fHist.fLHs->SetMarkerStyle(20);
  fHist.fLHs->SetMarkerSize(1);
  fHist.fLHs->SetLineColor(kRed);
  fHist.fLHs->SetMarkerColor(kRed);

  fHist.fLogLHb = new TH1D(Form("h_log_lhb_%s" ,GetName()),"LogLHb" ,NxLogLH,-450,50);
  fHist.fLogLHb->SetMarkerStyle(20);
  fHist.fLogLHb->SetMarkerSize(1);
  fHist.fLogLHb->GetYaxis()->SetRangeUser(5.e-12,1);
  fHist.fLogLHb->SetLineColor(kBlue);
  fHist.fLogLHb->SetMarkerColor(kBlue);

  fHist.fLogLHs = new TH1D(Form("h_log_lhs_%s" ,GetName()),"LogLHs" ,NxLogLH,-450,50);
  fHist.fLogLHs->SetMarkerStyle(20);
  fHist.fLogLHs->SetMarkerSize(1);
  fHist.fLogLHs->GetYaxis()->SetRangeUser(5.e-12,1);
  fHist.fLogLHs->SetLineColor(kRed);
  fHist.fLogLHs->SetMarkerColor(kRed);
  
  fHist.fLogLHr = new TH1D(Form("h_log_lhr_%s" ,GetName()),"LogLHr" ,NxLogLH,-450,50);
  fHist.fLogLHr->SetMarkerStyle(20);
  fHist.fLogLHr->SetMarkerSize(1);
  fHist.fLogLHs->SetLineColor(kBlue+1);
  fHist.fLogLHs->SetMarkerColor(kBlue+1);
//-----------------------------------------------------------------------------
// calculate factorials, do that only once
// assume MaxNx to be large enough, so having N! values up to MaxNx-1 included is enough
//-----------------------------------------------------------------------------
  fFactorial[0] = 1; for (int i=1; i<MaxNx; i++) { fFactorial[i] = fFactorial[i-1]*i; }

  for (int i=0; i<MaxNx; i++) { fData[i] = new MData_t(); }
  
  fNExp = 10000000;
}

void  TLogLHR::SetCL(double CL) {
  double alpha = 1-TMath::Erf(5./sqrt(2));
  
  if (CL > 0) fCL = CL;
  else        fCL = 1-alpha/2; // always two-sided: 1-5.7330314e-07

  fLog1mCL        = log(1-fCL);
}


void TLogLHR::InitPoissonDist(double MuB, double MuS, double N, double* Prob, int NMax) {
  // 'N' is the number of measured events - it constrains the background fluctuations
  // N<0 means no prior knowledge 
  // array 'Prob' should have at least NMax elements
  // declare N as double to be able to scan
  // 'nature': given by MuB+Mus - sampled distribution

  if (fDebugLevel > 0) printf(">>> InitPoissonDist: MuB=%12.5e MuS=%12.5e N=%12.5e,NMax=%5i\n",MuB,MuS,N,NMax);

  double mean = MuB+MuS;
  Prob[0] = TMath::Exp(-mean);
  if (N < 0) {
    for (int i=1; i<NMax; i++) {
      Prob[i] = Prob[i-1]*mean/i;
    }
  }
  else {
//-----------------------------------------------------------------------------
// background probability constrained by the measurement of N events (Zech'1989)
//-----------------------------------------------------------------------------
    double pbn = 0; for (int k=0; k<=N; k++) { pbn += TMath::Exp(-MuB)*pow(MuB,k)/fFactorial[k]; }
    
    double pb[NMax];
    for (int i=0; i<NMax; i++) {
      if (i <= N) pb[i] = TMath::Exp(-MuB)*pow(MuB,i)/fFactorial[i]/pbn;
      else        pb[i] = 0;
    }
					// 'i' - bin in the constrained Poisson distribution
    for (int i=0; i<NMax; i++) {
      double pi = 0;
					// an experiment observed N events, use pb[k]
      for (int k=0; k<=i; k++) {
	double ps = TMath::Exp(-MuS)*pow(MuS,i-k)/fFactorial[i-k];
	pi        = pi + pb[k]*ps;
      }
      Prob[i] = pi;
    }
  }
}


void TLogLHR::Init(double MuB, double MuS, double NMeas, double MuBest) {
  // 'N' is the 'measured number of events'

  fMuB = MuB;
  fMuS = MuS;
  if (fDebugLevel > 0) printf("Init: MuB=%12.5e MuS: %12.5e NMeas:%12.5e\n",MuB,MuS,NMeas);

  if (fDebugLevel > 0) printf("Init: init fLHb\n");
  InitPoissonDist(MuB,MuBest,NMeas,fLHb,MaxNx);           // mubest normally 0, not zero for a scan

  if (fDebugLevel > 0) printf("Init: init fLHs - best idea of the signal\n");
  InitPoissonDist(MuB,MuS   ,NMeas,fLHs,MaxNx);
//-----------------------------------------------------------------------------
// [re]-initialize 1D histograms with the probabilities and integral probabilities
//-----------------------------------------------------------------------------
  fHist.fLHb->Reset();
  fHist.fLHs->Reset();

  fHist.fLogLHb->Reset();
  fHist.fLogLHs->Reset();
  fHist.fLogLHr->Reset();
  
  double s  = 0;
  double sw = 0;
  
  for (int i=0; i<MaxNx; i++) {
    fHist.fLHb->SetBinContent(i+1,fLHb[i]);
    fHist.fLHs->SetBinContent(i+1,fLHs[i]);

    fHist.fLHb->SetBinError(i+1,fLHb[i]/10.);
    fHist.fLHs->SetBinError(i+1,fLHs[i]/10.);

    double log_lhb(1.e-15);
    if (fLHb[i] > 0) log_lhb = log(fLHb[i]);
    double log_lhs(1.e-15);
    if (fLHs[i] > 0) log_lhs = log(fLHs[i]);

    fHist.fLogLHb->Fill(log_lhb,fLHs[i]);
    fHist.fLogLHs->Fill(log_lhs,fLHs[i]);

    MData_t* d = fData[i];

    d->N       = i;
    d->mub     = MuB;
    d->lhb     = fLHb[i];
    d->log_lhb = log_lhb;

    d->mus     = MuS;
    d->lhs     = fLHs[i];
    d->log_lhs = log_lhs;

    fHist.fLogLHr->Fill(d->log_lhr(),fLHs[i]);
//-----------------------------------------------------------------------------
// calculate the log_lhr mean using only the bins with non-zero probability of
// null hypothesis
//-----------------------------------------------------------------------------
    if (fLHb[i] != 0) {
      s   += fLHs[i]*d->log_lhr();
      sw  += fLHs[i];
    }
  }

  fMeanLLHR           = s/sw;
//-----------------------------------------------------------------------------
// order data in log_lhs
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

  if (fDebugLevel > 10) PrintData(fData);
}


void TLogLHR::PrintData(MData_t** Data, int MaxInd) {

  int maxind = MaxInd;
  if (maxind < 0) maxind = MaxNx;

  printf("------------------------------------------------------------------------------------------------------------\n");
  printf("  i   N      mub         lhb          log_lhb        mus          lhs        log_lhs       log_lhr       prob\n");
  printf("------------------------------------------------------------------------------------------------------------\n");

  double sum_lhs = 0;
  for (int i=0; i<maxind; i++) {
    MData_t* d = Data[i];
    sum_lhs += d->lhs;
    printf("%3i %3i %12.5e %12.5e  %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e \n",
	   i,d->N,d->mub,d->lhb,d->log_lhb, d->mus,d->lhs,d->log_lhs,d->log_lhr(),sum_lhs);
  }
}


double TLogLHR::PTail(MData_t** Data, double LogLHr) {
//-----------------------------------------------------------------------------
// calculate probability for an event to have log(LHR) <= LogLHr
// Data are ordered in increasing log_lhr , perform the computation starting from the end -
// this is computationally safer
//-----------------------------------------------------------------------------  
  double prob(1);
  
  for (int i=MaxNx-1; i>-0; i--) {
    MData_t* d = Data[i];

    if (fDebugLevel > 10) {
      printf(" PTail i:%4i lhs[i]: %12.5e log_lhr[i]: %12.5e LogLHr: %12.5e prob:%12.5e\n",
	     i,d->lhs,d->log_lhr(),LogLHr, prob);
    }

    if (d->log_lhr() <= LogLHr) break;
    prob -= d->lhs;
  }

  if (fDebugLevel > 0) printf(" PTail prob = %12.5e\n",prob);

  return prob;
}


void TLogLHR::ConfInterval(double Bgr, double SMin, double SMax, int NPoints, double* Prob, double* MeanLLHR) {
  // 'Bgr' : background expectation
  // 'N'   : number of measured events, has to be set to -1 here
  // scan SMin-SMax range of signals taking NSteps = NPoints-1 steps
  // to determine signal 'S' such that 'S+B' and 'B' hypotheses are inconsitent at 'fCL' level

  double step(0);

  if (NPoints > 1) step = (SMax-SMin)/(NPoints-1);
//-----------------------------------------------------------------------------
// measurement result, computation-wise, this place can be optimized 
//-----------------------------------------------------------------------------
  if (fDebugLevel > 0) printf(">>> ConfInterval   : Bgr=%12.5e SMin=%12.5e\n",Bgr,SMin);
 
  for (int is=0; is<NPoints; is++) {
    double sig = SMin+is*step;
//-----------------------------------------------------------------------------
// likelihood corresponding to the measurement of N events 
// we have an estimate of B, and know that background didn't fluctuate above N
//-----------------------------------------------------------------------------
    if (fDebugLevel > 10) printf("is = %3i, sig=%12.5e Bgr = %12.5e\n",is,sig,Bgr);
    
    Init(Bgr,sig,-1,0);
//-----------------------------------------------------------------------------
// ptail - probability, integrated over the region logLHr(B/B+sig) < log_lhr
//-----------------------------------------------------------------------------
    Prob    [is]   = PTail(fData,fLog1mCL);
    MeanLLHR[is]   = fMeanLLHR;
    
    if (fDebugLevel > 10) printf("is, fLog1mCL, prob, <logLHr>: %3i %12.5e %12.5e %12.5e\n",is,fLog1mCL,Prob[is],MeanLLHR[is]);
  }

  if (fDebugLevel > 0) printf(" >>> TLogLHR::ConfInterval : END\n");
}



 void TLogLHR::MeasInterval(double Bgr, double NMeas, double SMin, double SMax, int NPoints, double* Prob, double* MeanLLHR) {
  // 'Bgr' : background expectation
  // 'N'   : number of measured events
  // assume that a measurement has been performed, so the "best measured signal" s = N-Bgr;

  double step(0);

  double best_sig = NMeas - Bgr;
  if (best_sig < 0) best_sig = 0;

  if (NPoints > 1) step = (SMax-SMin)/(NPoints-1);
//-----------------------------------------------------------------------------
// measurement result, computation-wise, this place can be optimized 
//-----------------------------------------------------------------------------
  if (fDebugLevel > 0) printf(">>> MeasInterval   : Bgr=%12.5e NMeas = %12.5e SMin=%12.5e\n",Bgr,NMeas,SMin);
 
  for (int is=0; is<NPoints; is++) {
    double sig = SMin+is*step;
//-----------------------------------------------------------------------------
// likelihood corresponding to the measurement of N events 
// we have an estimate of B, and know that background didn't fluctuate above N
//-----------------------------------------------------------------------------
    if (fDebugLevel > 10) printf("is = %3i, sig=%12.5e Bgr = %12.5e NMeas = %12.5e\n",is,sig,Bgr,NMeas);
    
    Init(Bgr,best_sig,NMeas,sig);
//-----------------------------------------------------------------------------
// ptail - ratio of probabilites integrated over the region x < log(1-fCL)
//-----------------------------------------------------------------------------
    Prob[is]       = PTail(fData,fLog1mCL);
    MeanLLHR[is]   = fMeanLLHR;
    
    if (fDebugLevel > 10) printf("is, fLog1mCL, prob, <llhr>: %3i %12.5e %12.5e %12.5e\n",is,fLog1mCL,Prob[is],fMeanLLHR);
  }

  if (fDebugLevel > 0) printf(" >>> TLogLHR::ConfInterval : END\n");
}


void TLogLHR::DiscoveryProbCLb(double MuB, double SMin, double SMax, int NPoints, double* MuS, double* Prob) {
//-----------------------------------------------------------------------------
// in general, need to scan a range of signals, call this function multiple times
// calculate the probability for the background to be not consistent with the 'measurement'
//-----------------------------------------------------------------------------
  double lhb[MaxNx];

  double step = (NPoints > 1) ? (SMax-SMin)/(NPoints-1) : 0;
  
  TH1D* hist = new TH1D("hist","hist",1000,-90,10);

  for (int ix=0; ix<NPoints; ix++) {
    MuS[ix] = SMin+ix*step;
    double tot = MuB+MuS[ix];
//-----------------------------------------------------------------------------
// ndisc: number of "discovery experiments", pseudoexperiments in which NULL
// hypothesis is excluded at (1-fCL) level
//-----------------------------------------------------------------------------
    long int ndisc = 0;			
    hist->Reset();

    InitPoissonDist(MuB,0,-1,lhb,MaxNx);
    
    for (int i=0; i<fNExp; i++) {
      int nmeas = fRn.Poisson(tot);
//-----------------------------------------------------------------------------
// 'discovery experiment': likelihood of the background-only hypothesis < (1-CL)
//-----------------------------------------------------------------------------
      double p = lhb[nmeas];

      // double log_p = log(p);
      // hist->Fill(log_p);

      if (p < 1-fCL) ndisc += 1;
    }
    Prob[ix] = double(ndisc)/double(fNExp);
  }
}


void TLogLHR::DiscoveryProbCLs(double MuB, double SMin, double SMax, int NPoints, double* MuS, double* Prob) {
//-----------------------------------------------------------------------------
// in general, need to scan a range of signals, call this function multiple times
// just calculate the probability for background to be not conssitent with the measurement
//-----------------------------------------------------------------------------
  double step = (NPoints > 1) ? (SMax-SMin)/(NPoints-1) : 0;
  
  TH1D* hist = new TH1D("h_cls","log LHR cls",1000,-90,10);

  for (int ix=0; ix<NPoints; ix++) {
    MuS[ix] = SMin+ix*step;
    double tot = MuB+MuS[ix];
//-----------------------------------------------------------------------------
// ndisc: number of "discovery experiments", pseudoexperiments in which NULL
// hypothesis is excluded at (1-fCL) level
//-----------------------------------------------------------------------------
    long int ndisc = 0;			
    hist->Reset();

    Init(MuB,MuS[ix],-1,0);
    
    for (int i=0; i<fNExp; i++) {
      int nmeas = fRn.Poisson(tot);
//-----------------------------------------------------------------------------
// calculate likelihood of the background-only hypothesis, compare it to 
//-----------------------------------------------------------------------------
      double ps     = fLHs[nmeas];
      if (ps > 0) {
	double p   = fLHb[nmeas]/ps;
      	if (p < 1-fCL) ndisc += 1;
      }
    }
    Prob[ix] = double(ndisc)/double(fNExp);
  }
}

// void TLogLHR::DiscoveryProbCLs2(double MuB, double SMin, double SMax, int NPoints, double* MuS, double* Prob) {
// //-----------------------------------------------------------------------------
// // in general, need to scan a range of signals, call this function multiple times
// // just calculate the probability for background to be not conssitent with the measurement
// //-----------------------------------------------------------------------------
//   double step = (NPoints > 1) ? (SMax-SMin)/(NPoints-1) : 0;
  
//   TH1D* hist = new TH1D("h_cls","log LHR cls",1000,-90,10);

//   for (int ix=0; ix<NPoints; ix++) {
//     MuS[ix] = SMin+ix*step;
//     double tot = MuB+MuS[ix];
// //-----------------------------------------------------------------------------
// // ndisc: number of "discovery experiments", pseudoexperiments in which NULL
// // hypothesis is excluded at (1-fCL) level
// //-----------------------------------------------------------------------------
//     long int ndisc = 0;			
//     hist->Reset();

//     for (int i=0; i<fNExp; i++) {
//       int nmeas = fRn.Poisson(tot);
// //-----------------------------------------------------------------------------
// // calculate likelihood of the background-only hypothesis, compare it to 
// //-----------------------------------------------------------------------------
//       Init(MuB,MuS[ix],nmeas,0);
// //-----------------------------------------------------------------------------
// // calculate likelihood of the background-only hypothesis, compare it to 
// //-----------------------------------------------------------------------------
//       double ps     = fLHs[nmeas];
//       if (ps > 0) {
// 	double p   = fLHb[nmeas]/ps;
//       	if (p < 1-fCL) ndisc += 1;
//       }
//     }
//     Prob[ix] = double(ndisc)/double(fNExp);
//   }
// }
