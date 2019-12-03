// my implementation of Feldman-Cousins algorithm
// 

#include "murat/alg/TFeldmanCousinsA.hh"

ClassImp(TFeldmanCousinsA)

//-----------------------------------------------------------------------------
// CL > 0: CL has the meaning of a probability - 0.9, 0.95 .. .0.99 etc
//    < 0: CL is the "discovery probability", probability corresponding
//         to the 5 gaussian sigma level
//-----------------------------------------------------------------------------
TFeldmanCousinsA::TFeldmanCousinsA(const char* Name, double CL, int DebugLevel):
  TNamed(Name,Name),
  fRn()
{
  if (CL > 0) fCL = CL;
  else        fCL = 1-5.733e-7;

  fDebugLevel = DebugLevel;
//-----------------------------------------------------------------------------
// zero probabilities
//-----------------------------------------------------------------------------
  // for (int i=0; i<MaxNx; i++) {
  //   fProb[i] = 0.;
  // }
//-----------------------------------------------------------------------------
// calculate factorials - just once
//-----------------------------------------------------------------------------
  fFactorial[0] = 1;
  for (int i=1; i<MaxNx; i++) {
    fFactorial[i] = fFactorial[i-1]*i;
  }
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  fBgProbHist    = new TH1D(Form("h_bg_prob_%s",GetName()),"h_bg_prob",MaxNx,0,MaxNx);
  fBsProbHist    = new TH1D(Form("h_bs_prob_%s",GetName()),"h_bs_prob",MaxNx,0,MaxNx);

  fProbHist      = new TH1D(Form("h_prob_2D_%s",GetName()),"h prob 2D",MaxNx,-0.5,MaxNx-0.5);
  fBeltHist      = new TH1D(Form("h_belt_2D_%s",GetName()),"h belt 2D",MaxNx,-0.5,MaxNx-0.5);

  fNPE           = 10000000;
}

//-----------------------------------------------------------------------------
// the length of 'Prob' should be at least N
//-----------------------------------------------------------------------------
void TFeldmanCousinsA::InitPoissonDist(double Mean, double* Prob, int N) {

  for (int i=0; i<N; i++) {
    Prob[i] = TMath::Power(Mean,i)*TMath::Exp(-Mean)/fFactorial[i];
  }
}

//-----------------------------------------------------------------------------
void TFeldmanCousinsA::Init(double Bgr, double Sig) {
  //  fSigMean = SigMean;

  fMeanBgr = Bgr;
  fMeanSig = Sig;
  
  InitPoissonDist(fMeanBgr         , fBgProb, MaxNx);
  InitPoissonDist(fMeanBgr+fMeanSig, fBsProb, MaxNx);
//-----------------------------------------------------------------------------
// re-initialize 1D histograms
//-----------------------------------------------------------------------------
  for (int i=0; i<MaxNx; i++) {
    fBgProbHist->SetBinContent(i+1,fBgProb[i]);
    fBsProbHist->SetBinContent(i+1,fBsProb[i]);
  }
}

//-----------------------------------------------------------------------------
// construct Feldman-Cousins belt
//-----------------------------------------------------------------------------
int TFeldmanCousinsA::ConstructInterval(double Bgr, double Sig) {
  int rc(0);				// return code

  Init(Bgr,Sig);
  
  for (int ix=0; ix<MaxNx; ix++) {
    double sbest = ix-fMeanBgr;

    if (sbest <= 0) sbest = 0;

    double sb = sbest+fMeanBgr;

    fBestProb[ix] = TMath::Power(sb,ix)*TMath::Exp(-sb)/fFactorial[ix];

    fBestSig[ix]  = sbest;
    fLhRatio[ix]  = fBsProb[ix]/fBestProb[ix];
  }
//-----------------------------------------------------------------------------
// sort ranks
//-----------------------------------------------------------------------------
  double rmax;
    
  for (int i=0; i<MaxNx; i++) {
    fRank[i] = i;
  }

  for (int i1=0; i1<MaxNx; i1++) {
    rmax = fLhRatio[fRank[i1]];
    for (int i2=i1+1; i2<MaxNx; i2++) {
      double r2 = fLhRatio[fRank[i2]];
      if (r2 > rmax) {
	rmax      = r2;
	int i     = fRank[i1];
	fRank[i1] = fRank[i2];
	fRank[i2] = i;
      }
    }
  }
//-----------------------------------------------------------------------------
// build confidence interval corresponding to the probability fCL - 0.9 , 0.95, etc
//-----------------------------------------------------------------------------
  fIMin = MaxNx;
  fIMax = -1;

  int covered = 0;
  
  fProb = 0;
  
  for (int ix=0; ix<MaxNx; ix++) {
    int ind = fRank[ix];
    fProb    += fBsProb[ind];
    if (ind < fIMin) fIMin = ind;
    if (ind > fIMax) fIMax = ind;

    if(fDebugLevel > 0) {
      printf(" ix ind=fRank[ix] fBsProb[ind] fProb, 1-fProb fIMin fIMax :%3i %3i %10.3e %10.3e %10.3e %3i %3i\n",
	     ix,ind,fBsProb[ind],fProb,1-fProb,fIMin,fIMax);
    }
    
    if (fProb > fCL) {
      covered = 1;
      break;
    }
  }

  if (covered == 0) {
    printf("trouble ! interval not covered : prob = %12.5e , 1-CL = %12.5e\n",fProb,1-fCL);
  }
//-----------------------------------------------------------------------------
// suppose, everything is OK, update the belt histogram
//-----------------------------------------------------------------------------
  for (int i=fIMin; i<=fIMax; i++) {
    fBeltHist->SetBinContent(i+1,1);
  }
//-----------------------------------------------------------------------------
// finally, build the probability histogram
//-----------------------------------------------------------------------------
  for (int i=0; i<MaxNx; i++) {
    fProbHist->SetBinContent(i+1,fBsProb[i]);
  }

  return rc;
}

//-----------------------------------------------------------------------------
void TFeldmanCousinsA::PrintData(const char* Title, char DataType, void* Data, int MaxInd) {

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
void TFeldmanCousinsA::PrintProbs(int N) {

  printf(" <bgr> = %10.4f <sig> = %10.4f\n",fMeanBgr,fMeanSig);
  printf(" IMin, IMax, Prob = %3i %3i %12.5f\n",fIMin,fIMax,fProb);

  PrintData("BgProb"  ,'d',fBgProb  ,N);
  PrintData("BestSig" ,'d',fBestSig ,N);
  PrintData("BsProb"  ,'d',fBsProb  ,N);
  PrintData("BestProb",'d',fBestProb,N);
  PrintData("LhRatio" ,'d',fLhRatio ,N);
  PrintData("Rank   " ,'i',fRank    ,N);
}

//-----------------------------------------------------------------------------
// in general, need to scan a range of signals, call this function multiple times
//-----------------------------------------------------------------------------
void TFeldmanCousinsA::DiscoveryProb(double Bgr, double Signal) {

//-----------------------------------------------------------------------------
// construct FC CL confidence interval (covering fCL) assuming no signal
//-----------------------------------------------------------------------------
  ConstructInterval(Bgr,0);

  double tot = Bgr+Signal;
  long int ndisc = 0;
  for (int i=0; i<fNPE; i++) {
    int rn = fRn.Poisson(tot);
//-----------------------------------------------------------------------------
// definition of discovery: rn > fIMax, i.e  the  probability to observe'rn' is
// less than 1-fCL
//-----------------------------------------------------------------------------
    if (rn > fIMax) ndisc ++;
  }

  double prob = double(ndisc)/double(fNPE);

  printf (">>> ndisc: %10li disc prob = %12.5e\n",ndisc,prob);
}
