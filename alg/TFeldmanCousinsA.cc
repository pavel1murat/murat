// my implementation of Feldman-Cousins algorithm
// 

#include "murat/alg/TFeldmanCousinsA.hh"

#include "TCanvas.h"

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

  fCumBsProb[0] = fBsProb[0];
  for (int i=1; i<MaxNx;i++) {
    fCumBsProb[i] = fCumBsProb[i-1]+fBsProb[i];
  }
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

  PrintData("BgProb"  ,'d',fBgProb   ,N);
  PrintData("BestSig" ,'d',fBestSig  ,N);
  PrintData("BsProb"  ,'d',fBsProb   ,N);
  PrintData("CuBsProb",'d',fCumBsProb,N);
  PrintData("BestProb",'d',fBestProb ,N);
  PrintData("LhRatio" ,'d',fLhRatio  ,N);
  PrintData("Rank   " ,'i',fRank     ,N);
}

//-----------------------------------------------------------------------------
int TFeldmanCousinsA::ConstructBelt(double Bgr, double SMin, double SMax, int NSteps) {

  double step = (SMax-SMin)/NSteps;

  for (int i=0; i<NSteps; i++) {
    double s = SMin+(i+0.5)*step;
    ConstructInterval(Bgr,s);
    for (int ix=0; ix<MaxNx; ix++) {
      if ((ix < fIMin) || (ix > fIMax)) fBelt[i][ix] = 0;
      else                              fBelt[i][ix] = 1;
    }
  }
  return 0;
}


//-----------------------------------------------------------------------------
// constructing the FC belt assumes dividing the signal interval into NSteps,
// make that a parameter
// calculation of the median for low expected backgrounds doesn't have much sense,
// so the mean is a preference
//-----------------------------------------------------------------------------
double TFeldmanCousinsA::UpperLimit(double Bgr, double SMin, double SMax, int NSteps) {

  ConstructBelt(Bgr, SMin, SMax, NSteps);

  double step = (SMax-SMin)/NSteps;
  
  TH2D* h_belt = new TH2D("h2","FC belt",
			  TFeldmanCousinsA::MaxNx,0,TFeldmanCousinsA::MaxNx,
			  NSteps,SMin-step/2,SMax-step/2);
  
  for (int iy=0; iy<NSteps; iy++) {
    for (int ix=0; ix<TFeldmanCousinsA::MaxNx; ix++) {
      h_belt->SetBinContent(ix+1,iy+1,fBelt[iy][ix]);
    }
  }

  TCanvas* c_belt = new TCanvas("c_belt","plot h_belt",1000,800);
  c_belt->cd();

  h_belt->Draw("box");

  double sbelt[2][TFeldmanCousinsA::MaxNx];
  
  for (int ix=0; ix<TFeldmanCousinsA::MaxNx; ix++) {
    int iymin     = 0;
    int iymax     = 0;
    int inside    = 0;
    for (int iy=0; iy<NSteps; iy++) {
      if (fBelt[iy][ix] > 0) {
	if (inside == 0) {
	  iymin     = iy;
	  inside    = 1;
	}
      }
      else {
	if (inside == 1) {
	  iymax     = iy;
	  inside    = 0;
	  break;
	}
      }
    }

    sbelt[0][ix] = (iymin+0.5)*step;
    sbelt[1][ix] = (iymax+0.5)*step;

    printf("ix, smin, smax : %3i %12.5f %12.5f\n",ix,sbelt[0][ix],sbelt[1][ix]);
  }
//-----------------------------------------------------------------------------
// now generate bgr-only pseudo-experiments and plot the distribution
// for the excluded signal strength
//-----------------------------------------------------------------------------
  TH1D* h1 = new TH1D("h1","excluded",2000,0,20);
  TRandom3 trn;
  for (int i=0; i<1000000; i++) {
    int rn = trn.Poisson(Bgr);
    if (rn < TFeldmanCousinsA::MaxNx) {
      h1->Fill(sbelt[1][rn]);
    }
    else {
      printf("trouble, Bill Robertson: rn = %3i\n",rn);
    }
  }

  TCanvas* c2 = new TCanvas("c2","plot h1",1000,800);
  c2->cd();
  h1->Draw();

  double fc_upper_limit = h1->GetMean();
  printf(" mean excluded value : %12.5f\n",fc_upper_limit);
  return fc_upper_limit;

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
