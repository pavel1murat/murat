// my implementation of the Feldman-Cousins algorithm
// 
#include "murat/alg/TFeldmanCousinsA.hh"
#include "TCanvas.h"

ClassImp(TFeldmanCousinsA)
//-----------------------------------------------------------------------------
// CL > 0: CL has the meaning of a probability - 0.9, 0.95 .. .0.99 etc
//    < 0: CL is the "discovery probability", probability corresponding
//         to 5 gaussian sigma level, 
//-----------------------------------------------------------------------------
TFeldmanCousinsA::TFeldmanCousinsA(const char* Name, double CL, int DebugLevel):
  TNamed(Name,Name),
  fRn()
{
  if (CL > 0) fCL = CL;
  else        fCL = 1-5.733e-7/2.; // always one-sided

  fDebugLevel = DebugLevel;
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
  fHist.fBgProb    = new TH1D(Form("h_bg_prob_%s"   ,GetName()),"h_bg_prob",MaxNx,0,MaxNx);
  fHist.fBsProb    = new TH1D(Form("h_bs_prob_%s"   ,GetName()),"h_bs_prob",MaxNx,0,MaxNx);
  fHist.fCumBgProb = new TH1D(Form("h_cumbg_prob_%s",GetName()),"h_cumbg_prob",MaxNx,0,MaxNx);
  fHist.fCumBsProb = new TH1D(Form("h_cumbs_prob_%s",GetName()),"h_cumbs_prob",MaxNx,0,MaxNx);

  fHist.fProb      = new TH1D(Form("h_prob_2D_%s"    ,GetName()),"h prob 2D" ,MaxNx,-0.5,MaxNx-0.5);
  fHist.fBelt      = nullptr;

  fNExp            = 10000000;
}

//-----------------------------------------------------------------------------
// the length of 'Prob' should be at least N
//-----------------------------------------------------------------------------
void TFeldmanCousinsA::InitPoissonDist(double Mean, double* Prob, double* CumProb, int N) {

  for (int i=0; i<N; i++) {
    Prob[i] = TMath::Power(Mean,i)*TMath::Exp(-Mean)/fFactorial[i];
  }
					// integral
  CumProb[0] = Prob[0];
  for (int i=1; i<MaxNx;i++) {
    CumProb[i] = CumProb[i-1]+Prob[i];
  }

}

//-----------------------------------------------------------------------------
void TFeldmanCousinsA::Init(double Bgr, double Sig) {
  //  fSigMean = SigMean;

  fMeanBgr = Bgr;
  fMeanSig = Sig;
  
  InitPoissonDist(fMeanBgr         , fBgProb, fCumBgProb, MaxNx);
  InitPoissonDist(fMeanBgr+fMeanSig, fBsProb, fCumBsProb, MaxNx);
//-----------------------------------------------------------------------------
// [re]-initialize 1D histograms with the probabilities and integral probabilities
//-----------------------------------------------------------------------------
  for (int i=0; i<MaxNx; i++) {
    fHist.fBgProb->SetBinContent(i+1,fBgProb[i]);

    fHist.fBsProb->SetBinContent(i+1,fBsProb[i]);
    fHist.fProb  ->SetBinContent(i+1,fBsProb[i]);
    
    fHist.fCumBgProb->SetBinContent(i+1,fCumBgProb[i]);
    fHist.fCumBsProb->SetBinContent(i+1,fCumBsProb[i]);
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

  for (int i1=0; i1<MaxNx-1; i1++) {
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

    if (fDebugLevel > 0) {
      PrintData("Rank   " ,'i',fRank     ,18);
    }
  }
//-----------------------------------------------------------------------------
// build confidence interval corresponding to the probability fCL - 0.9 , 0.95, etc
//-----------------------------------------------------------------------------
  fIMin       = MaxNx;
  fIMax       = -1;

  int covered = 0;
  fProb       = 0;
  
  for (int ix=0; ix<MaxNx; ix++) {
    int ind = fRank[ix];
    fProb  += fBsProb[ind];
    if (ind < fIMin) fIMin = ind;
    if (ind > fIMax) fIMax = ind;

    if (fDebugLevel > 0) {
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
// vary signal from SMin to SMax in NSteps, construct FC belt, fill belt histogram
// fBelt is the FC belt histogram
//-----------------------------------------------------------------------------
int TFeldmanCousinsA::ConstructBelt(double Bgr, double SMin, double SMax, int NSteps) {

  double step = (SMax-SMin)/NSteps;

  if (fHist.fBelt) delete fHist.fBelt;

  fHist.fBelt = new TH2D(Form("h_belt_%s",GetName()),"FC belt",TFeldmanCousinsA::MaxNx,0,TFeldmanCousinsA::MaxNx,
 	               NSteps,SMin-step/2,SMax-step/2);
  
  for (int iy=0; iy<NSteps; iy++) {
    double s = SMin+(iy+0.5)*step;
    ConstructInterval(Bgr,s);
    for (int ix=0; ix<MaxNx; ix++) {
      if ((ix < fIMin) || (ix > fIMax)) fBelt[iy][ix] = 0;
      else                              fBelt[iy][ix] = 1;

      fHist.fBelt->SetBinContent(ix+1,iy+1,fBelt[iy][ix]);
    }
  }
//-----------------------------------------------------------------------------
// for convenience, for each N, number of measured events, define the belt boundaries - fSBelt
//-----------------------------------------------------------------------------
  for (int ix=0; ix<MaxNx; ix++) {
    int iymin     = 0;
    int iymax     = 0;
    int inside    = 0;
    for (int iy=0; iy<NSteps; iy++) {
      if (fBelt[iy][ix] > 0) {
	if (inside == 0) {
	  iymin   = iy;
	  inside  = 1;
	}
      }
      else {
	if (inside == 1) {
	  iymax   = iy;
	  inside  = 0;
	  break;
	}
      }
    }

    fSBelt[0][ix] = (iymin+0.5)*step;
    fSBelt[1][ix] = (iymax+0.5)*step;

    printf("ix, smin, smax : %3i %12.5f %12.5f\n",ix,fSBelt[0][ix],fSBelt[1][ix]);
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
//-----------------------------------------------------------------------------
// generate background-only pseudo-experiments and plot the distribution
// for the excluded signal strength
//-----------------------------------------------------------------------------
  TH1D* h1 = new TH1D(Form("h_excluded_%s",GetName()),"excluded",2500,0,50);
  TRandom3 trn;
  for (int i=0; i<1000000; i++) {
    int rn = trn.Poisson(Bgr);

    if (fDebugLevel > 0) {
      printf("Bgr, rn, fSBelt[0][rn], fSBelt[1][rn] = %12.5e %5i %12.5e %12.5e\n",
	Bgr,rn,fSBelt[0][rn], fSBelt[1][rn]);
    }

    if (rn < MaxNx) {
      h1->Fill(fSBelt[1][rn]);
    }
    else {
      printf("trouble, Bill Robertson: rn = %3i\n",rn);
    }
  }

  TCanvas* c2 = new TCanvas("c2","plot h1",1000,800);
  c2->cd();
  h1->Draw();

  TCanvas* c_belt = new TCanvas("c_belt","Belt Histogram",1000,800);
  c_belt->cd();
  fHist.fBelt->Draw("box");

  double fc_upper_limit = h1->GetMean();
  printf(" mean excluded value : %12.5f\n",fc_upper_limit);
  return fc_upper_limit;
}

//-----------------------------------------------------------------------------
// in general, need to scan a range of signals, call this function multiple times
//-----------------------------------------------------------------------------
void TFeldmanCousinsA::DiscoveryProb(double Bgr, double SMin, double SMax, int NSteps) {
//-----------------------------------------------------------------------------
// construct FC CL confidence interval (covering fCL) assuming Signal=0
//-----------------------------------------------------------------------------
  ConstructInterval(Bgr,0);

  double step = (SMax-SMin)/NSteps;

  double x[10000], y[10000];

  int nx = NSteps+1;

  for (int ix=0; ix<nx; ix++) {
    double sig = SMin+ix*step;
    double tot = Bgr+sig;
//-----------------------------------------------------------------------------
// ndisc: number of "discovery experiments", pseudoexperiments in which NULL
// hypothesis is excluded at (1-fCL) level
//-----------------------------------------------------------------------------
    long int ndisc = 0;			
    for (int i=0; i<fNExp; i++) {
      int rn = fRn.Poisson(tot);
//-----------------------------------------------------------------------------
// definition of the discovery:
// rn > fIMax, i.e  the  probability to observe'rn' is less than 1-fCL
//-----------------------------------------------------------------------------
      if (rn > fIMax) ndisc ++;
    }
    
    double prob = double(ndisc)/double(fNExp);

    x[ix] = sig;
    y[ix] = prob;
  }

  TGraph* gr = new TGraph(nx,x,y);

  gr->SetName(Form("gr_bgr_%06i_%s",int(Bgr*1000),GetName()));
  gr->SetTitle(Form("FC discovery prob for bgr=%5.3f events",Bgr));

  gr->SetMarkerStyle(20);
  gr->Draw("alp");
}
