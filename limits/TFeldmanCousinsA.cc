// my implementation of the Feldman-Cousins algorithm
// 
#include "murat/limits/TFeldmanCousinsA.hh"

#include "TCanvas.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TROOT.h"
#include "Math/QuantFuncMathCore.h"
#include "Math/DistFunc.h"
//-----------------------------------------------------------------------------
// CL > 0: CL has the meaning of a probability - 0.9, 0.95 .. .0.99 etc
//    < 0: CL is the "discovery probability", probability corresponding
//         to 5 gaussian sigma level, 
//-----------------------------------------------------------------------------
ClassImp(murat::TFeldmanCousinsA)

namespace murat {

  TFeldmanCousinsA::TFeldmanCousinsA(const char* Name, double CL, int DebugLevel):
  TNamed(Name,Name),
  fRn()
{
  SetCL(CL);

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

  fHist.fProb      = new TH1D(Form("h_prob_2D_%s"    ,GetName()),"h prob 2D" ,MaxNx,-0.5,float(MaxNx)-0.5);
  fHist.fBelt      = nullptr;

  fNExp            = 10000000;
					// make sure uninitialzied values don't make sense
  fBelt.fNy        = -1;
  fBelt.fSMin      = 1e6;
  fBelt.fSMax      = -1e6;
}

TFeldmanCousinsA::~TFeldmanCousinsA() {
  delete fHist.fBgProb;
  delete fHist.fBsProb;
  delete fHist.fCumBgProb;
  delete fHist.fCumBsProb;
  delete fHist.fProb;
  if (fHist.fBelt) delete fHist.fBelt;
}

void  TFeldmanCousinsA::SetCL(double CL) {
  double alpha = 1-TMath::Erf(5./sqrt(2));
  
  if (CL > 0) fCL = CL;
  else        fCL = 1-alpha/2; // always two-sided: 1-5.7330314e-07

  fLog1mCL        = log(1-fCL);
}

void TFeldmanCousinsA::InitPoissonDist(double Mean, double* Prob, double* CumProb, int N) {
//-----------------------------------------------------------------------------
// the length of 'Prob' should be at least N
// CumProb[N]: probability, for a given Mean, to have rn <= N
//-----------------------------------------------------------------------------

  Prob[0] = TMath::Exp(-Mean);
  for (int i=1; i<N; i++) {
    Prob[i] = Prob[i-1]*Mean/i;
  }
					// integral
  CumProb[0] = 1;
  for (int i=1; i<N;i++) {
    CumProb[i] = 0;
    for (int k=i; k<N; k++) {
      CumProb[i] += Prob[k];
    }
  }
}

//-----------------------------------------------------------------------------
  void TFeldmanCousinsA::Init(double MuB, double MuS, int FillHist) {
  //  fSigMean = SigMean;

  fMeanBgr = MuB;
  fMeanSig = MuS;
  
  InitPoissonDist(MuB    ,fBgProb,fCumBgProb,MaxNx);
  InitPoissonDist(MuB+MuS,fBsProb,fCumBsProb,MaxNx);
//-----------------------------------------------------------------------------
// [re]-initialize 1D histograms with the probabilities and integral probabilities
//-----------------------------------------------------------------------------
  if (FillHist) {
    for (int i=0; i<MaxNx; i++) {
      fHist.fBgProb->SetBinContent(i+1,fBgProb[i]);

      fHist.fBsProb->SetBinContent(i+1,fBsProb[i]);
      fHist.fProb  ->SetBinContent(i+1,fBsProb[i]);
    
      fHist.fCumBgProb->SetBinContent(i+1,fCumBgProb[i]);
      fHist.fCumBsProb->SetBinContent(i+1,fCumBsProb[i]);
    }
  }
}

int TFeldmanCousinsA::ConstructInterval(double MuB, double MuS) {
//-----------------------------------------------------------------------------
// output: [fIMin,fIMax] : a CL interval constructed using FC ordering for given MuB and MuS
//-----------------------------------------------------------------------------
  int rc(0);				// return code

  Init(MuB,MuS);
  
  for (int ix=0; ix<MaxNx; ix++) {
    double sbest = ix-MuB;

    if (sbest <= 0) sbest = 0;

    double sb = sbest+MuB;

    fBestProb[ix] = TMath::Power(sb,ix)*TMath::Exp(-sb)/fFactorial[ix];
    fBestSig [ix] = sbest;
    fLhRatio [ix] = fBsProb[ix]/fBestProb[ix];
  }
//-----------------------------------------------------------------------------
// sort ranks
//-----------------------------------------------------------------------------
  double rmax;
    
  for (int i=0; i<MaxNx; i++) fRank[i] = i;

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
// vary signal from SMin to SMax in NPoints, construct FC belt, fill belt histogram
// fBelt is the FC belt histogram
// avoid multiple useless re-initializations
//-----------------------------------------------------------------------------
int TFeldmanCousinsA::ConstructBelt(double Bgr, double SMin, double SMax, int NPoints) {

  if ((fBelt.fBgr == Bgr) and (fBelt.fNy == NPoints) and (fBelt.fSMin == SMin) and (fBelt.fSMax == SMax)) return 0;

  fBelt.fBgr  = Bgr;
  fBelt.fSMin = SMin;
  fBelt.fSMax = SMax;
  fBelt.fNy   = NPoints;
  
  fBelt.fDy   = (NPoints > 1) ? (SMax-SMin)/(NPoints-1) : 1;

  // if (fHist.fBelt) delete fHist.fBelt;

  // fHist.fBelt = new TH2D(Form("h_belt_%s",GetName()),"FC belt",TFeldmanCousinsA::MaxNx,0,TFeldmanCousinsA::MaxNx,
  // 			 NPoints,SMin-step/2,SMax-step/2);
  
  for (int iy=0; iy<NPoints; iy++) {
    double mus = SMin+(iy+0.5)*fBelt.fDy;
    ConstructInterval(Bgr,mus);
    for (int ix=0; ix<MaxNx; ix++) {
      if ((ix < fIMin) || (ix > fIMax)) fBelt.fCont[iy][ix] = 0;
      else                              fBelt.fCont[iy][ix] = 1;
    }
  }
//-----------------------------------------------------------------------------
// for convenience, for each N, number of measured events, define the belt boundaries - fSBelt
//-----------------------------------------------------------------------------
  for (int ix=0; ix<MaxNx; ix++) {
    int iymin     = 0;
    int iymax     = 0;
    int inside    = 0;
    for (int iy=0; iy<NPoints; iy++) {
      if (fBelt.fCont[iy][ix] > 0) {
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
					// probably want step = 10^-3
    fBelt.fSign[ix][0] = SMin+(iymin+0.5)*fBelt.fDy;
    fBelt.fSign[ix][1] = SMin+(iymax+0.5)*fBelt.fDy;
    fBelt.fIndx[ix][0] = iymin;
    fBelt.fIndx[ix][1] = iymax;

    if (fDebugLevel > 0) {
      printf("ix, smin, smax : %3i %12.5f %12.5f\n",ix,fBelt.fSign[ix][0],fBelt.fSign[ix][1]);
    }
  }

  return 0;
}

//-----------------------------------------------------------------------------
// in general, need to scan a range of signals, call this function multiple times
//-----------------------------------------------------------------------------
void TFeldmanCousinsA::DiscoveryProb(double MuB, double SMin, double SMax, int NPoints, double* MuS, double* Prob) {
//-----------------------------------------------------------------------------
// construct FC CL confidence interval (covering fCL) assuming Signal=0
//-----------------------------------------------------------------------------
  ConstructInterval(MuB,0);

  double step = (NPoints > 1) ? (SMax-SMin)/(NPoints-1) : 0;

  for (int ix=0; ix<NPoints; ix++) {
    MuS[ix]   = SMin+ix*step;
    double tot = MuB+MuS[ix];
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
    Prob[ix] = double(ndisc)/double(fNExp);
  }
}

//-----------------------------------------------------------------------------
// in general, need to scan a range of signals, call this function multiple times
//-----------------------------------------------------------------------------
  void TFeldmanCousinsA::DiscoveryProb(murat::analysis* Analysis, double SMin, double SMax, int NPoints, double* MuS, double* Prob) {
//-----------------------------------------------------------------------------
// construct FC CL confidence interval (covering fCL) assuming Signal=0
//-----------------------------------------------------------------------------
  double step = (NPoints > 1) ? (SMax-SMin)/(NPoints-1) : 0;

  for (int ix=0; ix<NPoints; ix++) {
    MuS[ix]   = SMin+ix*step;
//-----------------------------------------------------------------------------
// ndisc: number of "discovery experiments", pseudoexperiments in which NULL
// hypothesis is excluded at (1-fCL) level
//-----------------------------------------------------------------------------
    long int ndisc = 0;			
    for (int i=0; i<fNExp; i++) {
//-----------------------------------------------------------------------------
// when integrating over the nuisance parameters, redefine the interval each time
//-----------------------------------------------------------------------------
      double bgr = Analysis->FluctuateBackground();
      double tot = bgr+MuS[ix];
      int rn = fRn.Poisson(tot);
//-----------------------------------------------------------------------------
// definition of the discovery:
// rn > fIMax, i.e  the  probability to observe'rn' is less than 1-fCL
//-----------------------------------------------------------------------------
      ConstructInterval(bgr,0);

      if (rn > fIMax) ndisc ++;
    }
    Prob[ix] = double(ndisc)/double(fNExp);
  }
}


void TFeldmanCousinsA::DiscoveryProbMean(double MuB, double SMin, double SMax, int NPoints, double* MuS, double* Prob) {
//-----------------------------------------------------------------------------
// in general, need to scan a range of signals, call this function multiple times
// watch for 5 
// construct FC CL confidence interval (covering fCL) for MuS=0
//-----------------------------------------------------------------------------
  double mus;
  ConstructInterval(MuB,mus=0);

  double step = (NPoints > 1) ? (SMax-SMin)/(NPoints-1) : 0;

  for (int ix=0; ix<NPoints; ix++) {
    MuS[ix]    = SMin+ix*step;
    double tot = MuB+MuS[ix];
//-----------------------------------------------------------------------------
// ndisc: number of "discovery experiments", pseudoexperiments in which NULL
// hypothesis is excluded at (1-fCL) level
//-----------------------------------------------------------------------------
    double sum  = 0;
    double sumn = 0;
    for (int i=0; i<fNExp; i++) {
      int rn = fRn.Poisson(tot);
//-----------------------------------------------------------------------------
// define probability for the background to fluctuate above rn
//-----------------------------------------------------------------------------
      double p;
      if (rn > 0) p   = ROOT::Math::poisson_cdf_c(rn-1,MuB);
      else        p   = ROOT::Math::poisson_cdf_c(0   ,MuB);

      double      sig = ROOT::Math::gaussian_quantile_c(p,1);

      sum      += sig;
      sumn     += 1;
      if (fDebugLevel > 0) {
	printf("i, rn, p, sig, sum : %3i %3i %15.8e %12.5e %12.5e\n",i,rn, p,sig, sum);
      }
    }
    
    Prob[ix]  = sum/sumn; // ROOT::Math::gaussian_quantile(pp,1);

    if (fDebugLevel > 0) {
      printf("ix, sum, sumn, MuS[ix], Prob[ix] : %3i %12.5e %10.3e %12.5e %12.5e\n",
	     ix,sum,sumn,MuS[ix],Prob[ix]);
    }
  }
}

void TFeldmanCousinsA::PlotDiscoveryProbMean(double MuB, double MuS) {
  // double mus;
  // ConstructInterval(MuB,mus=0);

  TH1F* h_sigm = new TH1F("h_sigm","h1 sigm"     ,2000,    0,20);
  TH1F* h_prob = new TH1F("h_prob","h1 log(prob)",2000, - 99, 1);
//-----------------------------------------------------------------------------
// ndisc: number of "discovery experiments", pseudoexperiments in which NULL
// hypothesis is excluded at (1-fCL) level
//-----------------------------------------------------------------------------
  if (fDebugLevel > 0) {
    printf("  i  rn cum_prob[rn]        p          sig          sum       \n");
    printf("--------------------------------------------------------------\n");
  }
  for (int i=0; i<fNExp; i++) {
    int rn = fRn.Poisson(MuB+MuS);
    //    int rn = fRn.Poisson(MuS);
//-----------------------------------------------------------------------------
// define probability for the background to fluctuate above rn
//-----------------------------------------------------------------------------
    double p;
    if (rn > 0) p   = ROOT::Math::poisson_cdf_c(rn-1,MuB);
    else        p   = ROOT::Math::poisson_cdf_c(0   ,MuB);
    
    double      sig = ROOT::Math::gaussian_quantile_c(p,1);
      
    h_sigm->Fill(sig);
    h_prob->Fill(log(p));
  }

  const char* name = "c_fc_PlotDiscoveryProbMean";
  
  TCanvas* c = (TCanvas*) gROOT->GetListOfCanvases()->FindObject(name);
  if (c == nullptr) {
    c = new TCanvas("c_fc_PlotDiscoveryProbMean","c",1500,500);
    c->Divide(2,1);
  }
  c->cd(1);
  h_sigm->Draw();
  c->cd(2);
  h_prob->Draw();
}


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

void TFeldmanCousinsA::PrintProbs(int N) {
//-----------------------------------------------------------------------------
// print results
//-----------------------------------------------------------------------------
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


int TFeldmanCousinsA::SolveFor(double Val, const double* X, const double* Y, int NPoints, double* XVal) {
//-----------------------------------------------------------------------------
// fiven NPoints (X,Y) representing a smooth curve f(X), find XVal,
// such that f(XVal) = Val
//-----------------------------------------------------------------------------
  int i0 = -1;
  for (int i=0; i< NPoints; i++) {
    if (fDebugLevel > 0) {
      printf("i, X[i], Y[i] : %2i %10.3e %12.5e\n",i,X[i],Y[i]);
    }
    if (Y[i] < Val) continue;
    i0 = i-1;
    break;
  }

  if (i0 == -1) {
    printf("TFeldmanCousinsA::solve_for in trouble: wrong range, bail out\n");
    return -1;
  }
  
  int i1 = i0+1;
  int i2 = i0+2;
  
  TMatrixD m(3,3);

  m(0,0) = X[i0]*X[i0]; m(0,1) = X[i0]; m(0,2) = 1;
  m(1,0) = X[i1]*X[i1]; m(1,1) = X[i1]; m(1,2) = 1;
  m(2,0) = X[i2]*X[i2]; m(2,1) = X[i2]; m(2,2) = 1;

  TMatrixD minv(TMatrixD::kInverted,m);

  TVectorD vy(3), vp(3);
  vy(0) = Y[i0];
  vy(1) = Y[i1];
  vy(2) = Y[i2];

  vp = minv*vy;

  double a = vp(0);
  double b = vp(1);
  double c = vp(2);
//-----------------------------------------------------------------------------
// finally, the solution for XVal
//-----------------------------------------------------------------------------
  if (a < 0) *XVal = -b/(2*a) - sqrt(b*b/(4*a*a)-(c-Val)/a);
  else       *XVal = -b/(2*a) + sqrt(b*b/(4*a*a)-(c-Val)/a);

  if (fDebugLevel > 0) {
    printf ("a,b,c,x = %12.5e %12.5e %12.5e %12.5e \n",a,b,c,*XVal);
  }
  return 0;
}



void TFeldmanCousinsA::UpperLimit(double MuB, double SMin, double SMax, int NPoints, double* S, double* Prob) {
//-----------------------------------------------------------------------------
// construct FC CL confidence belt
// constructing the FC belt assumes dividing the signal interval into (NPoints-1) steps,
// use 1000 points, assume signal in the range [0,10]
//-----------------------------------------------------------------------------
  ConstructBelt(MuB,0,10,1001);
//-----------------------------------------------------------------------------
// ndisc: number of "discovery experiments", pseudoexperiments in which NULL
// hypothesis is excluded at (1-fCL) level
//-----------------------------------------------------------------------------
  int nsteps = (NPoints > 1) ? NPoints -1 : 1;
  double step = (SMax-SMin)/nsteps;
  
  for (int is=0; is<NPoints; is++) {
    S[is] = SMin+is*step;
    long int nexcl = 0;			
    for (int i=0; i<fNExp; i++) {
      int rn = fRn.Poisson(MuB);
//-----------------------------------------------------------------------------
// definition of the discovery:
// rn > fIMax, i.e  the  probability to observe'rn' is less than 1-fCL
//-----------------------------------------------------------------------------
      if (S[is] > fBelt.fSign[rn][1]) nexcl++;
      if (fDebugLevel > 0) {
	printf(" MuB, rn S[is] fBelt.fSign[rn][1] : %12.5e %3i %12.5e %12.5e\n",
	       MuB, rn, S[is],fBelt.fSign[rn][1]);
      }
    }
    Prob[is] = double(nexcl)/double(fNExp);
  }
}
}
