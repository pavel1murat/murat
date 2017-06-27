///////////////////////////////////////////////////////////////////////////////
//
// beam bins: 0: 3 GeV, 1: 5 GeV, 2: 8 GeV, 3: 12 GeV
//
// angular histograms : 0.15 - 2.15 radians, 10 bins, 0.2 rad bin
// momentum histograms - variable binning:
//
// 10 momentum slices:
//
// 0.10-0.15
// 0.15-0.20
// 0.20-0.25
// 0.25-0.30
// 0.30-0.35
// 0.35-0.40
// 0.45-0.50
// 0.50-0.60
// 0.60-0.70
// 0.70-0.80
//
// 9 angular slices:
//
// 0.35-0.55
// 0.55-0.75
// 0.75-0.95
// 0.95-1.15
// 1.15-1.35
// 1.35-1.55
// 1.55-1.75
// 1.75-1.95
// 1.95-2.15
//-----------------------------------------------------------------------------
#include "TH1.h"
#include "TString.h"
#include "TNtuple.h"

class HarpDataset {
public:
  enum { kNBeamBins       =  4 };
  enum { kNMomentumSlices = 11 };
  enum { kNThetaSlices    =  9 };

  struct Hist_t {
    TH1F*    fXsVsTheta    [kNBeamBins][kNMomentumSlices]; // momentum slices
    TH1F*    fXsVsMomentum [kNBeamBins][kNThetaSlices]; // theta    slices
  };

public:
  TString  fFn;				   // filename
  TString  fTarget;			   // 
  TString  fParticle;			   // 
  TNtuple* fNtuple;			   // ntuple with the data

  Hist_t   fHist;
					   // data
  float    tmin;
  float    tmax;
  float    pmin;
  float    pmax;
  float    xs [4];			   // cross sections (3,5,8,12 GeV)
  float    exs[4];			   // cross section errors

  float    fTheta   [kNThetaSlices+1];
  float    fMomentum[kNMomentumSlices+1];

  HarpDataset();
  HarpDataset(const char* Fn, const char* Target, const char* Particle);
  
  void  BookHistograms();
  void  FillMomentumHistograms();
  void  FillThetaHistograms   ();

  int   GetMomentumSlice(float PMin,  float PMax );
  int   GetThetaSlice   (float ThMin, float ThMax);
  int   GetBinNumber    (TH1* Hist, float X);

  void  InitLimits();
  
  void  PlotMomentumHist(int BeamBin, int ThetaSlice   );
  void  PlotThetaHist   (int BeamBin, int MomentumSlice);
};

//-----------------------------------------------------------------------------
HarpDataset::HarpDataset() {
}

//-----------------------------------------------------------------------------  
int HarpDataset::GetMomentumSlice(float PMin, float PMax) {
  int ind(-1);

  for (int i=0; i<kNMomentumSlices; i++) {
    if (fabs(PMin-fMomentum[i]) < 0.01) {
      ind = i;
      break;
    }
  }

  return ind;
}
  
//-----------------------------------------------------------------------------  
int HarpDataset::GetThetaSlice(float ThMin, float ThMax) {
  int ind(-1);
  
  for (int i=0; i<kNThetaSlices; i++) {
    if (fabs(ThMin-fTheta[i]) < 0.01) {
      ind = i;
      break;
    }
  }

  return ind;
}

//-----------------------------------------------------------------------------
int HarpDataset::GetBinNumber(TH1* Hist, float X) {
  int bin(-1);
  float xmin, xmax;
  int nbins = Hist->GetNbinsX();
  for (int i=1; i<=nbins; i++) {
    xmin = Hist->GetBinLowEdge(i);
    xmax = xmin+Hist->GetBinWidth(i);
    if ((X>=xmin) && (X<xmax)) {
      bin = i;
      break;
    }
  }
  return bin;
}

//-----------------------------------------------------------------------------
void HarpDataset::FillMomentumHistograms() {
  // 1. get theta slice - index of the dsigma/dp histograms

  int slice = GetThetaSlice(tmin,tmax);
  
  if (slice < 0) {
    printf(">>> ERROR: theta_slice = %i\n",slice);
  }
  else {
    float p = (pmin+pmax)/2;
    for (int ibeam=0; ibeam<4; ibeam++) {
      int bin = GetBinNumber(fHist.fXsVsMomentum[ibeam][0],p);
      if (bin > 0) {
	fHist.fXsVsMomentum[ibeam][slice]->SetBinContent(bin,xs [ibeam]);
	fHist.fXsVsMomentum[ibeam][slice]->SetBinError  (bin,exs[ibeam]);
      }
    }
  }
}

//-----------------------------------------------------------------------------
// xsection vs theta
//-----------------------------------------------------------------------------
void HarpDataset::FillThetaHistograms() {
  
  int slice = GetMomentumSlice(pmin,pmax);

  if (slice < 0) {
    printf(">>> ERROR: mom_slice = %i\n",slice);
  }
  else {
    float theta = (tmin+tmax)/2;
    for (int ibeam=0; ibeam<4; ibeam++) {
      int bin = GetBinNumber(fHist.fXsVsTheta[ibeam][0],theta);
      if (bin > 0) {
	fHist.fXsVsTheta[ibeam][slice]->SetBinContent(bin,xs [ibeam]);
	fHist.fXsVsTheta[ibeam][slice]->SetBinError  (bin,exs[ibeam]);
      }
    }
  }
}

//-----------------------------------------------------------------------------
void HarpDataset::InitLimits() {
  float theta_limit[kNThetaSlices+1] = {
    0.35, 0.55, 0.75, 0.95, 1.15, 1.35, 1.55, 1.75, 1.95, 2.15
  };

  float momentum_limit[kNMomentumSlices+1] = {
    0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.60, 0.70, 0.80
  };


  for (int i=0; i<kNThetaSlices+1; i++) {
    fTheta[i] = theta_limit[i];
  }
  
  for (int i=0; i<kNMomentumSlices+1; i++) {
    fMomentum[i] = momentum_limit[i];
  }
  
}

//-----------------------------------------------------------------------------
void HarpDataset::BookHistograms() {
  char name[200], title[200];
  
  float mom_lower[] = {
    0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.60, 0.70, 0.80
  };
  
  TH1::AddDirectory(0);

  for (int ibeam=0; ibeam<4; ibeam++) {
    for (int i=0; i<kNThetaSlices; i++) {
      sprintf(name ,"d2sigma_dp_dth_vs_mom_%i_%02i",ibeam,i);
      sprintf(title,"d2#sigma/dpd#theta vs P, %5.2f < #theta < %5.2f",fTheta[i],fTheta[i+1]);
      fHist.fXsVsMomentum[ibeam][i] = new TH1F(name,title,13,mom_lower);
      fHist.fXsVsMomentum[ibeam][i]->SetMarkerStyle(20);
      fHist.fXsVsMomentum[ibeam][i]->GetXaxis()->SetTitle("P, GeV/c");
    }
  }

  for (int ibeam=0; ibeam<4; ibeam++) {
    for (int i=0; i<kNMomentumSlices; i++) {
      sprintf(name,"d2sigma_dp_dth_vs_th_%i_%02i",ibeam,i);
      sprintf(title,"d^{2}#sigma/dpd#theta vs #theta, %4.2f < P < %4.2f",fMomentum[i],fMomentum[i+1]);
      fHist.fXsVsTheta[ibeam][i] = new TH1F(name,title,16,-0.05,3.15);
      fHist.fXsVsTheta[ibeam][i]->SetMarkerStyle(20);
      fHist.fXsVsTheta[ibeam][i]->GetXaxis()->SetTitle("#theta, rad");
    }
  }
}

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
HarpDataset::HarpDataset(const char* Fn, const char* Target, const char* Particle) {
  fFn       = Fn;
  fTarget   = Target;
  fParticle = Particle;
  fNtuple   = new TNtuple("harp","harp","tmin:tmax:pmin:pmax:x3:ex3:x5:ex5:x8:ex8:x12:ex12");
  fNtuple->ReadFile(Fn);

  InitLimits();
  BookHistograms();
  
  float* var;

  int n = fNtuple->GetEntries();
  for (int i=0; i<n; i++) {
    fNtuple->GetEntry(i);
    var    = fNtuple->GetArgs();
    //    printf("%10.3f %10.3f %10.3f\n",var[0],var[1],var[2]);
    tmin   = var[0];
    tmax   = var[1];
    pmin   = var[2];
    pmax   = var[3];
					// 4 cross sections - at 3,5,8, and 12 GeV/c
    xs [0] = var[ 4];
    exs[0] = var[ 5];
    xs [1] = var[ 6];
    exs[1] = var[ 7];
    xs [2] = var[ 8];
    exs[2] = var[ 9];
    xs [3] = var[10];
    exs[3] = var[11];

    FillMomentumHistograms();
    FillThetaHistograms   ();
  }
}

HarpDataset *PimTa, *PipTa, *PimPb, *PipPb;
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
void init() {
  const char* HarpDataDir = "./murat/data";
  PimTa = new HarpDataset(Form("%s/harp_piminus_Ta_data.txt",HarpDataDir),"Ta","pi-");
  PipTa = new HarpDataset(Form("%s/harp_piplus_Ta_data.txt" ,HarpDataDir),"Ta","pi+");
  PimPb = new HarpDataset(Form("%s/harp_piminus_Pb_data.txt",HarpDataDir),"Pb","pi-");
  PipPb = new HarpDataset(Form("%s/harp_piplus_Pb_data.txt" ,HarpDataDir),"Pb","pi+");
}
