//-----------------------------------------------------------------------------
// 'cerc' : CE, rad-corrected
//-----------------------------------------------------------------------------
#include "TH1.h"
#include "TMath.h"

#include "math.h"
#include "Stntuple/val/stntuple_val_functions.hh"

class cerc  {
public:
  static double alpha;
  static double mmu;
  // recoil: mu^2/(2*M_Al)  M_al = 25133
  static double me;
  static double pi;

  TH1F*  h1;
  TH1F*  h2;

  cerc() {}
  ~cerc() {}

  static double f_rad_corr(double* X, double* P);
  
  TF1*   fun_rad_corr;

  void   rad_corr();
  void   print();
};

double cerc::alpha = 1./137.036;
double cerc::mmu   = 104.971;    // use Ee(max) instead of mmu (105.658);
double cerc::me    =   0.511;
double cerc::pi    = TMath::Pi();

//-----------------------------------------------------------------------------
// formula comes from mu2e-7615
//-----------------------------------------------------------------------------
double cerc::f_rad_corr(double* X, double* P) {
  double dg_de(0);
  double ee = X[0];
  if      ((ee < me ) || (ee >= mmu)) dg_de = 0;
  else {
    dg_de = alpha/(2*pi*mmu)*(log(4*(ee*ee)/(me*me))-2)*(ee*ee+mmu*mmu)/mmu/(mmu-ee);
    if (dg_de < 0) dg_de = 0;
  }
  
  return dg_de;
}

//-----------------------------------------------------------------------------
void cerc::rad_corr() {
  double x,y;
  double xmin(0),xmax(110.);

  double p[10]; // used only for interfacing

  int    nbx = 11000; // 11000;
  double bin = (xmax-xmin)/nbx;
  
  h1 = new TH1F("th","th",nbx,xmin,xmax);


  for (int i=1; i<=nbx; i++) {
    x = (i-0.5)*bin;
    y = cerc::f_rad_corr(&x,p);
    h1->SetBinContent(i,y);
  }

  double a0, a1, a2, a3, a01;
  double b1, b2, b3;

  printf("bin, nbx = %10.5f, %6i\n",bin,nbx);

  // find the last non-empty bin

  int ibmax = (mmu-xmin)/bin;
  

  a0  = h1->Integral()*bin;
  a01 = h1->Integral(0,ibmax-1)*bin;
//-----------------------------------------------------------------------------
// normalize the integral to one
//-----------------------------------------------------------------------------
  h1->SetBinContent(ibmax,(1-a01)/bin);
  
  a1 = 1-h1->Integral(0,(mmu-0.1)/bin)*bin;
  a2 = 1-h1->Integral(0,(mmu-1.5)/bin)*bin;
  a3 = 1-h1->Integral(0,(mmu-2.0)/bin)*bin;

  b1 = h1->Integral((mmu-0.1)/bin+1,nbx)*2*bin;
  b2 = h1->Integral((mmu-1.5)/bin+1,nbx)*2*bin;
  b3 = h1->Integral((mmu-2.0)/bin+1,nbx)*2*bin;

  printf("a0, a01 = %10.6f %10.6f\n",a0,a01);

  printf( " a1, a2, a3: %10.6lf %10.6lf %10.6lf\n",a1,a2,a3);
  printf( " b1, b2, b3: %10.6lf %10.6lf %10.6lf\n",b1,b2,b3);
//-----------------------------------------------------------------------------
// distribution from mu2e-7615, rebin to have bin = 20 keV
//-----------------------------------------------------------------------------
  h1->SetLineColor(2);
  h1->SetLineWidth(2);
  h1->Rebin(2);
  h1->Draw("");
//-----------------------------------------------------------------------------
// read PHOTOS histogram
//-----------------------------------------------------------------------------
  const char* photos_ce_hist = "/home/murat/hist/photos/muon_decay/photos_2.hist";

  h2 = (TH1F*) gh1(photos_ce_hist,"PhotosAna","genp_0/e_0")->Clone("photos");
  h2->Rebin(2);
  h2->SetMarkerStyle(20);
  h2->SetMarkerSize(0.5);
  h2->DrawNormalized("sames",h1->Integral());
}

//-----------------------------------------------------------------------------
void cerc::print() {
  double x, x1, y, y1;
  double xmin(0),xmax(110.);

  double p[10]; // used only for interfacing

  int    nbx = 105;
  double bin = 0.1;

  double integral = 0;
  for (double x=0; x <=105; x+=bin) {
    if (x < 104.8999) {
      y         = cerc::f_rad_corr(&x,p);
      x1        = x+bin;
      y1        = cerc::f_rad_corr(&x1,p);
    }
    else {
      y         = (1-integral)/bin;
      y1        = y;
    }
    integral += (y+y1)/2*bin;
    
    printf("%10.5f  %12.5e %12.5e \n",x,y,integral);
  }

}
