//


#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraph.h"

#include "plot/smooth.hh"
#include "plot/TMuonDecay.hh"

namespace {
  //  double emu(105.194), mAl(25133.), mmu(105.6583715), ;
  double me(0.511) ;
  //  double z(13.), alpha(1/137.036);

  //  TCanvas* c;

  //  TH1D* Hist(0);

  //  TF1   *f1(0), *f2(0);
  
  //  char fn_eps[200],fn_pdf[200],cmd[500];
};


double  TMuonDecay::fMAl   = 25133.;
double  TMuonDecay::fMmu   =   105.6583715;
double  TMuonDecay::fEmu   =   105.194;
double  TMuonDecay::fMe    =     0.511;
double  TMuonDecay::fZAl   =    13.;
double  TMuonDecay::fAlpha =  1/137.036;
double  TMuonDecay::fAPi   = TMuonDecay::fAlpha/M_PI;

ClassImp(TMuonDecay)

//-----------------------------------------------------------------------------
TMuonDecay::TMuonDecay() : TNamed("","") {
}

//-----------------------------------------------------------------------------
double TMuonDecay::ff(double* X, double* P) {
  double x = X[0];
  
  double f = log(1-x)/x;
  return -f;
}



//-----------------------------------------------------------------------------
double TMuonDecay::Li2(double X) {

  static int         first_time(1);
  static TH1F*       h1;
  static smooth*     smf;
  static TF1*        f_ff;

  //  static TGraph      *gr;

  double    x, y;
  int const np(10000);

  if (first_time) {
    first_time = 0;
    
    f_ff = new TF1("f_ff",ff,0,1,0);

    f_ff->SetNpx(np);

    double step = 1./np;

    h1 = new TH1F("h1","h1",np,0,1);
    
    for (int i=0; i<np; i++) {
      x = (i+0.5)*step;
      y = f_ff->Integral(0,x);
      h1->SetBinContent(i,y);
    }

    smf = new smooth(h1,0,1);
  
    // gr = new TGraph(np,x,y);
    // gr->Draw();
  }

  double f = smf->GetFunc()->Eval(X);

  return f;
}


//-----------------------------------------------------------------------------
double TMuonDecay::f_li2(double* X, double* P) {
  double x = X[0];
  return Li2(x);
}

//-----------------------------------------------------------------------------
double TMuonDecay::f_val_li2(double* X, double* P) {
  double x = X[0];
  return Li2(x)+Li2(1-x)+log(x)*log(1-x);
}

//-----------------------------------------------------------------------------
void TMuonDecay::plot_li2() {
  TF1* f = new TF1("f_li2",f_li2,0,1,0);
  f->Draw();
}


//-----------------------------------------------------------------------------
void TMuonDecay::plot_val_li2() {
  TF1* f = new TF1("f_val_li2",f_val_li2,0,1,0);
  f->Draw();
}


//-----------------------------------------------------------------------------
// mu -> e nu nu decay spectrum in the leading order
// skip (G^2/12pi^2)
// from D.Ross, Nuovo Cimento v10A,3,475(1972)
//-----------------------------------------------------------------------------
double TMuonDecay::muon_spectrum_lo_ross(double E) {

  double f(0), E2, me2, mmu, mmu2, emax;

  E2   = E*E;
  me2  = me*me;
  mmu  = TMuonDecay::fMmu;
  mmu2 = mmu*mmu;
  emax = TMuonDecay::fMmu/2 + me2/(2*mmu);

  if ((E >= me) && (E < emax)) f = sqrt(E2-me2)*(E*(3*(mmu2+me2)-4*E*mmu)+2*mmu*me2);

  //  printf("E = %12.3f f = %12.3f  P[0] = %12.5f\n",E,f,P[0]);
  
  return f;
}

//-----------------------------------------------------------------------------
double TMuonDecay::f_muon_spectrum_lo_ross(double* X, double* P) {

  double e, f;

  e    = X[0];
  f    = P[0]*muon_spectrum_lo_ross(e);
  return f;
}


//-----------------------------------------------------------------------------
// mu -> e nu nu decay spectrum, Born approximation
// from http://arxiv.org/abs/hep-ph/0202102 (Arbuzov, Czarnecki, Gaponenko)
//-----------------------------------------------------------------------------
double TMuonDecay::muon_spectrum_lo_arbuzov(double E) {

  double f, emax, me2, mmu, r, r2, v, z, z0;

  mmu  = TMuonDecay::fMmu;
  me2  = me*me;
  //  mmu2 = mmu*mmu;

  r    = TMuonDecay::fMe/TMuonDecay::fMmu;
  r2   = r*r;
  emax = mmu/2*(1+r2);
  z    = E/emax;
  z0   = 2*r/(1.+r2);
  v    = sqrt(1.-me2/(E*E));
  
  if ((E >= me) && (E < emax)) {
    f = 6*pow(1+r2,4)*v*z*(z*(1-z)+(2./9)*(3./4)*(4*z*z-3*z-z0*z0));
  }
  else                         f = 0;

  //  printf("E = %12.3f f = %12.3f  P[0] = %12.5f\n",E,f,P[0]);
  
  return f;
}

//-----------------------------------------------------------------------------
// mu -> e nu nu decay spectrum, Born approximation
// from http://arxiv.org/abs/hep-ph/0202102 (Arbuzov, Czarnecki, Gaponenko)
//-----------------------------------------------------------------------------
double TMuonDecay::f_muon_spectrum_lo_arbuzov(double* X, double* P) {

  double f;

  f = P[0]*muon_spectrum_lo_arbuzov(X[0]);

  //  printf("E = %12.3f f = %12.3f  P[0] = %12.5f\n",E,f,P[0]);
  
  return f;
}

//-----------------------------------------------------------------------------
// mu -> e nu nu decay spectrum, Born approximation
// from http://arxiv.org/abs/hep-ph/0202102 (Arbuzov, Czarnecki, Gaponenko)
//-----------------------------------------------------------------------------
double TMuonDecay::R1Arbuzov(double Z) {
  double f, logz, log1mz;

  logz   = log(Z);
  log1mz = log(1-Z);

  f = -2*Li2(1-Z)+logz*log1mz-2*logz*logz-log1mz/Z-5./4;

  return f;
}


//-----------------------------------------------------------------------------
double TMuonDecay::muon_spectrum_nlo_arbuzov(double E) {

  double   a2pi, ax, bx, hx, f, f0, f1LL, emax, L, me2, mmu, mmu2, R1, r, r2, z, z2, z3;
  double   logz, log1mz;

  //  E2   = E*E;
  mmu    = TMuonDecay::fMmu;
  me2    = me*me;
  mmu2   = mmu*mmu;
  L      = log(mmu2/me2);

  r      = me/mmu;
  r2     = r*r;
  emax   = mmu/2*(1+r2);
  z      = E/emax;
  z2     = z*z;
  z3     = z2*z;
  logz   = log(z);
  log1mz = log(1-z);
  a2pi   = TMuonDecay::fAlpha/(2*M_PI);
  
  //  z0   = 2*r/(1.+r2);
  //  v    = sqrt(1.-me2/E2);
  
  if ((E >= me) && (E < emax)) {
    f0   = muon_spectrum_lo_arbuzov(E);

    f1LL = 5./6+2*z-4*z2+8./3*z3+2*z2*(3-2*z)*log((1-z)/z);

    R1   = -2*Li2(1-z)+logz*log1mz-2*logz*logz-log1mz/z-5./4;

    //    f1   = (L-1)*f1LL + 2*z2*(3-2*z)*R1 + (1-z)/6.*((10+34*z-32*z2)*log(z)+(5-27*z+34*z2));

    ax   = -f1LL + 2*z2*(3-2*z)*R1 + (1-z)/6.*((10+34*z-32*z2)*log(z)+(5-27*z+34*z2));

    bx   = f1LL;
      
    hx   = ax + L*bx;

    f    = f0*(1 + a2pi*hx/f0);
  }
  else                         f = 0;

  //  printf("E = %12.3f f = %12.3f  P[0] = %12.5f\n",E,f,P[0]);
  
  return f;
}

//-----------------------------------------------------------------------------
double TMuonDecay::f_muon_spectrum_nlo_arbuzov(double* X, double* P) {
  double f;

  f = P[0]*muon_spectrum_nlo_arbuzov(X[0]);

  return f;
}


//-----------------------------------------------------------------------------
// mu -> e nu nu decay spectrum in the NLO
// http://link.springer.com/article/10.1134%2FS1547477109050033
//-----------------------------------------------------------------------------
double TMuonDecay::f_muon_spectrum_nlo_kuraev(double* X, double* P) {

  double           f, E, E2, L, logx, log1mx, me2, mmu, mmu2, emax, x;
  double           p_lo_1, p_lo_2;
  double           a1, a2, a3, ax, bx, hx;
  double           a2pi; 
  double           pi = M_PI;
  
  E      = X[0];
  E2     = E*E;
  me2    = me*me;
  mmu    = TMuonDecay::fMmu;
  mmu2   = mmu*mmu;
  a2pi   = fAlpha/(2*pi);
					// the paper has a typo - they use ln(mmu^2/me^) instead of ln(mmu/me)
  L      = log(mmu/me);

  x      = 2*E/mmu;
  logx   = log(x);
  log1mx = log(1-x);
  
  emax   = mmu/2*(1 + me2/mmu2);

  if ((E >= me) && (E < emax)) {
    
    p_lo_1 = sqrt(E2-me2);
    p_lo_2 = (E*(3*(mmu2+me2)-4*E*mmu)+2*mmu*me2);

    a1     = 4*Li2(x)-2*pi*pi/3.-4+2*logx*(3*log1mx-2*logx+1)-2*(1+1/x)*log1mx;
    a2     = (1-x)*(5+17*x-16*x*x)/(3*x*x*(3-2*x))*log(x);
    a3     = (1-x)*(-22*x+34*x*x)/(3*x*x*(3-2*x));
    
    ax     = a1+a2+a3;

    bx     = 3.+4*(log1mx-logx)+(1-x)*(5+17*x-34*x*x)/(3*x*x*(3-2*x));

    hx     = ax+L*bx;
    
    f      = P[0]*p_lo_1*p_lo_2*(1+a2pi*hx);
    
  }
  else   f = 0;

  //  printf("E = %12.3f f = %12.3f  P[0] = %12.5f\n",E,f,P[0]);
  
  return f;
}


//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
void TMuonDecay::plot_muon_spectrum_lo(const char* Opt) {
  
  TF1* f_lo = new TF1("f_muon_spectrum_lo",f_muon_spectrum_lo_arbuzov,me,55.,1);
  f_lo->SetParameter(0,1.);

  f_lo->SetNpx(1000);

  double qint = f_lo->Integral(0.,55.);
  printf("qint = %12.5f\n",qint);
  
  double p0 = f_lo->GetParameter(0);
  printf("p0 = %12.5g\n",p0);
  
  f_lo->SetParameter(0,1./qint);

  f_lo->SetLineWidth(1);
  f_lo->SetLineColor(kBlue);
  f_lo->Draw(Opt);
}

//-----------------------------------------------------------------------------
void TMuonDecay::plot_muon_spectrum_nlo(const char* Name, const char* Opt) {

  TF1* f_nlo(NULL); 
  char name[100];

  sprintf(name,"f_muon_spectrum_nlo_name_%s",Name);

  if (strcmp(Name,"kuraev") == 0) {
    f_nlo = new TF1(name,f_muon_spectrum_nlo_kuraev,me,55.,1);
  }
  else if (strcmp(Name,"arbuzov") == 0) {
    f_nlo = new TF1(name,f_muon_spectrum_nlo_arbuzov,me,55.,1);
  }

  f_nlo->SetParameter(0,1.);

  f_nlo->SetNpx(1000);

  double qint = f_nlo->Integral(0.,55.);
  printf("qint = %12.5f\n",qint);
  
  double p0 = f_nlo->GetParameter(0);
  printf("p0 = %12.5g\n",p0);
  
  f_nlo->SetParameter(0,1./qint);

  f_nlo->SetLineWidth(1);
  f_nlo->Draw(Opt);
}
