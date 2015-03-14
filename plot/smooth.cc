///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "cmath"
#include "plot/smooth.hh"

//-----------------------------------------------------------------------------
double smooth::Eval(double* X) {

  int    i0, i1;
  double x, x0, dx0, dx1, f0, f1, f, bin, w0, w1;

  if (fHist == 0) {
    f = 0;
  }
  else {
    x  = X[0];

    bin = fHist->GetBinWidth(1);
    i0  = int ((x-fHist->GetXaxis()->GetXmin())/bin)+1 ;

    x0  = fHist->GetBinCenter(i0);

    dx0 = x-x0;
    if (dx0 >= 0) {
      i1  = i0+1;
      dx1 = dx0-bin;
      w0  = fabs(dx1)/bin;
      w1  = fabs(dx0)/bin;
    }
    else         {
      i1  = i0-1;
      dx1 = bin+dx0;
      if (i0 >= 2) {
	w0  = fabs(dx1)/bin;
	w1  = fabs(dx0)/bin;
      }
      else {
	w0  = 1;
	w1  = 0;
      }
    }
    
    f0  = fP0[i0] + dx0*fP1[i0] + dx0*dx0*fP2[i0];
    f1  = fP0[i1] + dx1*fP1[i1] + dx1*dx1*fP2[i1];
    
    //    f   = (f0*fabs(dx1)+f1*fabs(dx0))/bin;
    f   = f0*w0+f1*w1;

//   printf("%10.3f %3i %10.3f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n",
// 	 x,i0,x0,fP0[i0],fP1[i0],fP2[i0],dx0,f0,dx1,f1);
  }
  return f;
}


//-----------------------------------------------------------------------------
// before calling f_smooth (or the corresponding function make sure that 
// smooth::gVar points to the data structure
//-----------------------------------------------------------------------------
double smooth::func(double* X, double* P) {

  double f(0);
					// a hack...
  smooth* o = (smooth*) ((long int)P[0]);

  if (o) f = o->Eval(X);
  else   f = 0;

  return f;
}

//-----------------------------------------------------------------------------
// before calling set smooth::gVar to the address of an actual structure
//-----------------------------------------------------------------------------
smooth::smooth(const TH1* Hist, double XMin, double XMax) {

  if (Hist == 0) {
    fP0   = 0;
    fP1   = 0;
    fP2   = 0;
    fHist = 0;
    fFunc = 0;
  }
  else {

    int nx = Hist->GetNbinsX();

    double x1,x2,x3,y1,y2,y3;

    fP0   = new double[nx];
    fP1   = new double[nx];
    fP2   = new double[nx];
    fHist = (TH1*) Hist->Clone("Clone");

    if (XMin > XMax) {
      XMin = fHist->GetXaxis()->GetXmin();
      XMax = fHist->GetXaxis()->GetXmax();
    }

//   for (int i=1; i<=nx; i++) {
//     printf("%3i %10.5f %10.5f \n",i,Hist->GetBinCenter(i),Hist->GetBinContent(i));
//   }

    double bin = Hist->GetBinWidth(1);
    for (int i=2; i<=nx-1; i++) {
      x2 = Hist->GetBinCenter(i);
      x1 = x2-bin;
      x3 = x2+bin;
      
      y1 = Hist->GetBinContent(i-1);
      y2 = Hist->GetBinContent(i);
      y3 = Hist->GetBinContent(i+1);
      
      fP0[i] = y2;
      fP1[i] = (y3-y1)/2./bin;
      fP2[i] = (y1+y3-2*y2)/2/bin/bin;
      
    //    printf("%2i %10.3f %10.5f %10.5f %10.5f \n",i,x2,fP0[i],fP1[i],fP2[i]);
    }
//-----------------------------------------------------------------------------
// parametrization is finished, how to use it? 
//-----------------------------------------------------------------------------
//  if (gVar->func) delete gVar->func;
    fFunc = new TF1("hist_smooth_func",smooth::func,XMin,XMax,1);

    if (fHist == 0) {
      fFunc->SetParameter(0,0);
    }
    else {
      double x = (double) ((long int) this);
      fFunc->SetParameter(0,x);
    }
  }
}

//-----------------------------------------------------------------------------
smooth::~smooth() {

  if (fHist) {
    delete fP0;
    delete fP1;
    delete fP2;

    delete fFunc;
    delete fHist;
  }
}
