//

/* --

Double_t CrystalBall(Double_t *x,Double_t *par) {

  Double_t t = (x[0]-par[2])/par[3];
  if (par[0] < 0) t = -t;

  Double_t absAlpha = fabs((Double_t)par[0]);

  if (t >= -absAlpha) {
    return par[4]*exp(-0.5*t*t);
  }
  else {
    Double_t a =  TMath::Power(par[1]/absAlpha,par[1])*exp(-0.5*absAlpha*absAlpha);
    Double_t b= par[1]/absAlpha - absAlpha; 

    return par[4]*(a/TMath::Power(b - t, par[1]));
  }
}

-- */

namespace {
  TH1F* h; 
  TF1*  f; 
};


//-----------------------------------------------------------------------------
// A1*exp((x-x0)^2/2sig^2 for (x-x0)/sig > -alp
// A1*A2*(B-(x-x0)/sig)^{-n}
// A2 = 
// P[0] = N
// P[1] = x[0]
// P[2] = sigma
// P[3] = alpha  (> 0)
// P[4] = n
// P[5] = tail fraction
// P[6] = tail exponential slope
//-----------------------------------------------------------------------------
double f_crystal_ball_save_v0(double* X, double* P) {
  double f, alpha, abs_alpha, n, a, b;

  double dx0 = (X[0]-P[1]);
  double dx = (X[0]-P[1])/P[2];

  alpha     = P[3];
  abs_alpha = fabs(alpha);
  n         = P[4];

  if (dx0 > -abs_alpha) {
    if (dx < 0) {
      f = P[0]*TMath::Exp(-dx*dx/2.);
    }
    else if (dx0 < 0.5) {
      f = P[0]*(1-P[5])*TMath::Exp(-dx*dx/2.);
    }
    else {
      f = P[0]*((1-P[5])*TMath::Exp(-dx*dx/2.)+P[5]*TMath::Exp(-dx/P[6]));
    }
  }
  else {
    a = TMath::Power((n/abs_alpha),n)*TMath::Exp(-(alpha*alpha)/2);
    b = n/abs_alpha-abs_alpha;
    f = P[0]*a*TMath::Power((b-dx),-n);
  }
  
  return f;
}

//-----------------------------------------------------------------------------
// try two gaussians - an exponential tail extending only in one direction
// doesn't make much sense
//-----------------------------------------------------------------------------
double f_crystal_ball_save(double* X, double* P) {
  double f, alpha, abs_alpha, n, a, b;
  double dx = (X[0]-P[1])/P[2];

  alpha     = P[3];
  abs_alpha = fabs(alpha);
  n         = P[4];

  if (dx > -abs_alpha) {
    double dx1 = fabs((X[0]-P[1])/P[6]);
    f = P[0]*((1-P[5])*TMath::Exp(-dx*dx/2.)+P[5]*TMath::Exp(-dx1));
  }
  else {
    a = TMath::Power((n/abs_alpha),n)*TMath::Exp(-(alpha*alpha)/2);
    b = n/abs_alpha-abs_alpha;
    f = P[0]*a*TMath::Power((b-dx),-n);
  }
  
  return f;
}



//-----------------------------------------------------------------------------
// try two gaussians - an exponential tail extending only in one direction
// doesn't make much sense
//-----------------------------------------------------------------------------
double f_crystal_ball(double* X, double* P) {
  double f, alpha, abs_alpha, n, a, b;

  double dx  = (X[0]-P[1])/P[2];
  double dx1 = (X[0]-P[1])/P[6];

  alpha     = P[3];
  abs_alpha = fabs(alpha);
  n         = P[4];

  if (dx > -abs_alpha) {
    f = P[0]*((1-P[5])*TMath::Exp(-dx*dx/2.)+P[5]*TMath::Exp(-dx1*dx1/2.));
  }
  else {
    a = TMath::Power((n/abs_alpha),n)*TMath::Exp(-(alpha*alpha)/2);
    b = n/abs_alpha-abs_alpha;
    f = P[0]*a*TMath::Power((b-dx),-n);
  }
  
  return f;
}

//-----------------------------------------------------------------------------
void cb_init_parameters(TF1* F, 
			double P0, double P1, double P2, double P3, 
			double P4, double P5, double P6) 
{
  F->SetParameter(0,P0);
  F->SetParameter(1,P1);
  F->SetParameter(2,P2);
  F->SetParameter(3,P3);
  F->SetParameter(4,P4);
  F->SetParameter(5,P5);
  F->SetParameter(6,P6);
}

//-----------------------------------------------------------------------------
void cb_fit(TH1F* Hist, TF1* F, double XMin, double XMax) {
  Hist->GetXaxis()->SetRangeUser(XMin,XMax);
  Hist->Fit(F,"L","",XMin,XMax);
}

//-----------------------------------------------------------------------------
void create_fit_function(TF1*& F, double X0, double XMin, double XMax) {

  F = new TF1("f_crystal_ball",f_crystal_ball,XMin,XMax,7);

  F->SetParName(0,"anorm");
  F->SetParName(1,"mean");
  F->SetParName(2,"sigma");
  F->SetParName(3,"alpha");
  F->SetParName(4,"N");
  F->SetParName(5,"tail_frac");
  F->SetParName(6,"tail_lamb");

  F->SetParLimits(2,0.,1.);
  F->SetParLimits(5,0.,0.5);
  F->SetParLimits(6,0.,2.);
  F->SetLineWidth(1);
 
}

//-----------------------------------------------------------------------------
void fit_crystal_ball(const char* File, const char* Module, const char* Hist, 
		      double X0, double XMin, double XMax) {

  TH1F* h = (TH1F*) gh1(File,Module,Hist)->Clone("h_fit");
  double anorm = h->GetEntries()/10;

  TF1*  f;
  create_fit_function(f,X0,XMin,XMax);
  //  cb_init_parameters (f,anorm,X0,0.15,3.,1.,0.01,0.5);
  cb_init_parameters (f,anorm,X0,0.140,3,2.,0.05,1.);
  cb_fit             (h,f,XMin,XMax);
}

//-----------------------------------------------------------------------------
void fit_crystal_ball(TH1* Hist, double X0, double XMin, double XMax) {

  TH1F* h = (TH1F*) Hist;
  double anorm = h->GetEntries()/10;

  TF1*  f;
  create_fit_function(f,X0,XMin,XMax);
  cb_init_parameters (f,anorm,X0,0.120,2.,1.,0.01,0.5);
  cb_fit             (h,f,XMin,XMax);
}


//-----------------------------------------------------------------------------
// fit a histogram, normalized to a unit area
//-----------------------------------------------------------------------------
void fit_crystal_ball_norm(const char* File, const char* Module, const char* Hist, 
			   double X0, double XMin, double XMax) {

  TH1F* h = (TH1F*) gh1(File,Module,Hist)->Clone("h_fit");
  h->Scale(1./h->Integral());
  double anorm = h->Integral()/10;

  TF1*  f;
  create_fit_function(f,X0,XMin,XMax);
  cb_init_parameters (f,anorm,X0,0.15,3.,1.,0.01,0.1);
  cb_fit             (h,f,XMin,XMax);
}

//-----------------------------------------------------------------------------
// fit a histogram, normalized to a unit area
//-----------------------------------------------------------------------------
void fit_crystal_ball_norm(TH1* Hist, double X0, double XMin, double XMax) {

  TH1F* h = (TH1F*) Hist;
  h->Scale(1./h->Integral());
  double anorm = h->Integral()/10;

  TF1*  f;
  create_fit_function(f,X0,XMin,XMax);
  cb_init_parameters (f,anorm,X0,0.15,3.,1.,0.01,0.5);
  cb_fit             (h,f,XMin,XMax);
}


//-----------------------------------------------------------------------------
double cbe_f(double* X, double* P) {
  double x0, x1, f, sig, alpha, abs_alpha, delta, lam, n, a, b, dx0, dx1, xdx;

  x0        = P[1];
  sig       = P[2];
  alpha     = P[3];
  x1        = P[5];
  delta     = x1-x0;
  dx0       = X[0]-x0;
  xdx       = dx0/sig;
  lam       = sig*sig/delta; // *2.; //P[6];

  abs_alpha = fabs(alpha);
  n         = P[4];

  if (xdx > -abs_alpha) {
    dx1       = X[0]-x1;
    if (dx1 < 0) { 
					// gauss
      f = P[0]*TMath::Exp(-xdx*xdx/2.);
    }
    else {
      double d1 = delta/sig;
      f         = P[0]*TMath::Exp(-d1*d1/2.)*TMath::Exp(-pow(dx1,1.0)/(lam*1.0));
      //      f = P[0]*TMath::Exp(-xdx*xdx/2.);
    }
  }
  else {
    a = TMath::Power((n/abs_alpha),n)*TMath::Exp(-(alpha*alpha)/2);
    b = n/abs_alpha-abs_alpha;
    f = P[0]*a*TMath::Power((b-xdx),-n);
  }
  
  return f;
}

//-----------------------------------------------------------------------------
void cbe_init_parameters(TF1* F, double P0, double P1, double P2, double P3, 
			 double P4, double P5,double P6) {
  F->SetParameter(0,P0);
  F->SetParameter(1,P1);
  F->SetParameter(2,P2);
  F->SetParameter(3,P3);
  F->SetParameter(4,P4);
  F->SetParameter(5,P5);

//   F->SetParameter(6,P6);
//   F->SetParError (6,0.);
}

//-----------------------------------------------------------------------------
void cbe_fit(TH1F* Hist, TF1* F, double XMin, double XMax) {
  Hist->GetXaxis()->SetRangeUser(XMin,XMax);
  Hist->Fit(F,"L","",XMin,XMax);
}

//-----------------------------------------------------------------------------
void cbe_create_fit_function(TF1*& F, double X0, double XMin, double XMax) {

  F = new TF1("cbe_f",cbe_f,XMin,XMax,6);

  F->SetParName(0,"anorm");
  F->SetParName(1,"mean");
  F->SetParName(2,"sigma");
  F->SetParName(3,"alpha");
  F->SetParName(4,"N");
  F->SetParName(5,"x1");
  //  F->SetParName(6,"scale");

  F->SetParLimits(2,0.,1.);
  F->SetParLimits(5,0.,5.);
  //  F->SetParLimits(6,0.,5.);

  F->SetLineWidth(1);
}

//-----------------------------------------------------------------------------
void cbe_fit_crystal_ball(const char* File, const char* Module, const char* Hist, 
			  double X0, double XMin, double XMax, double Sigma = 0.18) {

  TH1F* h = (TH1F*) gh1(File,Module,Hist)->Clone("h_fit");

  double anorm = h->GetEntries()/10;

  printf("anorm = %12.5e\n",anorm);

  cbe_create_fit_function(f,X0,XMin,XMax);

  cbe_init_parameters    (f,anorm,X0,Sigma,1.,3.,X0+0.6,0.5);

  h->Draw();
  //  f->Draw("same");

  cbe_fit                (h,f,XMin,XMax);
}

//-----------------------------------------------------------------------------
// initial parameter setting assumes fit of delta(p) = P(reco)-P(MC) distribution
//-----------------------------------------------------------------------------
void cbe_fit_crystal_ball(TH1* Hist, double X0, double XMin, double XMax) {

  TF1*  f;
  TH1F* h      = (TH1F*) Hist;
  double anorm = h->GetEntries()/10;


  cbe_create_fit_function(f,X0,XMin,XMax);
  cbe_init_parameters    (f,anorm,X0,0.180,2.,1.,X0+0.6,0.5);
  cbe_fit                (h,f,XMin,XMax);
}


