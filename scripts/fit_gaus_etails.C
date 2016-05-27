// follow Souvik's idea 


namespace {
  int _Debug = 0;
}

//-----------------------------------------------------------------------------
double gauss_etails(double* X, double* P) {

  double f(0);
  double sig;

  double x  = X[0];
  double dx = x-P[1];

  double sig = fabs(P[2]);
  double dxs = dx/sig;

  if (_Debug != 0) {
    printf("x , dx, sig, dxs, P[3], P[4] =  : %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",
	   x , dx, sig, dxs,P[3],P[4]);
  }

  if (dx < -P[3]) {
    double d1 = P[3]/fabs(P[2]);             // supposed to be always positive
    f = P[0]*TMath::Gaus(-d1*d1/2)*TMath::Exp((-(P[1]-P[3]-x)/(sig*sig/P[3])));

    if (_Debug != 0) {
      printf(" branch 1: d1 = %12.5e, f = %12.5e\n", d1,f);
    }

  }
  else if (dx < P[4]) {
    f = P[0]*TMath::Gaus(-dxs*dxs/2);

    if (_Debug != 0) {
      printf(" branch 2: f = %12.5e\n", d1,f);
    }

  }
  else {
    double d2 = P[4]/fabs(P[2]);             // supposed to be always positive
    f = P[0]*TMath::Gaus(-d2*d2/2)*TMath::Exp((-(x-P[1]-P[4])/(sig*sig/P[4])));

    if (_Debug != 0) {
      printf(" branch 3: d2 = %12.5e, f = %12.5e\n", d2,f);
    }

  }

  return f;
}

//-----------------------------------------------------------------------------
int fit_gauss_etails(TH1* Hist, double A, double X0, double Sig, double D1, double D2, double XMin, double XMax) {

  TF1* fget = new TF1("fget",gauss_etails,XMin,XMax,5);

  fget->SetParameters(A,X0,Sig,D1,D2);

  fget->SetParLimits(2,0,100*Sig);
  fget->SetParLimits(3,0,100*D1 );
  fget->SetParLimits(4,0,100*D2 );

  Hist->Fit(fget,"","",XMin, XMax);

  return 0;
}

//-----------------------------------------------------------------------------
int fit_gauss_etails(const char* Fn, const char* Module, const char* Hist, 
		     double A, double X0, double Sig, double D1, double D2, double XMin, double XMax) {

  TH1F* h = (TH1F*) gh1(Fn,Module,Hist)->Clone("h_fit_gauss_etails");

  fit_gauss_etails(h,A,X0,Sig,D1,D2,XMin,XMax);

  return 0;
}


//-----------------------------------------------------------------------------
void test_gauss_etails() {

  _Debug = 1;

  TF1* fget = new TF1("fget",gauss_etails,100,106,5);

  fget->SetParameters(10000.,104.3,0.2,2,2);

  double x = 104.3;
  double f = fget->Eval(x);

  printf("fget(%12.5e) = %12.5e\n",x,f);
}
