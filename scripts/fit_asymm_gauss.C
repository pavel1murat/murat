///////////////////////////////////////////////////////////////////////////////
// fit distribution with a pseudo-gaussian function which has different values of
// gaussian sigma to the lft and to the right of maximum
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
double asymm_gauss(double* X, double* P) {

  double sig;

  double x  = X[0];
  double dx = x-P[1];

  if (dx <= 0) sig = P[2];
  else         sig = P[3];

  double f = P[0]*TMath::Exp(-(dx*dx)/(2*sig*sig));

  return f;
}

//-----------------------------------------------------------------------------
int fit_asymm_gauss(TH1* Hist, double A, double X0, double SigN, double SigP, double XMin, double XMax) {

  TF1* fag = new TF1("fag",asymm_gauss,-2*XMin,2*XMax,4);

  fag->SetParameters(A,X0,SigN,SigP);

  Hist->Fit(fag,"","",XMin, XMax);

  return 0;
}

//-----------------------------------------------------------------------------
int fit_asymm_gauss(const char* Fn, const char* Module, const char* Hist, 
		   double A, double X0, double SigLeft, double SigRight, double XMin, double XMax) {

  TH1F* h = (TH1F*) gh1(Fn,Module,Hist)->Clone("h_fit_asymm_gauss");

  fit_asymm_gauss(h,A,X0,SigLeft,SigRight,XMin,XMax);

  return 0;
}
