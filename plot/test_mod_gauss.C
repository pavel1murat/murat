//-----------------------------------------------------------------------------
// ROOT::Math::lognormal_pdf(double x, double  m, double s, double x0 = 0)		
//-----------------------------------------------------------------------------
double f_mod_gaus(double* X, double* P) {
  double x   = X[0];

  double mu  = P[0];
  double sig = P[1];

  double x1 = mu*(x + 10./x);

  double d = (x*x+1/(x*x)-2)/(2*sig*sig);

  double ff = exp(-d);

  printf("x, x1, mu, sig, d, ff : %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",x, x1, mu,sig,d, ff);
  return ff;
}


//-----------------------------------------------------------------------------
// parameterization of the pbar background - uncertainty 100%...
// sig^2/<x>^2 = exp(sigma^2)-1
// <x>         = exp(mu+sigma^2/2)
//-----------------------------------------------------------------------------
void test_mod_gauss(double Mean = 1.e-2, double SigmaOverMean = 1.) {
  
  //  double ROOT::Math::lognormal_pdf	(double x, double m, double s, double 	x0 = 0)

  double mu  = Mean;
  double sig = SigmaOverMean*mu;

  TF1* f = new TF1("f_mod_gaus", f_mod_gaus,0.,10,2);
  f->SetNpx(1000);
  f->SetParameters(mu,sig);
  f->Draw();
}
