//-----------------------------------------------------------------------------
// ROOT::Math::lognormal_pdf(double x, double  m, double s, double x0 = 0)		
//-----------------------------------------------------------------------------
double f_f(double* X, double* P) {
  double x   = X[0];
  double mu  = P[0];
  double sig = P[1];
  double x0  = P[2];

  return  ROOT::Math::lognormal_pdf(x,mu,sig,x0);
}


//-----------------------------------------------------------------------------
void test_lognormal_001(double Mu = 1, double Sigma = 0.1, double X0 = 0) {
  
  TF1* f = new TF1("f1",f_f,-25,25,3);
  f->SetParameters(Mu,Sigma,X0);
  f->SetNpx(50000);

  f->Draw();
}


//-----------------------------------------------------------------------------
// parameterization of the pbar background - uncertainty 100%...
// sig^2/<x>^2 = exp(sigma^2)-1
// <x>         = exp(mu+sigma^2/2)
//-----------------------------------------------------------------------------
void test_lognormal_pbar(double Mean = 1.e-2, double SigmaOverMean = 1., double X0 = 0) {
  
  //  double ROOT::Math::lognormal_pdf	(double x, double m, double s, double 	x0 = 0)

  double sig = sqrt(log(1+SigmaOverMean));
  double mu  = log(Mean)-sig*sig/2;

  TF1* f = new TF1("f_lognormal", f_f,-0.1,0.1,3);
  f->SetNpx(10000);
  f->SetParameters(mu,sig,X0);
  f->Draw();
}


//-----------------------------------------------------------------------------
// parameterization of the pbar background - uncertainty 100%...
// sig^2/<x>^2 = exp(sigma^2)-1
// <x>         = exp(mu+sigma^2/2)
// x_MP        = exp(mu-sigma^2);
// sigma^2     = exp(2mu+sigma^2)(exp(sigma^2)-1)
// sigma/x_MP  = exp(3*sigma^2/2)*sqrt(exp(sigma^2)-1)
//-----------------------------------------------------------------------------
void test_lognormal_pbar(double MP = 1.e-2, double SigmaOverMP = 1., double X0 = 0) {
  
  //  double ROOT::Math::lognormal_pdf	(double x, double m, double s, double 	x0 = 0)

  double exp_sig2 = 1/2+sqrt(SigmaOverMP*SigmaOverMP+1/4);
  doubel sig      = sqrt(log(exp_sig2));
  double mu       = log(MP)-sig*sig/2;

  TF1* f = new TF1("f_lognormal", f_f,-0.1,0.1,3);
  f->SetNpx(10000);
  f->SetParameters(mu,sig,X0);
  f->Draw();
}
