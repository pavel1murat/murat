//



//-----------------------------------------------------------------------------
double rad_corr_jpsi(double* X, double* P) {

  const double alpha(1./137066);
  const double me(0.511), mmu(105.6), me2(me*me);
  
  double beta, f;

  double x   = X[0];
  double q2  = P[0];
  
  beta = alpha/TMath::Pi()*(TMath::Log(q2/me2)-1);

  f    = beta*TMath::Power(1-x,beta-1);

  return f;
}



//-----------------------------------------------------------------------------
void plot_rad_corr() {

  double mmu(105.6);
  
  TF1* f = new TF1("f",rad_corr_jpsi,0,1,1);

  f->SetParameter(0,mmu*mmu);


  TCanvas*  c = new TCanvas("c","c",0,0,1000,700);

  f->Draw();
}
