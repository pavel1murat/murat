//


//-----------------------------------------------------------------------------
double kuraev_fadin(double* X, double* P) {

  double alpha = 1./137.036;
  double pi    = TMath::Pi();
  double q2    = 105.6*105.6;
  double m2    = 0.511*0.511;
  double x     = X[0];
  double beta  = (alpha/pi)*(log(q2/m2)-1.);
  double f     = beta*TMath::Power(1-x,beta-1);

  return f;
}


//-----------------------------------------------------------------------------
void plot_kuraev_fadin() {

  TF1* f_kf = new TF1("f_kf",kuraev_fadin,0,1,0);

  f_kf->SetNpx(10000);
  f_kf->Draw();
}
