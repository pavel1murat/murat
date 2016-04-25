///////////////////////////////////////////////////////////////////////////////
// fit efficiency vs straw hit emax
// data come from e31s5730 dataset, each fileset of which is 100K events 
// generated with a given EMAX threshold (000002 = emax = 0.002 MeV)
///////////////////////////////////////////////////////////////////////////////

//-----------------------------------------------------------------------------
double f(double* X, double* P) {
  double x = X[0];

  double fun = P[0]/(1+(P[2]*P[2])*(x-P[1])*(x-P[1])/(x));

  return fun;
}


void fit_eff_vs_emax() {

  int const npt(8);

  double emax [npt] = { 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009 };

  double eff  [npt] = { 0.10330, 0.11036, 0.11265, 0.11200, 0.11246, 0.11075, 0.10913, 0.10855};

  double ex[npt], ey[npt];


  for (int i=0; i<npt; i++) {
    ex[i] = 0;
    ey[i] = eff[i]/100.;
  }

  TGraphErrors* gr = new TGraphErrors(npt,emax,eff,ex,ey);

  gr->Draw();

  TF1* f = new TF1("f",f,0.001,0.01,3);

  f->SetParameters(0.11,0.005,0.1);

  gr->Fit(f,"","",0.002,0.009);

  //  f->Draw("APL");
  //  f->Draw("same");
}
