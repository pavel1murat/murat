///////////////////////////////////////////////////////////////////////////////
// calculate, in 2 different ways, integrals of the LO and NLO DIO spectra
// to compare their ratio with papers by Czarnecki et al
///////////////////////////////////////////////////////////////////////////////
namespace {
  double emu(105.194), mAl(25133.), mmu(105.6583715) ;
  double z(13.), alpha(1./137.036);
};

//-----------------------------------------------------------------------------
double dio(double E) {

  double a5(8.6434), a6(1.16874), a7(-1.87828e-2), a8(9.16327e-3);
  double de, de5, w;

  double eb = mmu*(z*alpha)*(z*alpha)/2.;
  double er = E*E/(2*mAl);

  double x = E;


  x5 = x*x*x*x*x;

  w   = 1.e-17*x5*(a5 + x*(a6+x*(a7+a8*x)));

  if (x < 0) w = 0;

  return w;
}


//-----------------------------------------------------------------------------
double dio_rc(double E) {

  double x = E/mmu;

  double x5 = x*x*x*x*x;
  
  double f = 1.e-4*(1.44*pow(x,0.023)-0.22)*x5/mmu;

  return f;
}



//-----------------------------------------------------------------------------
double f_dio(double* X, double* P) {
  return dio(X[0]);
}

//-----------------------------------------------------------------------------
double f_dio_rc(double* X, double* P) {
  return dio_rc(X[0]);
}


//-----------------------------------------------------------------------------
void integrate_dio_rc(double dmax = 1.) {

  TF1* f0 = new TF1("f0",f_dio,0,5,0);

  double qint0 = f0->Integral(0,dmax);

  printf("qint0() = %12.5e\n",qint0);

  TF1* f = new TF1("f1",f_dio_rc,0,5,0);

  double qint = f1->Integral(0,dmax);

  printf("qint(rc) = %12.5e\n",qint);

  printf("ratio = %12.5e\n",qint0/qint);

}
